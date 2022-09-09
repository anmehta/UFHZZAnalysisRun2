import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("UFHZZ4LAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.GlobalTag.globaltag='106X_mcRun2_asymptotic_preVFP_v11'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

myfilelist = cms.untracked.vstring(
                                #DUMMYFILELIST,
                                #'/store/mc/RunIISummer20UL16MiniAODAPVv2/DYJetsToLL_Pt-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/2540000/08DAC15C-CF12-D441-BEE3-272A66A632A6.root',
                                #'/store/mc/RunIISummer20UL16MiniAODAPVv2/DYJetsToLL_Pt-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/120000/01D23322-7F98-5645-A05D-09B262AE82C5.root ',
                                '/store/mc/RunIISummer20UL16MiniAODAPVv2/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/2520000/0190A0DE-0D28-A644-AE94-D28B790F02C2.root',
                                )

process.source = cms.Source("PoolSource",fileNames = myfilelist,
                                    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                                                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DUMMYFILENAME.root")
)

# clean muons by segments 
process.boostedMuons = cms.EDProducer("PATMuonCleanerBySegments",
                                      src = cms.InputTag("slimmedMuons"),
                                      preselection = cms.string("track.isNonnull"),
                                      passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                      fractionOfSharedSegments = cms.double(0.499),
                     )


# Kalman Muon Calibrations
process.calibratedMuons = cms.EDProducer("KalmanMuonCalibrationsProducer",
                                         muonsCollection = cms.InputTag("boostedMuons"),
                                         isMC = cms.bool(True),
                                         isSync = cms.bool(False),
                                         useRochester = cms.untracked.bool(True),
                                         year = cms.untracked.int32(-2016)##2016 
                                         #roccor.Run2.v5/RoccoR2016aUL.txt pre  VFP -2016
                                         #roccor.Run2.v5/RoccoR2016bUL.txt post VFP  2016
                                         #year = cms.untracked.int32(2016)
                                         )

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
                                         src = cms.InputTag("slimmedElectrons"),
                                         cut = cms.string("pt > 5 && abs(eta)<2.5")
                                         )

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(123456), # for crab
        engineName = cms.untracked.string('TRandom3')
    )
)

####### new added
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,
                       runVID=True,
                       eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16UL_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                       era='2016preVFP-UL')

process.load("RecoEgamma.EgammaTools.calibratedEgammas_cff")
process.calibratedPatElectrons.correctionFile = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2016_UltraLegacy_preVFP_RunFineEtaR9Gain_v3"
process.calibratedPatElectrons.src = cms.InputTag("slimmedElectrons")

# FSR Photons
process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

import os
# Jet Energy Corrections
from CondCore.DBCommon.CondDBSetup_cfi import *
#=========================particleNet implatement======================
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from RecoBTag.ONNXRuntime.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsAll
from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll
updateJetCollection(
     process,
     jetSource = cms.InputTag('slimmedJetsAK8'),
     #jetSource = cms.InputTag('packedPatJetsAK8PFPuppi'+reclusterAK8JetPostFix+'SoftDrop'), #recluster
     pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
     svSource = cms.InputTag('slimmedSecondaryVertices'),
     rParam = 0.8,
     jetCorrections = ('AK8PFPuppi', cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
     btagDiscriminators = _pfDeepBoostedJetTagsAll+_pfParticleNetJetTagsAll,
     postfix='AK8WithDeepTags',
     #postfix='AK8CleanedWithZWithPuppiDaughters', #recluster
     printWarning = False
)
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
#newAK8Puppijet = "selectedUpdatedPatJetsAK8"+reclusterAK8JetPostFix+"WithPuppiDaughters"
newAK8Puppijet = "selectedUpdatedPatJetsAK8WithDeepTags"

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")

process.jetCorrFactors = process.updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 
              'L2Relative', 
              'L3Absolute'],
    payload = 'AK4PFchs' ) 


process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )



from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL16APV     #(or _chsalgos_106X_UL16APV for APV samples)
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone( 
        jets=cms.InputTag("slimmedJets"),
        #inputIsCorrected=True,
        #applyJec=False,
        inputIsCorrected=False,
        applyJec=True,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos = cms.VPSet(_chsalgos_106X_UL16APV),
    )

process.slimmedJetsJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
process.slimmedJetsJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']



#QGTag
process.load("CondCore.CondDB.CondDB_cfi")
qgDatabaseVersion = 'cmssw8020_v2'
# for hpc
QGdBFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/QGL_"+qgDatabaseVersion+".db"
# for crab
#QGdBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/QGL_"+qgDatabaseVersion+".db"
process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(messageLevel = cms.untracked.int32(1)),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PFchs'),
            label  = cms.untracked.string('QGL_AK4PFchs')
        ),
      ),
      connect = cms.string('sqlite_file:'+QGdBFile)
)
process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag( 'slimmedJetsJEC' )
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')
process.QGTagger.srcVertexCollection=cms.InputTag("offlinePrimaryVertices")


# Recompute MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
            isData=False,
            )

# Recomput Puppi MET
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True )
process.puppiNoLep.useExistingWeights     = True
process.puppi.useExistingWeights         = True
process.puppi.useExp                     = True
process.puppiNoLep.useExp            = True
runMetCorAndUncFromMiniAOD(process,
                           isData= False,
                           metType="Puppi",
                           postfix="Puppi",
                           jetFlavor="AK4PFPuppi",
                        )

from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
       # TheJets = cms.InputTag("updatedPatJetsUpdatedJEC"), #this should be the slimmedJets collection with up to date JECs !
        TheJets = cms.InputTag("slimmedJetsJEC"), #this should be the slimmedJets collection with up to date JECs !
        DataEraECAL = cms.string("UL2016preVFP"),
        DataEraMuon = cms.string("2016preVFP"),
        UseJetEMPt = cms.bool(False),
        PrefiringRateSystematicUnctyECAL = cms.double(0.2),
        PrefiringRateSystematicUnctyMuon = cms.double(0.2)
        )

# STXS
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)
process.myGenerator = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator"),
    signalParticlePdgIds = cms.vint32(25)
)
process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
  HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
  LHERunInfo = cms.InputTag('externalLHEProducer'),
  ProductionMode = cms.string('AUTO'),
)
# HZZ Fiducial from RIVET
process.rivetProducerHZZFid = cms.EDProducer('HZZRivetProducer',
  HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
)



# Analyzer
process.Ana = cms.EDAnalyzer('UFHZZ4LAna',
                              metpuppi     = cms.untracked.InputTag("slimmedMETsPuppi","","UFHZZ4LAnalysis"),
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronUnSSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              electronSrc  = cms.untracked.InputTag("calibratedPatElectrons"),
                              muonSrc      = cms.untracked.InputTag("calibratedMuons"),
                              tauSrc      = cms.untracked.InputTag("slimmedTaus"),
                              jetSrc       = cms.untracked.InputTag("slimmedJetsJEC"),
                              #mergedjetSrc = cms.untracked.InputTag("corrJets"),
                              mergedjetSrc = cms.untracked.InputTag(newAK8Puppijet),
                              metSrc       = cms.untracked.InputTag("slimmedMETs","","UFHZZ4LAnalysis"),
                              vertexSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
                              beamSpotSrc  = cms.untracked.InputTag("offlineBeamSpot"),
                              conversionSrc  = cms.untracked.InputTag("reducedEgamma","reducedConversions"),
                              isMC         = cms.untracked.bool(True),
                              isSignal     = cms.untracked.bool(True),
                              mH           = cms.untracked.double(125.0),
                              #CrossSection = cms.untracked.double(DUMMYCROSSSECTION),
                              CrossSection = cms.untracked.double(1),
                              FilterEff    = cms.untracked.double(1),
                              weightEvents = cms.untracked.bool(True),
                              elRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              muRhoSrc     = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                              rhoSrcSUS    = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
                              pileupSrc     = cms.untracked.InputTag("slimmedAddPileupInfo"),
                              pfCandsSrc   = cms.untracked.InputTag("packedPFCandidates"),
                              fsrPhotonsSrc = cms.untracked.InputTag("boostedFsrPhotons"),
                              prunedgenParticlesSrc = cms.untracked.InputTag("prunedGenParticles"),
                              packedgenParticlesSrc = cms.untracked.InputTag("packedGenParticles"),
                              genJetsSrc = cms.untracked.InputTag("slimmedGenJets"),
                              generatorSrc = cms.untracked.InputTag("generator"),
                              lheInfoSrc = cms.untracked.InputTag("externalLHEProducer"),
                              reweightForPU = cms.untracked.bool(True),
                              triggerSrc = cms.InputTag("TriggerResults","","HLT"),
                              triggerObjects = cms.InputTag("selectedPatTrigger"),
                              doJER = cms.untracked.bool(True),
                              doJEC = cms.untracked.bool(True),
                              doTriggerMatching = cms.untracked.bool(False),
                              #doTriggerMatching = cms.untracked.bool(True),
                              triggerList = cms.untracked.vstring(     
                                  'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v',
                                  'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
                                  'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v',
                                  'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v',
                                  'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
                                  'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
                                  'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
                                  'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                  'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
                                  'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
                                  'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v',
                                  'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v',
                                  'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v',
                                  'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
                                  'HLT_TripleMu_12_10_5_v',
                                  'HLT_Ele25_eta2p1_WPTight_Gsf_v',
                                  'HLT_Ele27_WPTight_Gsf_v',
                                  'HLT_Ele27_eta2p1_WPLoose_Gsf_v',
                                  'HLT_Ele32_eta2p1_WPTight_Gsf_v',
                                  'HLT_IsoMu20_v',
                                  'HLT_IsoTkMu20_v',
                                  'HLT_IsoMu22_v',
                                  'HLT_IsoTkMu22_v',
                                  'HLT_IsoMu24_v',
                                  'HLT_IsoTkMu24_v',
                              ),
                              verbose = cms.untracked.bool(False),              
                              skimLooseLeptons = cms.untracked.int32(2),              
                              skimTightLeptons = cms.untracked.int32(2),               
                              bestCandMela = cms.untracked.bool(False),   # for differential measurements
                              year = cms.untracked.int32(-2016),####for year put 2016,2017, or 2018 to select correct Muon training and electron MVA
                             )
process.metSequence = cms.Sequence(process.fullPatMetSequence+process.fullPatMetSequencePuppi)
process.p = cms.Path(process.fsrPhotonSequence*
                     process.boostedMuons*
                     process.calibratedMuons*
                     process.egmGsfElectronIDSequence*
                     process.egmPhotonIDSequence*
                     process.metSequence*
                     process.egammaPostRecoSeq*
                     process.calibratedPatElectrons*
                     process.jetCorrFactors*
                     process.pileupJetIdUpdated*
                     process.slimmedJetsJEC*
                     process.QGTagger*
                     #process.AK8PFJetCorrFactors*
                     #process.slimmedJetsAK8JEC*
                     #process.fullPatMetSequence*
                     #process.corrJets*
                     process.mergedGenParticles*process.myGenerator*process.rivetProducerHTXS*#process.rivetProducerHZZFid*
                     process.prefiringweight*
                     process.Ana
                     )

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
task = getPatAlgosToolsTask(process)
process.p.associate(task)
