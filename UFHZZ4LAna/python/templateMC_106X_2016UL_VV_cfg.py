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
process.GlobalTag.globaltag='106X_mcRun2_asymptotic_v17'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring(#'/store/mc/RunIISummer20UL16MiniAOD/GluGluHToZZTo2L2Q_M1500_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v13-v1/2540000/35301A1B-B3A0-6647-AC7F-768E31BBECEB.root',
                                    #'/store/mc/RunIISummer20UL16MiniAOD/VBF_HToZZTo2L2Q_M1500_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v13-v1/2560000/64A7F033-6B5F-884A-A864-ADAD69DC2779.root',
                                    '/store/mc/RunIISummer20UL16MiniAODv2/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/280000/7073E121-5BD4-264B-BF63-8FDBE14E5973.root',
                                    #DUMMYFILELIST
                                    )#(DUMMYFILELIST)

process.source = cms.Source("PoolSource",fileNames = myfilelist,
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )

process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string("WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17-v1_32.root")
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
                                         year = cms.untracked.int32(2016)##2016
                                         #roccor.Run2.v5/RoccoR2016aUL.txt pre  VFP -2016
                                         #roccor.Run2.v5/RoccoR2016bUL.txt post VFP 2016
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
                       era='2016postVFP-UL')

process.load("RecoEgamma.EgammaTools.calibratedEgammas_cff")
process.calibratedPatElectrons.correctionFile = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2016_UltraLegacy_postVFP_RunFineEtaR9Gain_v1"
process.calibratedPatElectrons.src = cms.InputTag("slimmedElectrons")

# FSR Photons
process.load('UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff')

import os
# Jet Energy Corrections
from CondCore.DBCommon.CondDBSetup_cfi import *
#era = "Summer19UL16_V7_MC"
# for HPC
#dBFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/"+era+".db"
# for crab
#dBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/"+era+".db"
#process.jec = cms.ESSource("PoolDBESSource",
#                          CondDBSetup,
#                          connect = cms.string("sqlite_file:"+dBFile),
#                          toGet =  cms.VPSet(
#        cms.PSet(
#            record = cms.string("JetCorrectionsRecord"),
#            tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
#            label= cms.untracked.string("AK4PF")
#            ),
#        cms.PSet(
#            record = cms.string("JetCorrectionsRecord"),
#            tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
#            label= cms.untracked.string("AK4PFPuppi")
#            ),

#        cms.PSet(
#            record = cms.string("JetCorrectionsRecord"),
#            tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK8PFchs"),
#            label= cms.untracked.string("AK8PFchs")
#            ),
#        )
#)
#process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

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

#process.AK8PFJetCorrFactors = process.updatedPatJetCorrFactors.clone(
#    src = cms.InputTag("slimmedJetsAK8"),
#    levels = ['L2Relative',
#              'L3Absolute',
#              'L2L3Residual'],
#    payload = 'AK8PFPuppi' )

process.slimmedJetsJEC = process.updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("jetCorrFactors"))
    )

#process.slimmedJetsAK8JEC = process.updatedPatJets.clone(
#    jetSource = cms.InputTag("slimmedJetsAK8"),
#    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("AK8PFJetCorrFactors"))
#    )

### add pileup id and discriminant to patJetsReapplyJEC
#process.load("RecoJets.JetProducers.PileupJetID_cfi")
#process.pileupJetIdUpdated = process.pileupJetId.clone(
#    jets=cms.InputTag("slimmedJets"),
#    inputIsCorrected=False,
#    applyJec=True,
#    vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
#)


from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL16     #(or _chsalgos_106X_UL16APV for APV samples)
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
        jets=cms.InputTag("slimmedJets"),
        #inputIsCorrected=True,
        #applyJec=False,
        inputIsCorrected=False,
        applyJec=True,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos = cms.VPSet(_chsalgos_106X_UL16),
    )

process.slimmedJetsJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
process.slimmedJetsJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']


# JER
#process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
# for hpc
#dBJERFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Summer20UL16_JRV3_MC.db"
## for crab
#dBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Summer20UL16_JRV3_MC.db"
#process.jer = cms.ESSource("PoolDBESSource",
#        CondDBSetup,
#        connect = cms.string("sqlite_file:"+dBJERFile),
#        toGet = cms.VPSet(
#            cms.PSet(
#                record = cms.string('JetResolutionRcd'),
#                tag    = cms.string('JR_Summer20UL16_JRV3_MC_PtResolution_AK4PFchs'),
#                label  = cms.untracked.string('AK4PFchs_pt')
#                ),
#            cms.PSet(
#                record = cms.string('JetResolutionRcd'),
#                tag    = cms.string('JR_Summer20UL16_JRV3_MC_PhiResolution_AK4PFchs'),
#                label  = cms.untracked.string('AK4PFchs_phi')
#                ),
#            cms.PSet(
#                record = cms.string('JetResolutionScaleFactorRcd'),
#                tag    = cms.string('JR_Summer20UL16_JRV3_MC_SF_AK4PFchs'),
#                label  = cms.untracked.string('AK4PFPuppi')
#                )
#            )
#        )
#process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')


#QGTag
process.load("CondCore.CondDB.CondDB_cfi")
qgDatabaseVersion = 'cmssw8020_v2'
# for hpc
#QGdBFile = os.environ.get('CMSSW_BASE')+"/src/UFHZZAnalysisRun2/UFHZZ4LAna/data/QGL_"+qgDatabaseVersion+".db"
# for crab
QGdBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/QGL_"+qgDatabaseVersion+".db"
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

# compute corrected pruned jet mass
#process.corrJets = cms.EDProducer ( "CorrJetsProducer",
#                                    jets    = cms.InputTag( "slimmedJetsAK8JEC" ),
#                                    vertex  = cms.InputTag( "offlineSlimmedPrimaryVertices" ),
#                                    rho     = cms.InputTag( "fixedGridRhoFastjetAll"   ),
#                                    payload = cms.string  ( "AK8PFchs" ),
#                                    isData  = cms.bool    (  False ),
#                                    year    = cms.untracked.int32(2016))


# Recompute MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
            isData=False,
            )

from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
        #TheJets = cms.InputTag("updatedPatJetsUpdatedJEC"), #this should be the slimmedJets collection with up to date JECs !
        TheJets = cms.InputTag("slimmedJetsJEC"), #this should be the slimmedJets collection with up to date JECs !
        DataEraECAL = cms.string("UL2016postVFP"),
        DataEraMuon = cms.string("2016postVFP"),
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
                              photonSrc    = cms.untracked.InputTag("slimmedPhotons"),
                              electronUnSSrc  = cms.untracked.InputTag("slimmedElectrons"),
                              electronSrc  = cms.untracked.InputTag("calibratedPatElectrons"),
                              muonSrc      = cms.untracked.InputTag("calibratedMuons"),
                              tauSrc      = cms.untracked.InputTag("slimmedTaus"),
                              jetSrc       = cms.untracked.InputTag("slimmedJetsJEC"),
                              #mergedjetSrc = cms.untracked.InputTag("corrJets"),
                              #mergedjetSrc = cms.untracked.InputTag("slimmedJetsAK8JEC"),
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
                              skimLooseLeptons = cms.untracked.int32(1),
                              skimTightLeptons = cms.untracked.int32(1),
                              bestCandMela = cms.untracked.bool(False),   # for differential measurements
                              #-2016 pre VFP; 2016 post VFP
                              year = cms.untracked.int32(2016),####for year put 2016,2017, or 2018 to select correct Muon training and electron MVA
                             )

process.p = cms.Path(process.fsrPhotonSequence*
                     process.boostedMuons*
                     process.calibratedMuons*
                     process.egmGsfElectronIDSequence*
                     process.egmPhotonIDSequence*
                     process.egammaPostRecoSeq*
                     process.calibratedPatElectrons*
                     process.jetCorrFactors*
                     process.pileupJetIdUpdated*
                     process.slimmedJetsJEC*
                     process.QGTagger*
                     #process.AK8PFJetCorrFactors*
                     #process.slimmedJetsAK8JEC*
                     process.fullPatMetSequence*
                     #process.corrJets*
                     process.mergedGenParticles*process.myGenerator*process.rivetProducerHTXS*#process.rivetProducerHZZFid*
                     process.prefiringweight*
                     process.Ana
                     )

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
task = getPatAlgosToolsTask(process)
process.p.associate(task)
