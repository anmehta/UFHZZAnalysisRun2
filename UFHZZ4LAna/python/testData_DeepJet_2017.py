import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("JetAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('UFHZZ4LAna')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.GlobalTag.globaltag='94X_dataRun2_v11'

process.Timing = cms.Service("Timing",
                             summaryOnly = cms.untracked.bool(True)
                             )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

myfilelist = cms.untracked.vstring(
#'/store/data/Run2017E/SingleElectron/MINIAOD/31Mar2018-v1/80000/A0D6B0A3-5937-E811-AC52-0CC47AA53D86.root',
#'/store/data/Run2017E/SingleElectron/MINIAOD/31Mar2018-v1/90000/06D3E100-6E37-E811-A4B0-0CC47AA53D5A.root',
#'/store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/90000/FEC62083-1E39-E811-B2A1-0CC47A4D75F8.root'
#'/store/data/Run2017B/DoubleEG/MINIAOD/31Mar2018-v1/00000/000A6D14-8037-E811-A09B-0CC47A5FBDC1.root',
'/store/data/Run2017B/DoubleEG/MINIAOD/31Mar2018-v1/00000/322DEFC6-8737-E811-9A74-001E67E6F7F1.root '

        #DUMMYFILELIST
        )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DeepJetTest.root")
)

#================================================================================
#=========================DeepAK8 Jets===========================================
#================================================================================
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from RecoBTag.MXNet.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsAll, _pfDeepBoostedJetTagsProbs, _pfDeepBoostedJetTagsMetaDiscrs, _pfMassDecorrelatedDeepBoostedJetTagsProbs, _pfMassDecorrelatedDeepBoostedJetTagsMetaDiscrs
updateJetCollection(
      process,
      jetSource = cms.InputTag('slimmedJetsAK8'),
      pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
      svSource = cms.InputTag('slimmedSecondaryVertices'),
      rParam = 0.8,
      jetCorrections = ('AK8PFPuppi', cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
      btagDiscriminators = _pfDeepBoostedJetTagsAll,
      postfix='AK8WithDeepTags',
      printWarning = True
      )
#process.updatedPatJetsTransientCorrected.addTagInfos = cms.bool(True)
#process.updatedPatJetsTransientCorrected.addBTagInfo = cms.bool(True)
DeepAK8 = cms.InputTag('selectedUpdatedPatJetsAK8WithDeepTags')


process.JetAnalysis(
          deepjet      = cms.untracked.InputTag("DeepAK8"),
)

process.p = cms.Path(
          process.JetAnalysis,
)
