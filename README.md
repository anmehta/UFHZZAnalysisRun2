HZZ Analyzer for CMS Run2 UL

------

To install:

SCRAM_ARCH=slc7_amd64_gcc700; export SCRAM_ARCH

cmsrel CMSSW_10_6_20

cd CMSSW_10_6_20/src/

cmsenv

git cms-init

git clone -b 10_6_20_UL https://github.com/tjavaid/UFHZZAnalysisRun2.git

git cms-merge-topic asculac:Electron_XGBoost_MVA_16UL_17UL

git cms-addpkg GeneratorInterface/RivetInterface

git cms-addpkg SimDataFormats/HTXS

git cms-addpkg RecoEgamma/EgammaTools

git cms-addpkg RecoEgamma/PhotonIdentification

git cms-addpkg RecoEgamma/ElectronIdentification

git cms-merge-topic cms-egamma:EgammaPostRecoTools

mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.

git clone https://github.com/cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/

git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/

git cms-addpkg EgammaAnalysis/ElectronTools

git cms-addpkg  RecoJets/JetProducers

git cms-addpkg PhysicsTools/PatAlgos/

git clone -b v2.3.5 https://github.com/JHUGen/JHUGenMELA

sh JHUGenMELA/MELA/setup.sh -j 8

scramv1 b -j 8

voms-proxy-init --rfc --voms cms

cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_106X_2018UL_cfg.py

cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_106X_2017UL_cfg.py

==============================================================

# pre-existing

#cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_102X_2017_Legacy_cfg.py

#cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_102X_2016_Legacy_cfg.py

cp UFHZZAnalysisRun2/Utilities/crab/* .

voms-proxy-init --valid=168:00
#probably need "voms-proxy-init -voms cms -rfc"

source /cvmfs/cms.cern.ch/crab3/crab.sh

python SubmitCrabJobs.py -t "myTask_Data" -d datasets_2016ReReco.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_80X_M1703Feb_2l_cfg.py

or similary for MC:

python SubmitCrabJobs.py -t "myTask_MC" -d datasets_Summer16_25ns_MiniAOD.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_80X_M17_4l_cfg.py

You can use manageCrabTask.py to check the status, resubmit, or kill your task. E.g. after submitting:

nohup python -u manageCrabTask.py -t resultsAna_Data_M17_Feb19 -r -l >& managedata.log &

This will start an infinite loop of running crab resubmit on all of your tasks, then sleep for 30min. You should kill the process once all of your tasks are done. Once all of your tasks are done, you should run the following command to purge your crab cache so that it doesn't fill up:

python manageCrabTask.py -t resultsAna_Data_M17_Feb19 -p

UFHZZ4LAna/python/templateMC_102X_Legacy16_4l_cfg.py
UFHZZ4LAna/python/templateMC_102X_Legacy17_4l_cfg.py
UFHZZ4LAna/python/templateMC_102X_Legacy18_4l_cfg.py
UFHZZ4LAna/python/templateData_102X_Legacy16_3l_cfg.py
UFHZZ4LAna/python/templateData_102X_Legacy17_3l_cfg.py
UFHZZ4LAna/python/templateData_102X_Legacy18_3l_cfg.py
