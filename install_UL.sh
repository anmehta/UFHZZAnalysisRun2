source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700

scramv1 project CMSSW CMSSW_10_6_26
cd CMSSW_10_6_26/src/
#cmsenv
eval `scramv1 runtime -sh`

git cms-init

git clone -b UL_10_6_26 https://github.com/qyguo/UFHZZAnalysisRun2.git

##git cms-merge-topic asculac:Electron_XGBoost_MVA_16UL_17UL

git cms-addpkg GeneratorInterface/RivetInterface

git cms-addpkg SimDataFormats/HTXS

git cms-addpkg RecoEgamma/PhotonIdentification

git cms-addpkg RecoEgamma/ElectronIdentification

git cms-merge-topic cms-egamma:EgammaPostRecoTools

git cms-addpkg RecoEgamma/EgammaTools

git clone https://github.com/cms-egamma/EgammaPostRecoTools.git

mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/.

#git clone https://github.com/cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/

git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/

git cms-addpkg EgammaAnalysis/ElectronTools

git cms-addpkg  RecoJets/JetProducers

git cms-addpkg PhysicsTools/PatAlgos/

git clone -b v2.3.5 https://github.com/JHUGen/JHUGenMELA

sh JHUGenMELA/MELA/setup.sh -j 8

git clone https://github.com/bachtis/Analysis.git -b KaMuCa_V4 KaMuCa

git clone -b tmp_Ferrico https://github.com/ferrico/KinZfitter.git

scramv1 b -j 8

# local run and tests

# cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_106X_2018UL_cfg.py

# cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/Sync_106X_2017UL_cfg.py

# pre-existing

cp UFHZZAnalysisRun2/Utilities/crab/* .

source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init --rfc --voms cms

#python SubmitCrabJobs.py -t "DataUL18" -d SampleList_UL18_Data.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_106X_2018UL_cfg.py

#####
# cfg of Data and MC to be used for each year
#####

#UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_106X_2016ULAPV_cfg.py
#UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_106X_2016UL_cfg.py
#UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_106X_2017UL_cfg.py
#UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_106X_2018UL_cfg.py
#UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_106X_2016ULAPV_cfg.py
#UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_106X_2016UL_cfg.py
#UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_106X_2017UL_cfg.py
#UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_106X_2018UL_cfg.py

#####
# or similary for MC:
#####

#python SubmitCrabJobs.py -t "MC_UL18" -d SampleList_UL18_MC.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_106X_2018UL_cfg.py

#You can use manageCrabTask.py to check the status, resubmit, or kill your task. E.g. after submitting:

#nohup python -u manageCrabTask.py -t resultsAna_Data_M17_Feb19 -r -l >& managedata.log &

#This will start an infinite loop of running crab resubmit on all of your tasks, then sleep for 30min. You should kill the process once all of your tasks are done. Once all of your tasks are done, you should run the following command to purge your crab cache so that it doesn't fill up:

#python manageCrabTask.py -t resultsAna_Data_M17_Feb19 -p
