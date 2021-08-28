HZZ Analyzer for CMS Run2

------

To install:

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_6_12

cd CMSSW_10_6_12/src

cmsenv

git cms-init

git clone -b 106X_2l2q https://github.com/jialin-guo1/UFHZZAnalysisRun2.git

cp UFHZZAnalysisRun2/install*.sh .

./install_2.sh

cp UFHZZAnalysisRun2/Utilities/crab/* .

voms-proxy-init --valid=168:00
#probably need "voms-proxy-init -voms cms -rfc"

source /cvmfs/cms.cern.ch/crab3/crab.sh

test for cmsRun
note: change path of JER and QGTag database files to fit interactively run in cfg.py
i.e: vim UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_106X_Legacy18_2l_cfg.py
     comment out QGdBFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/QGL_"+qgDatabaseVersion+".db" and dBJERFile = "src/UFHZZAnalysisRun2/UFHZZ4LAna/data/Autumn18_V7_MC.db"
cmsRun UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_106X_Legacy18_2l_cfg.py

for crab job 
Data:
python SubmitCrabJobs.py -t "myTask_Data" -d 2018Data.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateData_106X_Legacy18_2l_cfg.py

or similary for MC:
python SubmitCrabJobs.py -t "myTask_MC" -d 2018MC.txt -c UFHZZAnalysisRun2/UFHZZ4LAna/python/templateMC_106X_Legacy18_2l_cfg.py

You can use manageCrabTask.py to check the status, resubmit, or kill your task. E.g. after submitting:

nohup python -u manageCrabTask.py -t resultsAna_myTask_Data -r -l >& managedata.log &

This will start an infinite loop of running crab resubmit on all of your tasks, then sleep for 30min. You should kill the process once all of your tasks are done. Once all of your tasks are done, you should run the following command to purge your crab cache so that it doesn't fill up:

python manageCrabTask.py -t resultsAna_myTask_Data -p

UFHZZ4LAna/python/templateMC_102X_Legacy16_4l_cfg.py
UFHZZ4LAna/python/templateMC_102X_Legacy17_4l_cfg.py
UFHZZ4LAna/python/templateMC_102X_Legacy18_4l_cfg.py
UFHZZ4LAna/python/templateData_102X_Legacy16_3l_cfg.py
UFHZZ4LAna/python/templateData_102X_Legacy17_3l_cfg.py
UFHZZ4LAna/python/templateData_102X_Legacy18_3l_cfg.py
