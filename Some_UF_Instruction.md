# Short recipe for UF computing 

[Web server for T2 Florida](http://tier2.ihepa.ufl.edu/)
[UF research computing](https://www.rc.ufl.edu)

--------------

#### A full list of ihepa machines
archer, alachua, melrose, newberry, or gainesville.ihepa.ufl.edu

```
ssh -X -Y ferrico@melrose.ihepa.ufl.edu
newberry.ihepa.ufl.edu
```

--------------

#### Copy files from UF Tier 2 storage
```
gfal-copy gsiftp://cmsio.rc.ufl.edu//cms/data/store/user/ferrico/<PATH_TO_FILE> file://<PATH_TO_COPY_TO>
gfal-copy gsiftp://cmsio.rc.ufl.edu//cms/data/store/user/t2/users/ferrico/Full_RunII/ggH/GluGluHToZZTo4L_M125_2017.root ./
tar xzvf heppyOutput_3.tgz

# Copy content in a folder to a existing folder on UF Tier 2
gfal-copy -r <PATH_TO_FOLDER> gsiftp://cmsio.rc.ufl.edu//cms/data/store/user/ferrico/<PATH_TO_FOLDER>
```
--------------

#### Interact with UF Tier 2 storage directly
```
uberftp cmsio.rc.ufl.edu "ls /cms/data/store/user/ferrico/"
uberftp cmsio.rc.ufl.edu "help" # To see available commands
```

#### To hadd files:
```
ssh -X -Y ferrico@hpg2.rc.ufl.edu
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_4_4/src/
cmsenv 
cd /cmsuf/data/store/user/t2/users/ferrico/
hadd -f new_path /cmsuf/data/store/user/ferrico/2018data/UFHZZAnalysisRun2/ #the second is the crab_path
```


#### Prefix path to access Tier2
```
    To access CMS data from IHEPA,
    please use root://cmsio2.rc.ufl.edu/cms/data/store/...
               gsiftp://cmsio.rc.ufl.edu/cms/data/store/...
```

#### Change password
In command line, do:
```
yppasswd
```
