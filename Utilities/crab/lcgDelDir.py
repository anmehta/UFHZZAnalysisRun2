import os, sys, fnmatch

def find_files(directory, pattern):
  for root, dirs, files in os.walk(directory):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = os.path.join(root, basename)
        yield filename

dirNames = [
'/cms/data/store/user/dsperka/UFHZZAnalysisRun2/Data_80X_2lskim_M17_Feb02/DoubleEG/crab_DoubleEG_Run2016G-23Sep2016-v1/170202_170810/'
]

for dirName in dirNames:
  if (not os.path.isdir(dirName)): continue
  print 'deleting',dirName
  for filename in find_files(dirName, '*'):
    #print 'Found C source:', filename
    os.system('lcg-del -b -l -D srmv2 srm://srm.ihepa.ufl.edu:8443/srm/v2/server?SFN='+str(filename))

  subdirs = [x[0] for x in os.walk(dirName)]

  for i in range(len(subdirs)):
    subdir = subdirs[len(subdirs)-1-i]
    print subdir
    os.system('lcg-del -b -l -d -D srmv2 srm://srm.ihepa.ufl.edu:8443/srm/v2/server?SFN='+subdir)
    
# to copy
#lcg-cp -v -b -D srmv2 file:/test   srm://srm.ihepa.ufl.edu:8443/srm/v2/server?SFN=/cms/data/store/user/dsperka/

