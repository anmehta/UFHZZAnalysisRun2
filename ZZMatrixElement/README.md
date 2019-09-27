HiggsAnalysis-ZZMatrixElement
=============================

[Wiki page for the MEM packages (MEMCalculators wrapper and constituents)](https://twiki.cern.ch/twiki/bin/viewauth/CMS/HZZ4lME)
[Wiki page for the MELA package](https://twiki.cern.ch/twiki/bin/view/CMS/MELAProject)
[Wiki page for the MEKD package](https://twiki.cern.ch/twiki/bin/view/CMS/HZZ4lMEKD)


Checking out on top of other CMSSW packages
-------------------------------------------

Assuming that you have set up a CMSSW area with packages from the cms-sw/cmssw git repository (`git cms-init`), you have to set up this repository so that its .git folder does not interfere with the cms-sw/cmssw folder:

```
git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
```

At this point, git commands apply to this repository if issued inside the ZZMatrixElement folder or to the cmssw repository if issued outside (and in this case ZZMatrixElement is reported as an untracked folder by git status).


Mixing with packages from cvs
-------------------------------------------

Just set up the git part first, then checkout the cvs packages as usual. 
They will appear as untracked folders in git status.

You can move the ZZAnalysis folder into an existing area if needed, it will not break git providing that ZZAnalysis/.git/ is moved as well.


Committing and pushing
-------------------------------------------

```
git pull
[edit files]
git add [files to be added]
git commit -m ["commit message"] [files to be added]
git push origin master
```
