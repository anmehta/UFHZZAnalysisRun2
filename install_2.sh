SCRAM_ARCH=slc7_amd64_gcc700; export SCRAM_ARCH
cmsrel CMSSW_10_6_26
cd CMSSW_10_6_26/src/
cmsenv
git cms-init

git clone -b 106X_2l2q git@github.com:jialin-guo1/UFHZZAnalysisRun2.git

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

voms-proxy-init --rfc --voms cms

