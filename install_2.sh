git cms-init
git cms-addpkg GeneratorInterface/RivetInterface
git cms-addpkg SimDataFormats/HTXS
git cms-addpkg RecoEgamma/EgammaTools
git cms-addpkg RecoEgamma/PhotonIdentification
git cms-addpkg RecoEgamma/ElectronIdentification
git cms-merge-topic cms-egamma:EgammaPostRecoTools
git cms-addpkg PhysicsTools/PatAlgos/
#(vi .git/info/sparse-checkout & git read-tree -mu HEAD)
#git clone https://github.com/cms-analysis/HiggsAnalysis-ZZMatrixElement.git ZZMatrixElement
#cd ZZMatrixElement
#git checkout -b from-v223 v2.2.3
#. setup.sh -j 12
#cd ../
git clone -b v2.3.5 https://github.com/JHUGen/JHUGenMELA
sh JHUGenMELA/MELA/setup.sh -j 8
git clone -b RunII_v2 https://github.com/VBF-HZZ/UFHZZAnalysisRun2-Accessary.git
#mv UFHZZAnalysisRun2-Accessary/* ./
#rm -rf UFHZZAnalysisRun2-Accessary
mv UFHZZAnalysisRun2-Accessary/RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Autumn18_ID_ISO_cff.py ./RecoEgamma/ElectronIdentification/python/Identification/
mv UFHZZAnalysisRun2-Accessary/RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Summer16_ID_ISO_cff.py  ./RecoEgamma/ElectronIdentification/python/Identification/
mv UFHZZAnalysisRun2-Accessary/RecoEgamma/ElectronIdentification/python/ElectronMVAValueMapProducer_cfi.py  RecoEgamma/ElectronIdentification/python/
mv UFHZZAnalysisRun2-Accessary/RecoEgamma/ElectronIdentification/python/ElectronIDValueMapProducer_cfi.py RecoEgamma/ElectronIdentification/python/
mv UFHZZAnalysisRun2-Accessary/RecoEgamma/ElectronIdentification/data/MVAWeightFiles RecoEgamma/ElectronIdentification/data/
mv UFHZZAnalysisRun2-Accessary/* ./
rm -rf UFHZZAnalysisRun2-Accessary
#mv UFHZZAnalysisRun2-Accessary ../
git clone -b tmp_Ferrico https://github.com/ferrico/KinZfitter.git
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_102X_v3
scram b -j 8
