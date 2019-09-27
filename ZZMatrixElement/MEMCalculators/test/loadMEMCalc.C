

{

  // These libs are already linked to libZZMatrixElementMEMCalculators.so
//   gSystem->Load("$CMSSW_BASE/src/ZZMatrixElement/MELA/data/$SCRAM_ARCH/libmcfm.so");
//   gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libZZMatrixElementMELA.so");
//   gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libZZMatrixElementMEKD.so");
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libZZMatrixElementMEMCalculators.so");

  // No way to generate a dictionary for MEMCalculators.h. You will have to use a compiled macro to
  // use the MEMs class.
  //  gROOT->LoadMacro("$CMSSW_BASE/src/ZZMatrixElement/MEMCalculators/interface/MEMCalculators.h+");

}
