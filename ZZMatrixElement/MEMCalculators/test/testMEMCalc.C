#include <vector>
#include "TLorentzVector.h"
#include "../interface/MEMCalculators.h"
#include "TLorentzRotation.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

using namespace MEMNames;
using namespace std;

vector<TLorentzVector> Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi,double pt4l, double Y4l);

//<<<<<<< testMEMCalc.C
//void testMEMCalc(TString fileName,int channel){
//int main(){
//TString fileName="../../../../../../MELA/test/ggtoHtoZZ4l_VariousKDs_0+m0+h_withDisc_withMELA.root";
//int channel=3;
//=======
//void testMEMCalc(TString fileName,int channel){
int main(){
//TString fileName="/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_5_2_5/src/ZZMatrixElement/MELA/test/4e/HZZ4lTree_powheg15jhuGenV3H126_withDisc_new_withMELA_new_3.root";
//TString fileName="/afs/cern.ch/work/x/xiaomeng/test/myWorkingArea/CMSSW_5_2_5/src/ZZMatrixElement/MELA/test/2mu2e/HZZ4lTree_powheg15jhuGenV3H126_withDisc_new_withMELA_new_2.root";
TString fileName="root://lxcms00://data3/2014/HZZ_stat/140211/PRODFSR_8TeV/2mu2e/HZZ4lTree_powheg15jhuGenV3-0PMH125.6.root";
//TString fileName="/afs/cern.ch/work/c/chmartin/private/Fa2stuff/MC_from_JHU/SMHiggsToZZTo4L_M-126_8TeV_powheg15-JHUgenV4-pythia6_false.root";
int channel=3;
//>>>>>>> 1.8
  TChain* t = new TChain("SelectedTree");
  t->Add(fileName);
  
  fileName.ReplaceAll(".root","_test.root");

  //TFile* file = new TFile(fileName,"RECREATE");
  TFile* file = new TFile("test.root","RECREATE");

  TTree* tree = (TTree*) t->CloneTree(0,"fast");

  MEMs test(8,126,"", false);

  double p0plus_mela_test,p2_mela_test,p2qqb_mela_test;
  double bkg_mela_test,p0plus_melaNorm_test;
  double p0hplus_mela_test,p1plus_mela_test;
  double p0minus_mela_test,p1_mela_test;

  double p0plus_VAJHU_test,p2_VAJHU_test,p2qqb_VAJHU_test;
  double p0hplus_VAJHU_test,p1plus_VAJHU_test;
  double p0minus_VAJHU_test,p1_VAJHU_test;

  double qqZZ_VAMCFM_test,ggZZ_VAMCFM_test,p0plus_VAMCFM_test,ggzz_p0plus_VAMCFM_test ;
  double qqZZ_VAMCFMNorm_test;

  double p0plus_mad, p0minus_mad, p2_mad, qqZZ_mad;

  double KDsigbkg, KD0minus, KD0hplus, KD2mplus;

  double dummy1,dummy2;
	
	double pg1g4_mela, pg1g4_VAJHU, pg1g2_mela, pg1g2_VAJHU, pg1g4_pi2_VAJHU, pg1g2_pi2_VAJHU;
	double ggzz_VAMCFM, ggzz_c5_VAMCFM, ggzz_c1_VAMCFM,ggzz_ci_VAMCFM;
  // new branches
  tree->Branch("p0plus_mela_test",&p0plus_mela_test,"p0plus_mela_test/D");
  tree->Branch("p0hplus_mela_test",&p0hplus_mela_test,"p0hplus_mela_test/D");
  tree->Branch("p0minus_mela_test",&p0minus_mela_test,"p0minus_mela_test/D");
  tree->Branch("p1_mela_test",&p1_mela_test,"p1_mela_test/D");
  tree->Branch("p1plus_mela_test",&p1plus_mela_test,"p1plus_mela_test/D");
  tree->Branch("p2_mela_test",&p2_mela_test,"p2_mela_test/D");
  tree->Branch("p2qqb_mela_test",&p2qqb_mela_test,"p2qqb_mela_test/D");

  tree->Branch("p0plus_melaNorm_test",&p0plus_melaNorm_test,"p0plus_melaNorm_test/D");
  tree->Branch("bkg_mela_test",&bkg_mela_test,"bkg_mela_test/D");

  tree->Branch("p0plus_VAMCFM_test",&p0plus_VAMCFM_test,"p0plus_VAMCFM_test/D");
  tree->Branch("ggzz_p0plus_VAMCFM_test",&ggzz_p0plus_VAMCFM_test,"ggzz_p0plus_VAMCFM_test/D");
  tree->Branch("qqZZ_VAMCFM_test",&qqZZ_VAMCFM_test,"qqZZ_VAMCFM_test/D");
  tree->Branch("qqZZ_VAMCFMNorm_test",&qqZZ_VAMCFMNorm_test,"qqZZ_VAMCFMNorm_test/D");
  tree->Branch("ggZZ_VAMCFM_test",&ggZZ_VAMCFM_test,"ggZZ_VAMCFM_test/D");

  tree->Branch("ggzz_VAMCFM_test",&ggzz_VAMCFM,"ggzz_VAMCFM_test/D");
  tree->Branch("ggzz_c5_VAMCFM_test",&ggzz_c5_VAMCFM,"ggzz_c5_VAMCFM_test/D");
  tree->Branch("ggzz_c1_VAMCFM_test",&ggzz_c1_VAMCFM,"ggzz_c1_VAMCFM_test/D");
  tree->Branch("ggzz_ci_VAMCFM_test",&ggzz_ci_VAMCFM,"ggzz_ci_VAMCFM_test/D");

  tree->Branch("p0plus_VAJHU_test",&p0plus_VAJHU_test,"p0plus_VAJHU_test/D");
  tree->Branch("p0minus_VAJHU_test",&p0minus_VAJHU_test,"p0minus_VAJHU_test/D");
  tree->Branch("p0hplus_VAJHU_test",&p0hplus_VAJHU_test,"p0hplus_VAJHU_test/D");
  tree->Branch("p1_VAJHU_test",&p1_VAJHU_test,"p1_VAJHU_test/D");
  tree->Branch("p1plus_VAJHU_test",&p1plus_VAJHU_test,"p1plus_VAJHU_test/D");
  tree->Branch("p2_VAJHU_test",&p2_VAJHU_test,"p2_VAJHU_test/D");
  tree->Branch("p2qqb_VAJHU_test",&p2qqb_VAJHU_test,"p2qqb_VAJHU_test/D");

  tree->Branch("p2_mad",&p2_mad,"p2_mad/D");
  tree->Branch("p0plus_mad",&p0plus_mad,"p0plus_mad/D");
  tree->Branch("p0minus_mad",&p0minus_mad,"p0minus_mad/D");
  tree->Branch("qqZZ_mad",&qqZZ_mad,"qqZZ_mad/D");

  tree->Branch("KDsigbkg",&KDsigbkg,"KDsigbkg/D");
  tree->Branch("KD0minus",&KD0minus,"KD0minus/D");
  tree->Branch("KD0hplus",&KD0hplus,"KD0hplus/D");
  tree->Branch("KD2mplus",&KD2mplus,"KD2mplus/D");
	double pg1g1prime2_VAJHU, p0_g1prime2_VAJHU, pzzzg_VAJHU, pzzgg_VAJHU, p0_Zgs_VAJHU, p0_gsgs_VAJHU;
  tree->Branch("pg1g4_mela_test",&pg1g4_mela,"pg1g4_mela_test/D");
  tree->Branch("pg1g4_VAJHU_test",&pg1g4_VAJHU,"pg1g4_VAJHU_test/D");
  tree->Branch("pg1g1prime2_VAJHU_test",&pg1g1prime2_VAJHU,"pg1g1prime2_VAJHU_test/D");
  tree->Branch("p0_g1prime2_VAJHU_test",&p0_g1prime2_VAJHU,"p0_g1prime2_VAJHU_test/D");
  tree->Branch("pg1g2_mela_test",&pg1g2_mela,"pg1g2_mela_test/D");
  tree->Branch("pg1g2_VAJHU_test",&pg1g2_VAJHU,"pg1g2_VAJHU_test/D");
  tree->Branch("pg1g2_pi2_VAJHU_test",&pg1g2_pi2_VAJHU,"pg1g2_pi2_VAJHU_test/D");
  tree->Branch("pg1g4_pi2_VAJHU_test",&pg1g4_pi2_VAJHU,"pg1g4_pi2_VAJHU_test/D");
	
  tree->Branch("pzzzg_VAJHU_test",&pzzzg_VAJHU,"pzzzg_VAJHU_test/D");
  tree->Branch("pzzgg_VAJHU_test",&pzzgg_VAJHU,"pzzgg_VAJHU_test/D");
  tree->Branch("p0_Zgs_VAJHU_test",&p0_Zgs_VAJHU,"p0_Zgs_VAJHU_test/D");
  tree->Branch("p0_gsgs_VAJHU_test",&p0_gsgs_VAJHU,"p0_gsgs_VAJHU_test/D");
	double p2h3_VAJHU_test, p2h3_qq_VAJHU_test,p2h3_decay_VAJHU_test;
	double p2h2_VAJHU_test, p2h2_qq_VAJHU_test,p2h2_decay_VAJHU_test;
	double p2h6_VAJHU_test, p2h6_qq_VAJHU_test,p2h6_decay_VAJHU_test;
	double p2h7_VAJHU_test, p2h7_qq_VAJHU_test,p2h7_decay_VAJHU_test;
	double p2h9_VAJHU_test, p2h9_qq_VAJHU_test,p2h9_decay_VAJHU_test;
	double p2h10_VAJHU_test, p2h10_qq_VAJHU_test,p2h10_decay_VAJHU_test;
  tree->Branch("p2h3_VAJHU_test",&p2h3_VAJHU_test,"p2h3_VAJHU_test/D");
  tree->Branch("p2h2_VAJHU_test",&p2h2_VAJHU_test,"p2h2_VAJHU_test/D");
  tree->Branch("p2h6_VAJHU_test",&p2h6_VAJHU_test,"p2h6_VAJHU_test/D");
  tree->Branch("p2h7_VAJHU_test",&p2h7_VAJHU_test,"p2h7_VAJHU_test/D");
  tree->Branch("p2h9_VAJHU_test",&p2h9_VAJHU_test,"p2h9_VAJHU_test/D");
  tree->Branch("p2h10_VAJHU_test",&p2h10_VAJHU_test,"p2h10_VAJHU_test/D");

  tree->Branch("p2h3_qq_VAJHU_test",&p2h3_qq_VAJHU_test,"p2h3_qq_VAJHU_test/D");
  tree->Branch("p2h2_qq_VAJHU_test",&p2h2_qq_VAJHU_test,"p2h2_qq_VAJHU_test/D");
  tree->Branch("p2h6_qq_VAJHU_test",&p2h6_qq_VAJHU_test,"p2h6_qq_VAJHU_test/D");
  tree->Branch("p2h7_qq_VAJHU_test",&p2h7_qq_VAJHU_test,"p2h7_qq_VAJHU_test/D");
  tree->Branch("p2h9_qq_VAJHU_test",&p2h9_qq_VAJHU_test,"p2h9_qq_VAJHU_test/D");
  tree->Branch("p2h10_qq_VAJHU_test",&p2h10_qq_VAJHU_test,"p2h10_qq_VAJHU_test/D");

  tree->Branch("p2h3_decay_VAJHU_test",&p2h3_decay_VAJHU_test,"p2h3_decay_VAJHU_test/D");
  tree->Branch("p2h2_decay_VAJHU_test",&p2h2_decay_VAJHU_test,"p2h2_decay_VAJHU_test/D");
  tree->Branch("p2h6_decay_VAJHU_test",&p2h6_decay_VAJHU_test,"p2h6_decay_VAJHU_test/D");
  tree->Branch("p2h7_decay_VAJHU_test",&p2h7_decay_VAJHU_test,"p2h7_decay_VAJHU_test/D");
  tree->Branch("p2h9_decay_VAJHU_test",&p2h9_decay_VAJHU_test,"p2h9_decay_VAJHU_test/D");
  tree->Branch("p2h10_decay_VAJHU_test",&p2h10_decay_VAJHU_test,"p2h10_decay_VAJHU_test/D");
	double p0_self_VAJHU;
  tree->Branch("p0_self_VAJHU",&p0_self_VAJHU,"p0_self_VAJHU/D");
	double pggHZZ_10;
  tree->Branch("pggHZZ_10",&pggHZZ_10,"pggHZZ_10/D");
  // old branches

  float mzz,m1,m2,h1,h2,hs,phi,phi1,pt4l,Y4l;

  t->SetBranchAddress("ZZMass",&mzz);
  t->SetBranchAddress("Z1Mass",&m1);
  t->SetBranchAddress("Z2Mass",&m2);
  t->SetBranchAddress("helcosthetaZ1",&h1);
  t->SetBranchAddress("helcosthetaZ2",&h2);
  t->SetBranchAddress("costhetastar",&hs);
  t->SetBranchAddress("helphi",&phi);
  t->SetBranchAddress("phistarZ1",&phi1);
  t->SetBranchAddress("ZZPt",&pt4l);
  t->SetBranchAddress("ZZRapidity",&Y4l);

  vector<TLorentzVector> p4;
  vector<int> id;

  // channel 1: 4e 2: 4mu 2: 2e2mu
  if(channel==1){
    id.push_back(11);
    id.push_back(-11);
    id.push_back(11);
    id.push_back(-11);  
  }
  if(channel==2){
    id.push_back(13);
    id.push_back(-13);
    id.push_back(13);
    id.push_back(-13);  
  }
  if(channel==3){
    id.push_back(11);
    id.push_back(-11);
    id.push_back(13);
    id.push_back(-13);  
  }

  //for(int i=0 ; i<t->GetEntries(); i++){
  for(int i=0 ; i<500; i++){

    t->GetEntry(i);
    //t->Show(i);
    if ( i % 100 == 0 ) cout << i << "/" << t->GetEntries() << endl;

    if(mzz>140) continue;

    p4 = Calculate4Momentum(mzz,m1,m2,acos(hs),acos(h1),acos(h2),phi1,phi,pt4l,Y4l);
    
    test.computeME(kSMHiggs,kAnalytical,p4,id,p0plus_mela_test);
    test.computeME(kSMHiggs,kMELA_HCP,p4,id,p0plus_melaNorm_test);
    test.computeME(kqqZZ,kMELA_HCP,p4,id,bkg_mela_test);

    test.computeME(k0hplus,kAnalytical,p4,id,p0hplus_mela_test);
    test.computeME(k0minus,kAnalytical,p4,id,p0minus_mela_test);
    test.computeME(k1minus,kAnalytical,p4,id,p1_mela_test);
    test.computeME(k1plus,kAnalytical,p4,id,p1plus_mela_test);
    test.computeME(k2mplus_gg,kAnalytical,p4,id,p2_mela_test);
    test.computeME(k2mplus_qqbar,kAnalytical,p4,id,p2qqb_mela_test);
   
    test.computeME(kggZZ,kMCFM,p4,id,ggZZ_VAMCFM_test);
    test.computeME(kSMHiggs,kMCFM,p4,id,p0plus_VAMCFM_test);
    test.computeME(kggZZ_SMHiggs,kMCFM,p4,id,ggzz_p0plus_VAMCFM_test);
    test.computeME(kqqZZ,kMCFM,p4,id,qqZZ_VAMCFM_test);

    test.computeME(kggHZZ_10,kMCFM,p4,id,pggHZZ_10);

    test.computeME(kSMHiggs,kJHUGen,p4,id,p0plus_VAJHU_test);
    test.computeME(k0hplus,kJHUGen,p4,id,p0hplus_VAJHU_test);
    test.computeME(k0minus,kJHUGen,p4,id,p0minus_VAJHU_test);
    test.computeME(k0_g1prime2,kJHUGen,p4,id,p0_g1prime2_VAJHU);
    test.computeME(k0_Zgs,kJHUGen,p4,id,p0_Zgs_VAJHU);
    test.computeME(k0_gsgs,kJHUGen,p4,id,p0_gsgs_VAJHU);
    test.computeME(k1minus,kJHUGen,p4,id,p1_VAJHU_test);
    test.computeME(k1plus,kJHUGen,p4,id,p1plus_VAJHU_test);
    test.computeME(k2mplus_gg,kJHUGen,p4,id,p2_VAJHU_test);
    test.computeME(k2mplus_qqbar,kJHUGen,p4,id,p2qqb_VAJHU_test);

    test.computeME(k2h2plus_gg,kJHUGen,p4,id,p2h2_VAJHU_test);
    test.computeME(k2h2plus_qqbar,kJHUGen,p4,id,p2h2_qq_VAJHU_test);
    test.computeME(k2h2plus_prodIndep,kJHUGen,p4,id,p2h2_decay_VAJHU_test);

    test.computeME(k2h2plus_gg,kJHUGen,p4,id,p2h2_VAJHU_test);
    test.computeME(k2h2plus_qqbar,kJHUGen,p4,id,p2h2_qq_VAJHU_test);
    test.computeME(k2h2plus_prodIndep,kJHUGen,p4,id,p2h2_decay_VAJHU_test);

    test.computeME(k2h3plus_gg,kJHUGen,p4,id,p2h3_VAJHU_test);
    test.computeME(k2h3plus_qqbar,kJHUGen,p4,id,p2h3_qq_VAJHU_test);
    test.computeME(k2h3plus_prodIndep,kJHUGen,p4,id,p2h3_decay_VAJHU_test);

    test.computeME(k2h6plus_gg,kJHUGen,p4,id,p2h6_VAJHU_test);
    test.computeME(k2h6plus_qqbar,kJHUGen,p4,id,p2h6_qq_VAJHU_test);
    test.computeME(k2h6plus_prodIndep,kJHUGen,p4,id,p2h6_decay_VAJHU_test);

    test.computeME(k2h7plus_gg,kJHUGen,p4,id,p2h7_VAJHU_test);
    test.computeME(k2h7plus_qqbar,kJHUGen,p4,id,p2h7_qq_VAJHU_test);
    test.computeME(k2h7plus_prodIndep,kJHUGen,p4,id,p2h7_decay_VAJHU_test);

    test.computeME(k2h9minus_gg,kJHUGen,p4,id,p2h9_VAJHU_test);
    test.computeME(k2h9minus_qqbar,kJHUGen,p4,id,p2h9_qq_VAJHU_test);
    test.computeME(k2h9minus_prodIndep,kJHUGen,p4,id,p2h9_decay_VAJHU_test);

    test.computeME(k2h10minus_gg,kJHUGen,p4,id,p2h10_VAJHU_test);
    test.computeME(k2h10minus_qqbar,kJHUGen,p4,id,p2h10_qq_VAJHU_test);
    test.computeME(k2h10minus_prodIndep,kJHUGen,p4,id,p2h10_decay_VAJHU_test);

    test.computeME(kSMHiggs,kMEKD,p4,id,p0plus_mad);
    test.computeME(k2mplus_gg,kMEKD,p4,id,p2_mad);
    test.computeME(k0minus,kMEKD,p4,id,p0minus_mad);
    test.computeME(kqqZZ,kMEKD,p4,id,qqZZ_mad);

    test.computeME_Interference(kg1g4,kAnalytical,p4,id,pg1g4_mela);
    test.computeME_Interference(kg1g4,kJHUGen,p4,id,pg1g4_VAJHU);
    test.computeME_Interference(kg1g4_pi_2,kJHUGen,p4,id,pg1g4_pi2_VAJHU);
    test.computeME_Interference(kg1g2_pi_2,kJHUGen,p4,id,pg1g2_pi2_VAJHU);
    test.computeME_Interference(kg1g2,kAnalytical,p4,id,pg1g2_mela);
    test.computeME_Interference(kg1g2,kJHUGen,p4,id,pg1g2_VAJHU);
    test.computeME_Interference(k_g1g1prime2,kJHUGen,p4,id,pg1g1prime2_VAJHU);
    test.computeME_Interference(kzzzg,kJHUGen,p4,id,pzzzg_VAJHU);
    test.computeME_Interference(kzzgg,kJHUGen,p4,id,pzzgg_VAJHU);
   	vector<complex<double> > *prod = new vector<complex<double> >;
   	vector<complex<double> >* decay = new vector<complex<double> >;
		decay->push_back(1.0);
		decay->push_back(0.0);
		decay->push_back(0.0);
		decay->push_back(0.0);
    test.computeME(kSpin0_gg,kJHUGen,p4,id,prod,decay, p0_self_VAJHU);
    qqZZ_VAMCFMNorm_test=test.qqZZ_MCFMNorm;

    //test.computeMEs(p4,id);

    vector<complex<double> > *coupling = new vector<complex<double> >;
    vector<complex<double> > *couplingprod= new vector<complex<double> >;
    test.computeME(MEMNames::kggZZ, MEMNames::kMCFM, p4, id, ggzz_VAMCFM);
    complex<double> coup(1.,0.);
    coupling->push_back(coup);
    //cout<<(*coupling).size()<<endl;

    test.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kJHUGen, p4, id, couplingprod,coupling, ggzz_c1_VAMCFM);
    coupling->clear();
    coup=5.;
    coupling->push_back(coup);
    test.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kJHUGen, p4, id, couplingprod,coupling, ggzz_c5_VAMCFM);
    coupling->clear();
    coup.real(0.);
    coup.imag(1.);
    coupling->push_back(coup);
    test.computeME(MEMNames::kggZZ_SMHiggs, MEMNames::kJHUGen, p4, id, couplingprod,coupling, ggzz_ci_VAMCFM);

    test.computeKD(kSMHiggs,kJHUGen,kqqZZ,kMCFM,&MEMs::probRatio,KDsigbkg,dummy1,dummy2);
    test.computeKD(kSMHiggs,kJHUGen,k0minus,kJHUGen,&MEMs::probRatio,KD0minus,dummy1,dummy2);
    test.computeKD(kSMHiggs,kJHUGen,k0hplus,kJHUGen,&MEMs::probRatio,KD0hplus,dummy1,dummy2);
    test.computeKD(kSMHiggs,kJHUGen,k2mplus_gg,kJHUGen,&MEMs::probRatio,KD2mplus,dummy1,dummy2);
    double KD;
    test.computeKD(kSMHiggs,kMEKD,kqqZZ,kMCFM,&MEMs::probRatio,KD,dummy1,dummy2);

    if(KD>=1.0 || KD<=0.0)
      cout << "KD: " << KD << endl;
    /*
    cout << "p0plusJHU: " << p0plus_VAJHU_test << endl;
    cout << "qqZZMCFM: " << qqZZ_VAMCFMNorm_test << endl;
    cout << "sigbkg: " << KDsigbkg << endl;
    */

    /*
    cout << "p0plus: " << p0plus_mad << endl;
    cout << "p0minus: " << p0minus_mad << endl;
    cout << "p2: " << p2_mad << endl;
    cout << "qqZZ: " << qqZZ_mad << endl;

    double KD,pA,pB;
    test.computeKD(kSMHiggs,k0minus,kAnalytical,p4,id,KD,pA,pB);
    cout << "pseudoMELA: " << KD << endl;
    test.computeKD(kSMHiggs,k2mplus_gg,kAnalytical,p4,id,KD,pA,pB);
    cout << "graviMELA: " << KD << endl;
    test.computeKD(kSMHiggs,kqqZZ,kMELA_HCP,p4,id,KD,pA,pB);
    cout << "MELA: " << KD << endl;
    */

    tree->Fill();

  } // end loop over events

  file->cd();
  tree->Write();
  file->Close();

}


vector<TLorentzVector> Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi,double pt4l, double Y4l)
{
  double phi1,phi2;
  phi1=TMath::Pi()-Phi1;
  phi2=Phi1+Phi;
    
    
  double gamma1,gamma2,beta1,beta2;
    
  gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
  gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
  beta1=sqrt(1-1/(gamma1*gamma1));
  beta2=sqrt(1-1/(gamma2*gamma2));
    
    
  //gluon 4 vectors
  TLorentzVector p1CM(0,0,Mx/2,Mx/2);
  TLorentzVector p2CM(0,0,-Mx/2,Mx/2);
    
  //vector boson 4 vectors
  TLorentzVector kZ1(gamma1*M1*sin(theta)*beta1,0, gamma1*M1*cos(theta)*beta1,gamma1*M1*1);   
  TLorentzVector kZ2(-gamma2*M2*sin(theta)*beta2,0, -gamma2*M2*cos(theta)*beta2,gamma2*M2*1);
    
  //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.
  
  TLorentzRotation Z1ToZ,Z2ToZ;
    
  Z1ToZ.Boost(0,0,beta1);
  Z2ToZ.Boost(0,0,beta2);
  Z1ToZ.RotateY(theta);
  Z2ToZ.RotateY(TMath::Pi()+theta);
    
    
  //fermions 4 vectors in vector boson rest frame
    
  TLorentzVector p3Z1((M1/2)*sin(theta1)*cos(phi1),(M1/2)*sin(theta1)*sin(phi1),(M1/2)*cos(theta1),(M1/2)*1);
       
  TLorentzVector p4Z1(-(M1/2)*sin(theta1)*cos(phi1),-(M1/2)*sin(theta1)*sin(phi1),-(M1/2)*cos(theta1),(M1/2)*1);
      
  TLorentzVector p5Z2((M2/2)*sin(theta2)*cos(phi2),(M2/2)*sin(theta2)*sin(phi2),(M2/2)*cos(theta2),(M2/2)*1);
    
  TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2),-(M2/2)*sin(theta2)*sin(phi2),-(M2/2)*cos(theta2),(M2/2)*1);
      
    
  // fermions 4 vectors in CM frame
    
  TLorentzVector p3CM,p4CM,p5CM,p6CM;
    
  p3CM=Z1ToZ*p3Z1;
  p4CM=Z1ToZ*p4Z1;
  p5CM=Z2ToZ*p5Z2;
  p6CM=Z2ToZ*p6Z2;

  // boost to lab frame:
  //TLorentzRotation labFrame;
  //labFrame.Boost(pt4l/Mx,0,sqrt(1.-1./(cosh(Y4l)*cosh(Y4l)));
		 
  vector<TLorentzVector> p;
  /*
*labFrame
*labFrame
*labFrame
*labFrame
  */
  p.push_back(p3CM); 
  p.push_back(p4CM);
  p.push_back(p5CM);
  p.push_back(p6CM);

  return p;
}

