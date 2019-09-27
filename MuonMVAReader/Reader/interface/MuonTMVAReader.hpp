#ifndef MuonMVAReader_hpp
#define MuonMVAReader_hpp

#include <vector>
#include <TROOT.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

using namespace std;

class MuonTMVAReader
{
   public:
      MuonTMVAReader( int );
      ~MuonTMVAReader(); 
  
      void Initialize();
      int isInitialized() const;
      float Get_MVA_value( float, float, float, float, float, float, float, float, float, float, float, float, float );
    
   private:
      TMVA::Reader *TMVAReader_pT_5;
      TMVA::Reader *TMVAReader_pT_10;
      
      edm::FileInPath fip_pT_5_;
      edm::FileInPath fip_pT_10_;
      
      int isInitialized_;
      TString weights_location_, weights_pT_5_, weights_pT_10_;
      
      float mu_pT_;
      float mu_eta_;
      float mu_N_hits_;
      float mu_N_pixel_hits_;
      float mu_N_tracker_hits_;
      float mu_chi_square_;
      float mu_pf_photon_iso_;
      float mu_pf_charged_had_iso_;
      float mu_pf_neutral_had_iso_;
      float mu_rho_;
      float mu_sip_;
      float mu_dxy_;
      float mu_dz_;
      
//      float is_tracker_mu_;
//      float is_global_mu_;
//      float mu_N_matched_stations_;
};

#endif