#include <MuonMVAReader/Reader/interface/MuonTMVAReader.hpp>


//================================
MuonTMVAReader::MuonTMVAReader( int year )
{
   TMVAReader_pT_5  = new TMVA::Reader("!Color:!Silent:Error" );
   TMVAReader_pT_10 = new TMVA::Reader("!Color:!Silent:Error" );

   if ( year == 2016 )
   {
      fip_pT_5_  = edm::FileInPath("MuonMVAReader/Reader/data/MVAWeightFiles/Summer_16_ID_ISO_SIP/pT_5/weights.xml");
      fip_pT_10_ = edm::FileInPath("MuonMVAReader/Reader/data/MVAWeightFiles/Summer_16_ID_ISO_SIP/pT_10/weights.xml");
   }
   
//   else if ( year == 2017 )
//   {
//      fip_pT_5_  = edm::FileInPath("MuonMVAReader/Reader/data/MVAWeightFiles/...");
//      fip_pT_10_ = edm::FileInPath("MuonMVAReader/Reader/data/MVAWeightFiles/...");
//   }
//
//   else if ( year == 2018 )
//   {
//      fip_pT_5_  = edm::FileInPath("MuonMVAReader/Reader/data/MVAWeightFiles/...");
//      fip_pT_10_ = edm::FileInPath("MuonMVAReader/Reader/data/MVAWeightFiles/...");
//   }
   
   else
   {
      std::cerr << "[ERROR] MuonMVAReader: Year " << year << " is currently not supported!" << std::endl;
   }
   
   weights_pT_5_  = fip_pT_5_.fullPath();
   weights_pT_10_ = fip_pT_10_.fullPath();
}
//================================


//=================================
MuonTMVAReader::~MuonTMVAReader()
{
}
//=================================


//================================
void MuonTMVAReader::Initialize()
{
   TMVAReader_pT_5->AddVariable("mu_N_hits",             &mu_N_hits_);
   TMVAReader_pT_5->AddVariable("mu_N_pixel_hits",       &mu_N_pixel_hits_);
   TMVAReader_pT_5->AddVariable("mu_N_tracker_hits",     &mu_N_tracker_hits_);
   TMVAReader_pT_5->AddVariable("mu_chi_square",         &mu_chi_square_);
   TMVAReader_pT_5->AddVariable("mu_pf_photon_iso",      &mu_pf_photon_iso_);
   TMVAReader_pT_5->AddVariable("mu_pf_charged_had_iso", &mu_pf_charged_had_iso_);
   TMVAReader_pT_5->AddVariable("mu_pf_neutral_had_iso", &mu_pf_neutral_had_iso_); 
   TMVAReader_pT_5->AddVariable("mu_rho",                &mu_rho_);
   TMVAReader_pT_5->AddVariable("mu_sip",                &mu_sip_);
   TMVAReader_pT_5->AddVariable("mu_dxy",                &mu_dxy_);
   TMVAReader_pT_5->AddVariable("mu_dz",                 &mu_dz_);
   
   TMVAReader_pT_10->AddVariable("mu_N_hits",             &mu_N_hits_);
   TMVAReader_pT_10->AddVariable("mu_N_pixel_hits",       &mu_N_pixel_hits_);
   TMVAReader_pT_10->AddVariable("mu_N_tracker_hits",     &mu_N_tracker_hits_);
   TMVAReader_pT_10->AddVariable("mu_chi_square",         &mu_chi_square_);
   TMVAReader_pT_10->AddVariable("mu_pf_photon_iso",      &mu_pf_photon_iso_);
   TMVAReader_pT_10->AddVariable("mu_pf_charged_had_iso", &mu_pf_charged_had_iso_);
   TMVAReader_pT_10->AddVariable("mu_pf_neutral_had_iso", &mu_pf_neutral_had_iso_); 
   TMVAReader_pT_10->AddVariable("mu_rho",                &mu_rho_);
   TMVAReader_pT_10->AddVariable("mu_sip",                &mu_sip_);
   TMVAReader_pT_10->AddVariable("mu_dxy",                &mu_dxy_);
   TMVAReader_pT_10->AddVariable("mu_dz",                 &mu_dz_);

//   TMVAReader->AddVariable("is_tracker_mu",         &is_tracker_mu_);
//   TMVAReader->AddVariable("is_global_mu",          &is_global_mu_);
//   TMVAReader->AddVariable("mu_N_matched_stations", &mu_N_matched_stations_);

   
   TMVAReader_pT_5->BookMVA("BDTG", weights_pT_5_);
   TMVAReader_pT_10->BookMVA("BDTG", weights_pT_10_);
   cout << "[INFO] MVA loaded with weigths: " << weights_pT_5_ << " and " << weights_pT_10_ << endl;
   
   isInitialized_ = true;
}
//================================


//====================================================================
int MuonTMVAReader::isInitialized() const { return isInitialized_; };
//====================================================================


//============================================================================================================================================================================
float MuonTMVAReader::Get_MVA_value( float mu_pT, float mu_eta, float mu_N_hits, float mu_N_pixel_hits, float mu_N_tracker_hits, float mu_chi_square, float mu_pf_photon_iso,
                                      float mu_pf_charged_had_iso, float mu_pf_neutral_had_iso, float mu_rho, float mu_sip, float mu_dxy, float mu_dz )
{
   if ( !isInitialized() )
   { 
      cout << "[Error] MuonMVAEstimator not properly initialized!" << endl; 
   }

   mu_N_hits_             = mu_N_hits;
   mu_N_pixel_hits_       = mu_N_pixel_hits;
   mu_N_tracker_hits_     = mu_N_tracker_hits;
   mu_chi_square_         = mu_chi_square;
   mu_pf_photon_iso_      = mu_pf_photon_iso;
   mu_pf_charged_had_iso_ = mu_pf_charged_had_iso;
   mu_pf_neutral_had_iso_ = mu_pf_neutral_had_iso;
   mu_rho_                = mu_rho;
   mu_sip_                = mu_sip;
   mu_dxy_                = mu_dxy;
   mu_dz_                 = mu_dz;
   
//   is_tracker_mu_         = is_tracker_mu;
//   is_global_mu_          = is_global_mu;
//   mu_N_matched_stations_ = mu_N_matched_stations;
   
   if ( mu_pT < 10)
   {
      return TMVAReader_pT_5->EvaluateMVA("BDTG");
   }
   else
   {
      return TMVAReader_pT_10->EvaluateMVA("BDTG");   
   }
}
//============================================================================================================================================================================