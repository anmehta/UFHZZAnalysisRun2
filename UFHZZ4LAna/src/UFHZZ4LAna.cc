// -*- C++ -*-
//
// Package:    UFHZZ4LAna
// Class:      UFHZZ4LAna
//
//

// system include files
#include <memory>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iomanip>
#include <set>
#include<TApplication.h>

#define PI 3.14159

// user include files
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TSpline.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
#include "TClonesArray.h"
#include "TCanvas.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

//HTXS
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"
//#include "SimDataFormats/HZZFiducial/interface/HZZFiducialVolume.h"

// PAT
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Reco
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

// KD's
#include "JHUGenMELA/MELA/interface/Mela.h"

//Helper
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LHelper.h"
//Muons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMuonAna.h"
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LMuonTree.h"
//Electrons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LElectronTree.h"
//Photons
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPhotonTree.h"
//Jets
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LJetTree.h"
//Final Leps
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LFinalLepTree.h"
//Sip
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LSipAna.h"
//PU
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LPileUp.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//GEN
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LGENAna.h"
//VBF Jets
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/HZZ4LJets.h"

//nJettiness
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/NJettiness.h"

// Jet energy correction
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include <vector>

// Kinematic Fit
#include "KinZfitter/KinZfitter/interface/KinZfitter.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyResolution.h"

// EWK corrections
#include "UFHZZAnalysisRun2/UFHZZ4LAna/interface/EwkCorrections.h"

// JEC related
#include "PhysicsTools/PatAlgos/plugins/PATJetUpdater.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

//JER related
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//BTag Calibration

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

//Muon MVA
//#include "MuonMVAReader/Reader/interface/MuonGBRForestReader.hpp"

// KalmanVertexFitter
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
// Rochester Corrections
#include "UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/src/RoccoR.cc"

#include "RecoVertex/KalmanVertexFit/interface/SingleTrackVertexConstraint.h"

//Auto MELA
#include "IvyBase.h"
#include <IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h>
#include <IvyFramework/IvyAutoMELA/interface/IvyMELAOutputStreamerExt.h>
#include <IvyFramework/IvyDataTools/interface/HostHelpersCore.h>

//
// class declaration
//
using namespace IvyStreamHelpers;
using namespace EwkCorrections;

class UFHZZ4LAna : public edm::EDAnalyzer{
public:
  explicit UFHZZ4LAna(const edm::ParameterSet&);
  ~UFHZZ4LAna();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static bool sortByPt( const reco::GenParticle &p1, const reco::GenParticle &p2 ){ return (p1.pt() > p2.pt()); };

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, const edm::EventSetup& iSetup);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup);

  void setMatrixElementList(std::vector<std::string> const& MElist, bool const& isGen);
  void setMatrixElementListFromFile(std::string fname, std::string const& MElistTypes, bool const& isGen); // MElistTypes is comman-separated


  RoccoR  *calibrator;
  float ApplyRoccoR(int Y, bool isMC, int charge, float pt, float eta, float phi, float genPt, float nLayers);

  //Helper Class
  HZZ4LHelper helper;
  //GEN
  HZZ4LGENAna genAna;
  //VBF
  HZZ4LJets jetHelper;
  //PU Reweighting
  edm::LumiReWeighting *lumiWeight;
  HZZ4LPileUp pileUp;
  //JES Uncertainties
  std::unique_ptr<JetCorrectionUncertainty> jecunc;
  std::unique_ptr<JetCorrectionUncertainty> jecmergedunc;
  // kfactors
  TSpline3 *kFactor_ggzz;
  std::vector<std::vector<float> > tableEwk;
  // data/MC scale factors
  TH2F *hElecScaleFac;
  TH2F *hElecScaleFac_Cracks;
  TH2F *hElecScaleFacGsf;
  TH2F *hElecScaleFacGsfLowET;
  TH2F *hMuScaleFac;
  TH2F *hMuScaleFacUnc;
  TH1D *h_pileup;
  TH1D *h_pileupUp;
  TH1D *h_pileupDn;
  std::vector<TH1F*> h_medians;
  TH2F *hbTagEffi;
  TH2F *hcTagEffi;
  TH2F *hudsgTagEffi;

  BTagCalibrationReader* reader;

  //Saved Events Trees
  TTree *passedEventsTree_All;
  void bookPassedEventTree(TString treeName, TTree *tree);
  void setTreeVariables( const edm::Event&, const edm::EventSetup&,
                        std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons,
                        std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons,
                        std::vector<pat::Jet> goodJets, std::vector<float> goodJetQGTagger,
                        std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult,
                        std::vector<pat::Jet> selectedMergedJets,
                        std::map<unsigned int, TLorentzVector> selectedFsrMap,
                        edm::Handle<pat::PackedCandidateCollection> &pfcands,
                        edm::Handle<edm::View<pat::Jet> > &jets);
  void setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                       edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                       edm::Handle<edm::View<reco::GenJet> > genJets);
  // -------------------------
  // RECO level information
  // -------------------------

  // Event Variables
  ULong64_t Run, Event, LumiSect;
  int nVtx, nInt;
  int finalState;
  std::string triggersPassed;
  bool passedTrig;

  float PV_x, PV_y, PV_z;
  float BS_x, BS_y, BS_z;
  float BS_xErr, BS_yErr, BS_zErr;
  float BeamWidth_x, BeamWidth_y;
  float BeamWidth_xErr, BeamWidth_yErr;

  // Event Weights
  float genWeight, pileupWeight, pileupWeightUp, pileupWeightDn, dataMCWeight, eventWeight, prefiringWeight, prefiringWeightECAL, prefiringWeightMuon;
  // pdf weights
  vector<float> qcdWeights;
  vector<float> nnloWeights;
  vector<float> pdfWeights;
  int posNNPDF;
  float pdfRMSup, pdfRMSdown, pdfENVup, pdfENVdown;
  // lepton variables
  vector<double> lep_pt_FromMuonBestTrack, lep_eta_FromMuonBestTrack, lep_phi_FromMuonBestTrack;
  vector<double> lep_position_x, lep_position_y, lep_position_z;
  vector<double> lep_pt_genFromReco;
  vector<double> lep_pt; vector<double> lep_pterr; vector<double> lep_pterrold;
  vector<double> lep_p; vector<double> lep_ecalEnergy; vector<int> lep_isEB; vector<int> lep_isEE;
  vector<double> lep_eta; vector<double> lep_phi; vector<double> lep_mass;
  vector<double> lepFSR_pt; vector<double> lepFSR_eta; vector<double> lepFSR_phi; vector<double> lepFSR_mass; vector<int> lepFSR_ID;

  vector<double> lep_errPre_Scale, lep_errPost_Scale, lep_errPre_noScale, lep_errPost_noScale;
  vector<double> lep_pt_UnS, lep_pterrold_UnS;

  vector<float> lep_dataMC; vector<float> lep_dataMCErr;

  vector<float> lep_d0BS;
  vector<float> lep_d0PV;
  vector<float> lep_numberOfValidPixelHits;
  vector<float> lep_trackerLayersWithMeasurement;
  vector<float> dataMC_VxBS; vector<float> dataMCErr_VxBS;
  vector<int> lep_genindex; //position of lepton in GENlep_p4 (if gen matched, -1 if not gen matched)
  vector<int> lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId; // gen matching even if not in GENlep_p4
  vector<int> lep_id;
  vector<float> lep_mva; vector<int> lep_ecalDriven;
  vector<int> lep_tightId; vector<int> lep_tightIdSUS; vector<int> lep_tightIdHiPt; //vector<int> lep_tightId_old;
  vector<float> lep_Sip; vector<float> lep_IP; vector<float> lep_isoNH; vector<float> lep_isoCH; vector<float> lep_isoPhot;
  vector<float> lep_isoPU; vector<float> lep_isoPUcorr;
  vector<float> lep_RelIso; vector<float> lep_RelIsoNoFSR; vector<float> lep_MiniIso;
  vector<float> lep_ptRatio; vector<float> lep_ptRel;
  vector<int> lep_missingHits;
  vector<string> lep_filtersMatched; // for each lepton, all filters it is matched to
  int nisoleptons;
  double muRho, elRho, rhoSUS;

  // Z candidate variables
  vector<float> singleBS_RecoLep_pt; vector<float> singleBS_RecoLep_ptError;  vector<float> singleBS_RecoLep_eta;  vector<float> singleBS_RecoLep_phi; vector<float> singleBS_RecoLep_mass; vector<float> singleBS_RecoLep_d0;

  // photon variables
  vector<double> pho_pt, pho_eta, pho_phi, photonCutBasedIDLoose;

  // MET
  float met; float met_phi;
  float met_jesup, met_phi_jesup, met_jesdn, met_phi_jesdn;
  float met_uncenup, met_phi_uncenup, met_uncendn, met_phi_uncendn;

  // Jets
  vector<int>    jet_iscleanH4l;
  int jet1index, jet2index, jet1index_2p5, jet2index_2p5;
  vector<double> jet_pt; vector<double> jet_eta; vector<double> jet_phi; vector<double> jet_mass; vector<double> jet_pt_raw;
  vector<float>  jet_pumva, jet_csvv2,  jet_csvv2_; vector<int> jet_isbtag;
  vector<int>    jet_hadronFlavour, jet_partonFlavour;
  vector<float>  jet_QGTagger, jet_QGTagger_jesup, jet_QGTagger_jesdn;
  vector<float> jet_axis2, jet_ptD; vector<int> jet_mult;
  vector<float>  jet_relpterr; vector<float>  jet_phierr;
  vector<float>  jet_bTagEffi;
  vector<float>  jet_cTagEffi;
  vector<float>  jet_udsgTagEffi;
  vector<int>    jet_jesup_iscleanH4l;
  vector<double> jet_jesup_pt; vector<double> jet_jesup_eta;
  vector<double> jet_jesup_phi; vector<double> jet_jesup_mass;
  vector<int>    jet_jesdn_iscleanH4l;
  vector<double> jet_jesdn_pt; vector<double> jet_jesdn_eta;
  vector<double> jet_jesdn_phi; vector<double> jet_jesdn_mass;
  vector<int>    jet_jerup_iscleanH4l;
  vector<double> jet_jerup_pt; vector<double> jet_jerup_eta;
  vector<double> jet_jerup_phi; vector<double> jet_jerup_mass;
  vector<int>    jet_jerdn_iscleanH4l;
  vector<double> jet_jerdn_pt; vector<double> jet_jerdn_eta;
  vector<double> jet_jerdn_phi; vector<double> jet_jerdn_mass;
  int njets_pt30_eta4p7; int njets_pt30_eta4p7_jesup; int njets_pt30_eta4p7_jesdn;
  int njets_pt30_eta4p7_jerup; int njets_pt30_eta4p7_jerdn;
  int njets_pt30_eta2p5; int njets_pt30_eta2p5_jesup; int njets_pt30_eta2p5_jesdn;
  int njets_pt30_eta2p5_jerup; int njets_pt30_eta2p5_jerdn;
  int nbjets_pt30_eta4p7; int nvjets_pt40_eta2p4;
  float pt_leadingjet_pt30_eta4p7;
  float pt_leadingjet_pt30_eta4p7_jesup; float pt_leadingjet_pt30_eta4p7_jesdn;
  float pt_leadingjet_pt30_eta4p7_jerup; float pt_leadingjet_pt30_eta4p7_jerdn;
  float pt_leadingjet_pt30_eta2p5;
  float pt_leadingjet_pt30_eta2p5_jesup; float pt_leadingjet_pt30_eta2p5_jesdn;
  float pt_leadingjet_pt30_eta2p5_jerup; float pt_leadingjet_pt30_eta2p5_jerdn;
  float absrapidity_leadingjet_pt30_eta4p7;
  float absrapidity_leadingjet_pt30_eta4p7_jesup; float absrapidity_leadingjet_pt30_eta4p7_jesdn;
  float absrapidity_leadingjet_pt30_eta4p7_jerup; float absrapidity_leadingjet_pt30_eta4p7_jerdn;
  float absdeltarapidity_hleadingjet_pt30_eta4p7;
  float absdeltarapidity_hleadingjet_pt30_eta4p7_jesup; float absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn;
  float absdeltarapidity_hleadingjet_pt30_eta4p7_jerup; float absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn;
  float DijetMass, DijetDEta, DijetFisher;

  // merged jets
  vector<int>   mergedjet_iscleanH4l;
  vector<float> mergedjet_pt, mergedjet_ptrow;
  vector<float> mergedjet_eta; vector<float> mergedjet_phi; vector<float> mergedjet_mass;
  vector<float> mergedjet_jer;
  vector<float> mergedjet_jer_pt, mergedjet_jer_eta, mergedjet_jer_phi, mergedjet_jer_mass;
  vector<float> mergedjet_jerup_pt, mergedjet_jerup_eta, mergedjet_jerup_phi, mergedjet_jerup_mass;
  vector<float> mergedjet_jerdn_pt, mergedjet_jerdn_eta, mergedjet_jerdn_phi, mergedjet_jerdn_mass;
  vector<float> mergedjet_jer_pterr, mergedjet_jer_phierr;

  vector<float> mergedjet_tau1; vector<float> mergedjet_tau2;
  vector<float> mergedjet_btag;
  vector<float> mergedjet_ZvsQCD;
  vector<float> mergedjet_ZbbvsQCD;
  vector<float> mergedjet_WvsQCD;
  vector<float> mergedjet_ZHbbvsQCD;
  vector<float> mergedjet_HbbvsQCD;
  vector<float> mergedjet_H4qvsQCD;

  vector<float> mergedjet_ZvsQCD_de; //de with decorrelated
  vector<float> mergedjet_ZbbvsQCD_de;
  vector<float> mergedjet_WvsQCD_de;
  vector<float> mergedjet_ZHbbvsQCD_de;
  vector<float> mergedjet_HbbvsQCD_de;
  vector<float> mergedjet_H4qvsQCD_de;

  vector<float> mergedjet_Net_Xbb_de; //particle net
  vector<float> mergedjet_Net_Xcc_de;
  vector<float> mergedjet_Net_Xqq_de;
  vector<float> mergedjet_Net_QCDbb_de;
  vector<float> mergedjet_Net_QCDcc_de;
  vector<float> mergedjet_Net_QCDb_de;
  vector<float> mergedjet_Net_QCDc_de;
  vector<float> mergedjet_Net_QCDother_de;

  //vector<float> mergedjet_L1;
  //vector<float> mergedjet_prunedmass;
  vector<float> mergedjet_softdropmass;

  vector<int> mergedjet_nsubjet;
  vector<float> mergedjet_subjet_softDropMass;
  vector<float> mergedjet_subjet_softDropMassUncorr;
  vector<vector<float> > mergedjet_subjet_pt; vector<vector<float> > mergedjet_subjet_eta;
  vector<vector<float> > mergedjet_subjet_phi; vector<vector<float> > mergedjet_subjet_mass;
  vector<vector<float> > mergedjet_subjetUncorr_pt; vector<vector<float> > mergedjet_subjetUncorr_eta;
  vector<vector<float> > mergedjet_subjetUncorr_phi; vector<vector<float> > mergedjet_subjetUncorr_mass;
  vector<vector<float> > mergedjet_subjet_btag;
  vector<vector<int> > mergedjet_subjet_partonFlavour, mergedjet_subjet_hadronFlavour;
  vector<int> mergedjet_nbHadrons, mergedjet_ncHadrons, mergedjet_hadronFlavour, mergedjet_partonFlavour;

  // FSR Photons
  int nFSRPhotons;
  vector<int> fsrPhotons_lepindex;
  vector<double> fsrPhotons_pt; vector<double> fsrPhotons_pterr;
  vector<double> fsrPhotons_eta; vector<double> fsrPhotons_phi;
  vector<double> fsrPhotons_mass;
  vector<float> fsrPhotons_dR; vector<float> fsrPhotons_iso;
  vector<float> allfsrPhotons_dR; vector<float> allfsrPhotons_pt; vector<float> allfsrPhotons_iso;

  // Event Category
  int EventCat;

  // -------------------------
  // GEN level information
  // -------------------------

  //Event variables
  int GENfinalState;

  // lepton variables
  vector<double> GENlep_pt; vector<double> GENlep_eta; vector<double> GENlep_phi; vector<double> GENlep_mass;
  vector<int> GENlep_id; vector<int> GENlep_status;
  vector<int> GENlep_MomId; vector<int> GENlep_MomMomId;
  int GENlep_Hindex[4];//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
  vector<float> GENlep_isoCH; vector<float> GENlep_isoNH; vector<float> GENlep_isoPhot; vector<float> GENlep_RelIso;

  // Z candidate variables
  vector<double> GENZ_pt; vector<double> GENZ_eta; vector<double> GENZ_phi; vector<double> GENZ_mass;
  vector<int> GENZ_DaughtersId; vector<int> GENZ_MomId;

  //hadZ and quark
  vector<double> GEN_Zq_pt, GEN_Zq_eta, GEN_Zq_phi, GEN_Zq_mass, GEN_Zq_MomMomid, GEN_Zq_Momid;
  vector<double> GEN_q_pt, GEN_q_eta, GEN_q_phi, GEN_q_mass;
  vector<int> GEN_Zq_id, GEN_q_id, GEN_q_status, GEN_q_Momid, GEN_q_MomMomid,GEN_q_nDaughters, GENH_isHard, GEN_id, GEN_status, GEN_Zq_isHard;
  vector<vector<int>> GEN_qdau_id, GEN_qdau_status;
  vector<vector<double>> GEN_qdau_pt, GEN_qdau_eta, GEN_qdau_phi, GEN_qdau_mass;
  // Higgs candidate variables (calculated using selected gen leptons)
  vector<double> GENH_pt; vector<double> GENH_eta; vector<double> GENH_phi; vector<double> GENH_mass;
  vector<double> GENH_Momid; vector<double> GENH_MomMomid; vector<double> GENH_status; vector<double> GENH_id;
  vector<double> GENH_nDaughters;
  vector<vector<double>> GENH_dau_pt, GENH_dau_eta,GENH_dau_phi,GENH_dau_mass;
  vector<vector<int>> GENH_dau_id, GENH_dau_status;
  float GENMH;

  //VBF quark
  vector<double> GEN_VBF_pt, GEN_VBF_eta, GEN_VBF_phi, GEN_VBF_mass;
  vector<int> GEN_VBF_id, GEN_VBF_Momid, GEN_VBF_MomMomid, GEN_VBF_status;

  // Jets
  vector<double> GENjet_pt; vector<double> GENjet_eta; vector<double> GENjet_phi; vector<double> GENjet_mass;
  int GENnjets_pt30_eta4p7; float GENpt_leadingjet_pt30_eta4p7;
  int GENnbjets_pt30_eta4p7;
  int GENnjets_pt30_eta2p5; float GENpt_leadingjet_pt30_eta2p5;
  float GENabsrapidity_leadingjet_pt30_eta4p7; float GENabsdeltarapidity_hleadingjet_pt30_eta4p7;
  int lheNb, lheNj, nGenStatus2bHad;

  // Global Variables but not stored in the tree
  vector<double> lep_ptreco;
  vector<int> lep_ptid; vector<int> lep_ptindex;
  vector<pat::Muon> recoMuons; vector<pat::Electron> recoElectrons; vector<pat::Electron> recoElectronsUnS;
  vector<pat::Tau> recoTaus; vector<pat::Photon> recoPhotons;
  vector<pat::PFParticle> fsrPhotons;
  TLorentzVector HVec, HVecNoFSR, Z1Vec, Z2Vec;
  TLorentzVector GENZ1Vec, GENZ2Vec;
  bool RecoFourMuEvent, RecoFourEEvent, RecoTwoETwoMuEvent, RecoTwoMuTwoEEvent;
  bool foundHiggsCandidate; bool foundZ1LCandidate; bool firstEntry;
  float jet1pt, jet2pt;
  float jet1pt2p5, jet2pt2p5;
  TLorentzVector Jet1, Jet2, Jet1_2p5, Jet2_2p5;

  vector<int> GENjet_id;

  // hist container
  std::map<std::string,TH1F*> histContainer_;

  //Input edm
  edm::EDGetTokenT<edm::View<pat::Electron> > elecSrc_;
  edm::EDGetTokenT<edm::View<pat::Electron> > elecUnSSrc_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonSrc_;
  edm::EDGetTokenT<edm::View<pat::Tau> > tauSrc_;
  edm::EDGetTokenT<edm::View<pat::Photon> > photonSrc_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetSrc_;
  edm::EDGetTokenT<edm::ValueMap<float> > qgTagSrc_;
  edm::EDGetTokenT<edm::ValueMap<float> > axis2Src_;
  edm::EDGetTokenT<edm::ValueMap<int> > multSrc_;
  edm::EDGetTokenT<edm::ValueMap<float> > ptDSrc_;
  edm::EDGetTokenT<edm::View<pat::Jet> > mergedjetSrc_;
  edm::EDGetTokenT<edm::View<pat::MET> > metSrc_;
  //edm::InputTag triggerSrc_;
  edm::EDGetTokenT<edm::TriggerResults> triggerSrc_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotSrc_;
  edm::EDGetTokenT<std::vector<reco::Conversion> > conversionSrc_;
  edm::EDGetTokenT<double> muRhoSrc_;
  edm::EDGetTokenT<double> elRhoSrc_;
  edm::EDGetTokenT<double> rhoSrcSUS_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSrc_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandsSrc_;
  edm::EDGetTokenT<edm::View<pat::PFParticle> > fsrPhotonsSrc_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedgenParticlesSrc_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedgenParticlesSrc_;
  edm::EDGetTokenT<edm::View<reco::GenJet> > genJetsSrc_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorSrc_;
  edm::EDGetTokenT<LHEEventProduct> lheInfoSrc_;
  edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken_;
  edm::EDGetTokenT<HTXS::HiggsClassification> htxsSrc_;
  //edm::EDGetTokenT<HZZFid::FiducialSummary> fidRivetSrc_;
  edm::EDGetTokenT< double > prefweight_token_;
  edm::EDGetTokenT< double > prefweightECAL_token_;
  edm::EDGetTokenT< double > prefweightMuon_token_;

  // Configuration
  const float Zmass;
  float mZ1Low, mZ2Low, mZ1High, mZ2High, m4lLowCut;
  float jetpt_cut, jeteta_cut;
  std::string elecID;
  bool isMC, isSignal;
  float mH;
  float crossSection;
  bool weightEvents;
  float isoCutEl, isoCutMu;
  double isoConeSizeEl, isoConeSizeMu;
  float sip3dCut, leadingPtCut, subleadingPtCut;
  float genIsoCutEl, genIsoCutMu;
  double genIsoConeSizeEl, genIsoConeSizeMu;
  float _elecPtCut, _muPtCut, _tauPtCut, _phoPtCut;
  float BTagCut;
  bool reweightForPU;
  std::string PUVersion;
  bool doFsrRecovery,bestCandMela, doMela, GENbestM4l, GENdoMela;
  //    bool GENdoMela;
  bool doPUJetID;
  int jetIDLevel;
  bool doJER;
  bool doJEC;
  bool doRefit;
  bool doTriggerMatching;
  bool checkOnlySingle;
  std::vector<std::string> triggerList;
  int skimLooseLeptons, skimTightLeptons;
  bool verbose;

  int year_,year;///use to choose Muon BDT
  bool isCode4l;

  // register to the TFileService
  edm::Service<TFileService> fs;

  // Counters
  float nEventsTotal;
  float sumWeightsTotal;
  float sumWeightsTotalPU;

  // JER
  JME::JetResolution resolution_pt, resolution_phi;
  JME::JetResolutionScaleFactor resolution_sf;

  JME::JetResolution ak8_resolution_pt, ak8_resolution_phi;
  JME::JetResolutionScaleFactor ak8_resolution_sf;

  string EleBDT_name_161718;
  string heepID_name_161718;

  // ME lists
  vector<string> lheMElist;
  vector<string> recoMElist;
  IvyMELAHelpers::GMECBlock MEblock;

};

UFHZZ4LAna::UFHZZ4LAna(const edm::ParameterSet& iConfig):
  histContainer_(),
  elecSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))),
  elecUnSSrc_(consumes<edm::View<pat::Electron> >(iConfig.getUntrackedParameter<edm::InputTag>("electronUnSSrc"))),
  muonSrc_(consumes<edm::View<pat::Muon> >(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc"))),
  tauSrc_(consumes<edm::View<pat::Tau> >(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc"))),
  photonSrc_(consumes<edm::View<pat::Photon> >(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc"))),
  jetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc"))),
  qgTagSrc_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
  axis2Src_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"))),
  multSrc_(consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult"))),
  ptDSrc_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"))),
  mergedjetSrc_(consumes<edm::View<pat::Jet> >(iConfig.getUntrackedParameter<edm::InputTag>("mergedjetSrc"))),
  metSrc_(consumes<edm::View<pat::MET> >(iConfig.getUntrackedParameter<edm::InputTag>("metSrc"))),
  triggerSrc_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerSrc"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  vertexSrc_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc"))),
  beamSpotSrc_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc"))),
  conversionSrc_(consumes<std::vector<reco::Conversion> >(iConfig.getUntrackedParameter<edm::InputTag>("conversionSrc"))),
  muRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("muRhoSrc"))),
  elRhoSrc_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("elRhoSrc"))),
  rhoSrcSUS_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("rhoSrcSUS"))),
  pileupSrc_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupSrc"))),
  pfCandsSrc_(consumes<pat::PackedCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCandsSrc"))),
  fsrPhotonsSrc_(consumes<edm::View<pat::PFParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("fsrPhotonsSrc"))),
  prunedgenParticlesSrc_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("prunedgenParticlesSrc"))),
  packedgenParticlesSrc_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("packedgenParticlesSrc"))),
  genJetsSrc_(consumes<edm::View<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genJetsSrc"))),
  generatorSrc_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("generatorSrc"))),
  lheInfoSrc_(consumes<LHEEventProduct>(iConfig.getUntrackedParameter<edm::InputTag>("lheInfoSrc"))),
  lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer",""))),
  htxsSrc_(consumes<HTXS::HiggsClassification>(edm::InputTag("rivetProducerHTXS","HiggsClassification"))),
  prefweight_token_(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"))),
  prefweightECAL_token_(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbECAL"))),
  prefweightMuon_token_(consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbMuon"))),
  //fidRivetSrc_(consumes<HZZFid::FiducialSummary>(edm::InputTag("rivetProducerHZZFid","FiducialSummary"))),
  Zmass(91.1876),
  mZ1Low(iConfig.getUntrackedParameter<double>("mZ1Low",40.0)),
  mZ2Low(iConfig.getUntrackedParameter<double>("mZ2Low",12.0)), // was 12
  mZ1High(iConfig.getUntrackedParameter<double>("mZ1High",120.0)),
  mZ2High(iConfig.getUntrackedParameter<double>("mZ2High",120.0)),
  m4lLowCut(iConfig.getUntrackedParameter<double>("m4lLowCut",70.0)),
  //     m4lLowCut(iConfig.getUntrackedParameter<double>("m4lLowCut",0.0)),
  jetpt_cut(iConfig.getUntrackedParameter<double>("jetpt_cut",10.0)),
  jeteta_cut(iConfig.getUntrackedParameter<double>("eta_cut",4.7)),
  elecID(iConfig.getUntrackedParameter<std::string>("elecID","NonTrig")),
  isMC(iConfig.getUntrackedParameter<bool>("isMC",true)),
  isSignal(iConfig.getUntrackedParameter<bool>("isSignal",false)),
  mH(iConfig.getUntrackedParameter<double>("mH",0.0)),
  crossSection(iConfig.getUntrackedParameter<double>("CrossSection",1.0)),
  weightEvents(iConfig.getUntrackedParameter<bool>("weightEvents",false)),
  isoCutEl(iConfig.getUntrackedParameter<double>("isoCutEl",9999.0)),
  isoCutMu(iConfig.getUntrackedParameter<double>("isoCutMu",0.35)),/////ios is applied to new Muon BDT //previous 0.35///Qianying
  isoConeSizeEl(iConfig.getUntrackedParameter<double>("isoConeSizeEl",0.3)),
  isoConeSizeMu(iConfig.getUntrackedParameter<double>("isoConeSizeMu",0.3)),
  sip3dCut(iConfig.getUntrackedParameter<double>("sip3dCut",4)),
  leadingPtCut(iConfig.getUntrackedParameter<double>("leadingPtCut",20.0)),
  subleadingPtCut(iConfig.getUntrackedParameter<double>("subleadingPtCut",10.0)),
  genIsoCutEl(iConfig.getUntrackedParameter<double>("genIsoCutEl",0.35)),
  genIsoCutMu(iConfig.getUntrackedParameter<double>("genIsoCutMu",0.35)),
  genIsoConeSizeEl(iConfig.getUntrackedParameter<double>("genIsoConeSizeEl",0.3)),
  genIsoConeSizeMu(iConfig.getUntrackedParameter<double>("genIsoConeSizeMu",0.3)),
  _elecPtCut(iConfig.getUntrackedParameter<double>("_elecPtCut",7.0)),
  _muPtCut(iConfig.getUntrackedParameter<double>("_muPtCut",5.0)),
  _tauPtCut(iConfig.getUntrackedParameter<double>("_tauPtCut",20.0)),
  _phoPtCut(iConfig.getUntrackedParameter<double>("_phoPtCut",10.0)),
  //BTagCut(iConfig.getUntrackedParameter<double>("BTagCut",0.4184)),/////2016: 0.6321; 2017: 0.4941; 2018: 0.4184
  reweightForPU(iConfig.getUntrackedParameter<bool>("reweightForPU",true)),
  PUVersion(iConfig.getUntrackedParameter<std::string>("PUVersion","Summer16_80X")),
  doFsrRecovery(iConfig.getUntrackedParameter<bool>("doFsrRecovery",true)),
  bestCandMela(iConfig.getUntrackedParameter<bool>("bestCandMela",true)),
  doMela(iConfig.getUntrackedParameter<bool>("doMela",false)),
  GENbestM4l(iConfig.getUntrackedParameter<bool>("GENbestM4l",false)),
  GENdoMela(iConfig.getUntrackedParameter<bool>("GENdoMela",true)),  // for GEN branches of mela based vairables
  doPUJetID(iConfig.getUntrackedParameter<bool>("doPUJetID",true)),
  jetIDLevel(iConfig.getUntrackedParameter<int>("jetIDLevel",2)),
  doJER(iConfig.getUntrackedParameter<bool>("doJER",true)),
  doJEC(iConfig.getUntrackedParameter<bool>("doJEC",true)),
  doRefit(iConfig.getUntrackedParameter<bool>("doRefit",true)),
  doTriggerMatching(iConfig.getUntrackedParameter<bool>("doTriggerMatching",!isMC)),
  checkOnlySingle(iConfig.getUntrackedParameter<bool>("checkOnlySingle",false)),
  triggerList(iConfig.getUntrackedParameter<std::vector<std::string>>("triggerList")),
  skimLooseLeptons(iConfig.getUntrackedParameter<int>("skimLooseLeptons",2)),
  skimTightLeptons(iConfig.getUntrackedParameter<int>("skimTightLeptons",2)),
  verbose(iConfig.getUntrackedParameter<bool>("verbose",false)),
  year_(iConfig.getUntrackedParameter<int>("year",2018)),////for year put 2016,2017, or 2018 to select correct training
  isCode4l(iConfig.getUntrackedParameter<bool>("isCode4l",false))
{
  if(!isMC){reweightForPU = false;}
  if(year_==2017||year_==2018||year_==2016)
  {
      year=year_;
  }
  if(year_==-2016)//-2016 -> pre VFP; 2016 -> post VFP
  {
      year=2016;
  }

  nEventsTotal=0.0;
  sumWeightsTotal=0.0;
  sumWeightsTotalPU=0.0;
  histContainer_["NEVENTS"]=fs->make<TH1F>("nEvents","nEvents in Sample",2,0,2);
  histContainer_["SUMWEIGHTS"]=fs->make<TH1F>("sumWeights","sum Weights of Sample",2,0,2);
  histContainer_["SUMWEIGHTSPU"]=fs->make<TH1F>("sumWeightsPU","sum Weights and PU of Sample",2,0,2);
  histContainer_["NVTX"]=fs->make<TH1F>("nVtx","Number of Vertices",36,-0.5,35.5);
  histContainer_["NVTX_RW"]=fs->make<TH1F>("nVtx_ReWeighted","Number of Vertices",36,-0.5,35.5);
  histContainer_["NINTERACT"]=fs->make<TH1F>("nInteractions","Number of True Interactions",61,-0.5,60.5);
  histContainer_["NINTERACT_RW"]=fs->make<TH1F>("nInteraction_ReWeighted","Number of True Interactions",61,-0.5,60.5);

  passedEventsTree_All = new TTree("passedEvents","passedEvents");

  //string elec_scalefac_Cracks_name_161718[3] = {"ElectronSF_Legacy_2016_Gap.root", "ElectronSF_Legacy_2017_Gap.root", "ElectronSF_Legacy_2018_Gap.root"};
  string elec_scalefac_Cracks_name_161718[3] = {"ElectronSF_UL2016postVFP_gap.root", "ElectronSF_UL2017_gap.root", "ElectronSF_UL2018_gap.root"};
  if (year_ == -2016) elec_scalefac_Cracks_name_161718[0]="ElectronSF_UL2016preVFP_gap.root";
  edm::FileInPath elec_scalefacFileInPathCracks(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_scalefac_Cracks_name_161718[year-2016]).c_str());
  TFile *fElecScalFacCracks = TFile::Open(elec_scalefacFileInPathCracks.fullPath().c_str());
  hElecScaleFac_Cracks = (TH2F*)fElecScalFacCracks->Get("EGamma_SF2D");

  //string elec_scalefac_name_161718[3] = {"ElectronSF_Legacy_2016_NoGap.root", "ElectronSF_Legacy_2017_NoGap.root", "ElectronSF_Legacy_2018_NoGap.root"};
  string elec_scalefac_name_161718[3] = {"ElectronSF_UL2016postVFP_nogap.root", "ElectronSF_UL2017_nogap.root", "ElectronSF_UL2018_nogap.root"};
  if (year_ == -2016)
      elec_scalefac_name_161718[0]="ElectronSF_UL2016preVFP_nogap.root";
  edm::FileInPath elec_scalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_scalefac_name_161718[year-2016]).c_str());
  TFile *fElecScalFac = TFile::Open(elec_scalefacFileInPath.fullPath().c_str());
  hElecScaleFac = (TH2F*)fElecScalFac->Get("EGamma_SF2D");

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaUL2016To2018#A_note_about_IDs_in_UL
  string elec_Gsfscalefac_name_161718[3] = {"egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root", "egammaEffi_ptAbove20.txt_EGM2D_UL2017.root", "egammaEffi_ptAbove20.txt_EGM2D_UL2018.root"};
  if (year_ == -2016)
      elec_Gsfscalefac_name_161718[0]="egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root";
  edm::FileInPath elec_GsfscalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_Gsfscalefac_name_161718[year-2016]).c_str());
  TFile *fElecScalFacGsf = TFile::Open(elec_GsfscalefacFileInPath.fullPath().c_str());
  hElecScaleFacGsf = (TH2F*)fElecScalFacGsf->Get("EGamma_SF2D");

  string elec_GsfLowETscalefac_name_161718[3]= {"egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root", "egammaEffi_ptBelow20.txt_EGM2D_UL2017.root", "egammaEffi_ptBelow20.txt_EGM2D_UL2018.root"};
  if (year_ == -2016)
      elec_GsfLowETscalefac_name_161718[0]="egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root";
  edm::FileInPath elec_GsfLowETscalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+elec_GsfLowETscalefac_name_161718[year-2016]).c_str());
  TFile *fElecScalFacGsfLowET = TFile::Open(elec_GsfLowETscalefacFileInPath.fullPath().c_str());
  hElecScaleFacGsfLowET = (TH2F*)fElecScalFacGsfLowET->Get("EGamma_SF2D");

  string mu_scalefac_name_161718[3] = {"final_HZZ_muon_SF_2016RunB2H_legacy_newLoose_newIso_paper.root", "final_HZZ_muon_SF_2017_newLooseIso_mupogSysts_paper.root", "final_HZZ_muon_SF_2018RunA2D_ER_newLoose_newIso_paper.root"};
  //string mu_scalefac_name_161718[3] = {"final_HZZ_SF_2016UL_mupogsysts_newLoose_pairdR0pt1_06072021.root", "final_HZZ_2017UL_SF_noTracking_Loose_syst_mupogSysts_pairdR0pt1_06072021.root", "final_HZZ_muon_SF_2018UL_POG_newLoose_JG_b4_JDG_b3_newIso_Syst_06072021.root"};
  edm::FileInPath mu_scalefacFileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+mu_scalefac_name_161718[year-2016]).c_str());
  TFile *fMuScalFac = TFile::Open(mu_scalefacFileInPath.fullPath().c_str());
  hMuScaleFac = (TH2F*)fMuScalFac->Get("FINAL");
  hMuScaleFacUnc = (TH2F*)fMuScalFac->Get("ERROR");

  string pileup_name_161718[3] = {"pileup_UL_2016.root", "pileup_UL_2017.root", "pileup_UL_2018.root"};
  edm::FileInPath pileup_FileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+pileup_name_161718[year-2016]).c_str());
  TFile *f_pileup = TFile::Open(pileup_FileInPath.fullPath().c_str());
  h_pileup = (TH1D*)f_pileup->Get("weights");
  h_pileupUp = (TH1D*)f_pileup->Get("weights_varUp");
  h_pileupDn = (TH1D*)f_pileup->Get("weights_varDn");

  string bTagEffi_name_161718[3] = {"bTagEfficiencies_2016.root", "bTagEfficiencies_2017.root", "bTagEfficiencies_2018.root"};
  edm::FileInPath BTagEffiInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+bTagEffi_name_161718[year-2016]).c_str());
  TFile *fbTagEffi = TFile::Open(BTagEffiInPath.fullPath().c_str());
  hbTagEffi = (TH2F*)fbTagEffi->Get("eff_b_M_ALL");
  hcTagEffi = (TH2F*)fbTagEffi->Get("eff_c_M_ALL");
  hudsgTagEffi = (TH2F*)fbTagEffi->Get("eff_udsg_M_ALL");

  //BTag calibration
  //string csv_name_161718[3] = {"DeepCSV_2016LegacySF_V1.csv", "DeepCSV_106XUL17SF_V2p1.csv", "DeepCSV_106XUL18SF.csv"};
  //string csv_name_161718[3] = {"DeepCSV_106XUL16postVFPSF_v2_v2.csv", "DeepCSV_106XUL17SF_V2p1.csv", "DeepCSV_106XUL18SF_V1p1.csv"};
  string csv_name_161718[3] = {"DeepCSV_106XUL16postVFPSF_v2_hzz.csv", "wp_deepCSV_106XUL17_v3_hzz.csv", "wp_deepCSV_106XUL18_v2_hzz.csv"};
  if (year_ == -2016)
      //csv_name_161718[0]="DeepCSV_106XUL16preVFPSF_v1_v1.csv";
      csv_name_161718[0]="DeepCSV_106XUL16preVFPSF_v1_hzz.csv";
  if (year==2016 && verbose)     std::cout<<"BTag calibration csv: "<<csv_name_161718[year-2016]<<std::endl;

  edm::FileInPath btagfileInPath(("UFHZZAnalysisRun2/UFHZZ4LAna/data/"+csv_name_161718[year-2016]).c_str());

  BTagCalibration calib("DeepCSV", btagfileInPath.fullPath().c_str());
  reader = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,  // operating point
                                     "central",             // central sys type
                                     {"up", "down"});      // other sys types
  // Btag information
  // https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation
  reader->load(calib,                // calibration instance
               BTagEntry::FLAV_B,    // btag flavour
               "comb");               // measurement type

  //if(year==2018)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Autumn18IdIsoValues"; BTagCut=0.4184; heepID_name_161718 = "heepElectronID-HEEPV70";}
  if(year==2018)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Summer18ULIdIsoValues"; BTagCut=0.4168; heepID_name_161718 = "heepElectronID-HEEPV70";}
  //if(year==2017)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Fall17IsoV2Values"; BTagCut=0.4941; heepID_name_161718 = "heepElectronID-HEEPV70";}
  if(year==2017)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Summer17ULIdIsoValues"; BTagCut=0.4506; heepID_name_161718 = "heepElectronID-HEEPV70";}
  //if(year==2016)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Summer16IdIsoValues"; BTagCut=0.6321; heepID_name_161718 = "heepElectronID-HEEPV70";}
  if(year_==2016)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Summer16ULIdIsoValues"; BTagCut=0.5847; heepID_name_161718 = "heepElectronID-HEEPV70";}
  if(year_==-2016)    {EleBDT_name_161718 = "ElectronMVAEstimatorRun2Summer16ULIdIsoValues"; BTagCut=0.6001; heepID_name_161718 = "heepElectronID-HEEPV70";}

  std::string DATAPATH = std::getenv( "CMSSW_BASE" );
  //std::string DATAPATH = std::getenv( "CMSSW_BASE" );
  //if(year == 2018)    DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2018.txt";
  if(year == 2018)    DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v5/RoccoR2018UL.txt";
  //if(year == 2017)    DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2017.txt";
  if(year == 2017)    DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v5/RoccoR2017UL.txt";
  //if(year == 2016)    DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v3/RoccoR2016.txt";
  if(year_ == 2016)    DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v5/RoccoR2016bUL.txt"; // for post VFP
  if(year_ == -2016)    DATAPATH+="/src/UFHZZAnalysisRun2/KalmanMuonCalibrationsProducer/data/roccor.Run2.v5/RoccoR2016aUL.txt"; // for pre VFP
  calibrator = new RoccoR(DATAPATH);


}

UFHZZ4LAna::~UFHZZ4LAna()
{
    //destructor --- don't do anything here
}


float KalmanEnergy(float px, float py, float pz, float mass){

    double E=TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);

    return E;
}



// ------------ method called for each event  ------------
void UFHZZ4LAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  using namespace pat;
  using namespace trigger;
  using namespace EwkCorrections;

  nEventsTotal += 1.0;

  Run = iEvent.id().run();
  Event = iEvent.id().event();
  LumiSect = iEvent.id().luminosityBlock();

  if (verbose) {
      cout<<"Run: " << Run << ",Event: " << Event << ",LumiSect: "<<LumiSect<<endl;
  }

  // ======= Get Collections ======= //
  if (verbose) {cout<<"getting collections"<<endl;}

  // trigger collection
  edm::Handle<edm::TriggerResults> trigger;
  iEvent.getByToken(triggerSrc_,trigger);
  const edm::TriggerNames trigNames = iEvent.triggerNames(*trigger);

  // trigger Objects
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);

  // vertex collection
  edm::Handle<reco::VertexCollection> vertex;
  iEvent.getByToken(vertexSrc_,vertex);

  // electron collection
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByToken(elecSrc_,electrons);
  if (verbose) cout<<electrons->size()<<" total electrons in the collection"<<endl;

  // electron before scale/smearing corrections
  edm::Handle<edm::View<pat::Electron> > electronsUnS;
  iEvent.getByToken(elecUnSSrc_,electronsUnS);

  // muon collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonSrc_,muons);
  if (verbose) cout<<muons->size()<<" total muons in the collection"<<endl;

  // photon collection
  edm::Handle<edm::View<pat::Photon> > photons;
  iEvent.getByToken(photonSrc_,photons);
  if (verbose) cout<<photons->size()<<" total photons in the collection"<<endl;

  // met collection
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByToken(metSrc_,mets);

  // Rho Correction
  edm::Handle<double> eventRhoMu;
  iEvent.getByToken(muRhoSrc_,eventRhoMu);
  muRho = *eventRhoMu;

  edm::Handle<double> eventRhoE;
  iEvent.getByToken(elRhoSrc_,eventRhoE);
  elRho = *eventRhoE;

  edm::Handle<double> eventRhoSUS;
  iEvent.getByToken(rhoSrcSUS_,eventRhoSUS);
  rhoSUS = *eventRhoSUS;

  // Conversions
  edm::Handle< std::vector<reco::Conversion> > theConversions;
  iEvent.getByToken(conversionSrc_, theConversions);

  // Beam Spot
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(beamSpotSrc_,beamSpot);
  const reco::BeamSpot BS = *beamSpot;

  // Particle Flow Cands
  edm::Handle<pat::PackedCandidateCollection> pfCands;
  iEvent.getByToken(pfCandsSrc_,pfCands);

  // FSR Photons
  edm::Handle<edm::View<pat::PFParticle> > photonsForFsr;
  iEvent.getByToken(fsrPhotonsSrc_,photonsForFsr);

  // Jets
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetSrc_,jets);

  if (!jecunc) {
      //ak4 uncertainty
      edm::ESHandle<JetCorrectorParametersCollection> jetCorrParameterSet;
      iSetup.get<JetCorrectionsRecord>().get("AK4PFchs", jetCorrParameterSet);
      const JetCorrectorParameters& jetCorrParameters = (*jetCorrParameterSet)["Uncertainty"];
      jecunc.reset(new JetCorrectionUncertainty(jetCorrParameters));

      //AK8 uncertainty
      iSetup.get<JetCorrectionsRecord>().get("AK8PFPuppi", jetCorrParameterSet);
      const JetCorrectorParameters& jetCorrParameters = (*jetCorrParameterSet)["Uncertainty"];
      jecmergedunc.reset(new JetCorrectionUncertainty(jetCorrParameters));
  }

  resolution_pt = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
  resolution_phi = JME::JetResolution::get(iSetup, "AK4PFchs_phi");
  resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");

  ak8_resolution_pt = JME::JetResolution::get(iSetup, "AK8PFPuppi_pt");
  ak8_resolution_phi = JME::JetResolution::get(iSetup, "AK8PFPuppi_phi");
  ak8_resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK8PFPuppi");

  edm::Handle<edm::ValueMap<float> > qgHandle;
  iEvent.getByToken(qgTagSrc_, qgHandle);

  edm::Handle<edm::ValueMap<float> > axis2Handle;
  iEvent.getByToken(axis2Src_, axis2Handle);

  edm::Handle<edm::ValueMap<int> > multHandle;
  iEvent.getByToken(multSrc_, multHandle);

  edm::Handle<edm::ValueMap<float> > ptDHandle;
  iEvent.getByToken(ptDSrc_, ptDHandle);

  edm::Handle<edm::View<pat::Jet> > mergedjets;
  iEvent.getByToken(mergedjetSrc_,mergedjets);

  // GEN collections
  edm::Handle<reco::GenParticleCollection> prunedgenParticles;
  iEvent.getByToken(prunedgenParticlesSrc_, prunedgenParticles);

  edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles;
  iEvent.getByToken(packedgenParticlesSrc_, packedgenParticles);

  edm::Handle<edm::View<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsSrc_, genJets);

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(generatorSrc_,genEventInfo);

  //vector<edm::Handle<LHEEventProduct> > lheInfos;
  //iEvent.getManyByType(lheInfos); // using this method because the label is not always the same (e.g. "source" in the ttH sample)

  edm::Handle<LHEEventProduct> lheInfo;
  iEvent.getByToken(lheInfoSrc_, lheInfo);

  if (isMC) {
      edm::Handle< double > theprefweight;
      iEvent.getByToken(prefweight_token_, theprefweight ) ;

      edm::Handle< double > theprefweightECAL;
      iEvent.getByToken(prefweightECAL_token_, theprefweightECAL ) ;

      edm::Handle< double > theprefweightMuon;
      iEvent.getByToken(prefweightMuon_token_, theprefweightMuon ) ;

      prefiringWeight = (*theprefweight);
      prefiringWeightECAL = (*theprefweightECAL);
      prefiringWeightMuon = (*theprefweightMuon);

  }
  else {
      prefiringWeight = 1.0;
      prefiringWeightECAL = 1.0;
      prefiringWeightMuon = 1.0;
  }

  // ============ Initialize Variables ============= //
  // Event Variables
  if (verbose) {cout<<"clear variables"<<endl;}
  nVtx = -1.0; nInt = -1.0;
  finalState = -1;
  triggersPassed="";
  passedTrig=false;

  // Event Weights
  genWeight=1.0; pileupWeight=1.0; pileupWeightUp=1.0; pileupWeightDn=1.0; dataMCWeight=1.0; eventWeight=1.0;

  qcdWeights.clear(); nnloWeights.clear(); pdfWeights.clear();
  pdfRMSup=1.0; pdfRMSdown=1.0; pdfENVup=1.0; pdfENVdown=1.0;

  //lepton variables
  lep_d0BS.clear();
  lep_d0PV.clear();
  lep_numberOfValidPixelHits.clear();
  lep_trackerLayersWithMeasurement.clear();

  lep_pt_FromMuonBestTrack.clear(); lep_eta_FromMuonBestTrack.clear(); lep_phi_FromMuonBestTrack.clear();
  lep_position_x.clear();	lep_position_y.clear();	lep_position_z.clear();
  lep_pt_genFromReco.clear();
  lep_pt_UnS.clear(); lep_pterrold_UnS.clear();
  lep_pt.clear(); lep_pterr.clear(); lep_pterrold.clear();
  lep_p.clear(); lep_ecalEnergy.clear(); lep_isEB.clear(); lep_isEE.clear();
  lep_errPre_Scale.clear(); lep_errPost_Scale.clear(); lep_errPre_noScale.clear(); lep_errPost_noScale.clear();
  lep_eta.clear(); lep_phi.clear(); lep_mass.clear();
  lepFSR_pt.clear(); lepFSR_eta.clear(); lepFSR_phi.clear(); lepFSR_mass.clear(); lepFSR_ID.clear();

  lep_genindex.clear(); lep_id.clear(); lep_dataMC.clear(); lep_dataMCErr.clear();
  dataMC_VxBS.clear(); dataMCErr_VxBS.clear();
  lep_matchedR03_PdgId.clear(); lep_matchedR03_MomId.clear(); lep_matchedR03_MomMomId.clear();
  lep_mva.clear(); lep_ecalDriven.clear();
  lep_tightId.clear(); lep_tightIdSUS.clear(); lep_tightIdHiPt.clear(); //lep_tightId_old.clear();
  lep_Sip.clear(); lep_IP.clear();
  lep_isoNH.clear(); lep_isoCH.clear(); lep_isoPhot.clear(); lep_isoPU.clear(); lep_isoPUcorr.clear();
  lep_RelIso.clear(); lep_RelIsoNoFSR.clear(); lep_MiniIso.clear();
  lep_ptRatio.clear(); lep_ptRel.clear();
  lep_missingHits.clear();
  lep_filtersMatched.clear();
  nisoleptons=0;

  singleBS_RecoLep_pt.clear();
  singleBS_RecoLep_ptError.clear();
  singleBS_RecoLep_eta.clear();
  singleBS_RecoLep_phi.clear();
  singleBS_RecoLep_mass.clear();
  singleBS_RecoLep_d0.clear();

  // photon variables
  pho_pt.clear(); pho_eta.clear(); pho_phi.clear(); photonCutBasedIDLoose.clear();

  // MET
  met=-1.0; met_phi=9999.0;
  met_jesup=-1.0; met_phi_jesup=9999.0; met_jesdn=-1.0; met_phi_jesdn=9999.0;
  met_uncenup=-1.0; met_phi_uncenup=9999.0; met_uncendn=-1.0; met_phi_uncendn=9999.0;

  // Jets
  jet_pt.clear(); jet_eta.clear(); jet_phi.clear(); jet_mass.clear(); jet_pt_raw.clear();
  jet_jesup_pt.clear(); jet_jesup_eta.clear(); jet_jesup_phi.clear(); jet_jesup_mass.clear();
  jet_jesdn_pt.clear(); jet_jesdn_eta.clear(); jet_jesdn_phi.clear(); jet_jesdn_mass.clear();
  jet_jerup_pt.clear(); jet_jerup_eta.clear(); jet_jerup_phi.clear(); jet_jerup_mass.clear();
  jet_jerdn_pt.clear(); jet_jerdn_eta.clear(); jet_jerdn_phi.clear(); jet_jerdn_mass.clear();
  jet_csvv2_.clear();
  jet_pumva.clear(); jet_csvv2.clear(); jet_isbtag.clear();
  jet_hadronFlavour.clear(); jet_partonFlavour.clear();
  jet_QGTagger.clear(); jet_QGTagger_jesup.clear(); jet_QGTagger_jesdn.clear();
  jet_relpterr.clear(); jet_phierr.clear();
  jet_bTagEffi.clear();
  jet_cTagEffi.clear();
  jet_udsgTagEffi.clear();
  jet_axis2.clear(); jet_ptD.clear(); jet_mult.clear();

  jet_iscleanH4l.clear();
  jet1index=-1; jet2index=-1; jet1index_2p5 = -1, jet2index_2p5 = -1;
  jet_jesup_iscleanH4l.clear(); jet_jesdn_iscleanH4l.clear();
  jet_jerup_iscleanH4l.clear(); jet_jerdn_iscleanH4l.clear();

  njets_pt30_eta4p7=0;
  njets_pt30_eta4p7_jesup=0; njets_pt30_eta4p7_jesdn=0;
  njets_pt30_eta4p7_jerup=0; njets_pt30_eta4p7_jerdn=0;

  njets_pt30_eta2p5=0;
  njets_pt30_eta2p5_jesup=0; njets_pt30_eta2p5_jesdn=0;
  njets_pt30_eta2p5_jerup=0; njets_pt30_eta2p5_jerdn=0;

  nbjets_pt30_eta4p7=0; nvjets_pt40_eta2p4=0;

  pt_leadingjet_pt30_eta4p7=-1.0;

  pt_leadingjet_pt30_eta4p7_jesup=-1.0; pt_leadingjet_pt30_eta4p7_jesdn=-1.0;
  pt_leadingjet_pt30_eta4p7_jerup=-1.0; pt_leadingjet_pt30_eta4p7_jerdn=-1.0;

  pt_leadingjet_pt30_eta2p5=-1.0;
  pt_leadingjet_pt30_eta2p5_jesup=-1.0; pt_leadingjet_pt30_eta2p5_jesdn=-1.0;
  pt_leadingjet_pt30_eta2p5_jerup=-1.0; pt_leadingjet_pt30_eta2p5_jerdn=-1.0;

  absrapidity_leadingjet_pt30_eta4p7=-1.0;
  absrapidity_leadingjet_pt30_eta4p7_jesup=-1.0; absrapidity_leadingjet_pt30_eta4p7_jesdn=-1.0;
  absrapidity_leadingjet_pt30_eta4p7_jerup=-1.0; absrapidity_leadingjet_pt30_eta4p7_jerdn=-1.0;

  absdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;
  absdeltarapidity_hleadingjet_pt30_eta4p7_jesup=-1.0; absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn=-1.0;
  absdeltarapidity_hleadingjet_pt30_eta4p7_jerup=-1.0; absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn=-1.0;

  DijetMass=-1.0; DijetDEta=9999.0; DijetFisher=9999.0;

  mergedjet_iscleanH4l.clear();
  mergedjet_pt.clear(); mergedjet_eta.clear(); mergedjet_phi.clear(); mergedjet_mass.clear();
  mergedjet_jer.clear();
  mergedjet_jer_pt.clear(); mergedjet_jer_eta.clear(); mergedjet_jer_phi.clear(); mergedjet_jer_mass.clear();
  mergedjet_jerup_pt.clear(); mergedjet_jerup_eta.clear(); mergedjet_jerup_phi.clear(); mergedjet_jerup_mass.clear();
  mergedjet_jerdn_pt.clear(); mergedjet_jerdn_eta.clear(); mergedjet_jerdn_phi.clear(); mergedjet_jerdn_mass.clear();
  mergedjet_ptrow.clear();
  mergedjet_jer_pterr.clear(); mergedjet_jer_phierr.clear();

  //mergedjet_L1.clear();
  mergedjet_softdropmass.clear();
  //mergedjet_prunedmass.clear();
  mergedjet_tau1.clear(); mergedjet_tau2.clear();
  mergedjet_btag.clear();

  mergedjet_ZvsQCD.clear();
  mergedjet_ZbbvsQCD.clear();
  mergedjet_ZHbbvsQCD.clear();
  mergedjet_WvsQCD.clear();
  mergedjet_ZHbbvsQCD.clear();
  mergedjet_HbbvsQCD.clear();
  mergedjet_H4qvsQCD.clear();

  mergedjet_ZvsQCD_de.clear();
  mergedjet_ZbbvsQCD_de.clear();
  mergedjet_ZHbbvsQCD_de.clear();
  mergedjet_WvsQCD_de.clear();
  mergedjet_ZHbbvsQCD_de.clear();
  mergedjet_HbbvsQCD_de.clear();
  mergedjet_H4qvsQCD_de.clear();

  mergedjet_Net_Xbb_de.clear();
  mergedjet_Net_Xcc_de.clear();
  mergedjet_Net_Xqq_de.clear();
  mergedjet_Net_QCDbb_de.clear();
  mergedjet_Net_QCDcc_de.clear();
  mergedjet_Net_QCDb_de.clear();
  mergedjet_Net_QCDc_de.clear();
  mergedjet_Net_QCDother_de.clear();

  mergedjet_nsubjet.clear();
  mergedjet_subjet_pt.clear(); mergedjet_subjet_eta.clear();
  mergedjet_subjet_phi.clear(); mergedjet_subjet_mass.clear();
  mergedjet_subjetUncorr_pt.clear(); mergedjet_subjetUncorr_eta.clear();
  mergedjet_subjetUncorr_phi.clear(); mergedjet_subjetUncorr_mass.clear(); mergedjet_subjet_softDropMassUncorr.clear();
  mergedjet_subjet_softDropMass.clear();
  mergedjet_subjet_btag.clear();
  mergedjet_subjet_partonFlavour.clear(); mergedjet_subjet_hadronFlavour.clear();
  mergedjet_nbHadrons.clear(); mergedjet_ncHadrons.clear(); mergedjet_hadronFlavour.clear(); mergedjet_partonFlavour.clear();

  // FSR Photons
  nFSRPhotons=0;
  fsrPhotons_lepindex.clear(); fsrPhotons_pt.clear(); fsrPhotons_pterr.clear();
  fsrPhotons_eta.clear(); fsrPhotons_phi.clear();
  fsrPhotons_dR.clear(); fsrPhotons_iso.clear();
  allfsrPhotons_dR.clear(); allfsrPhotons_pt.clear(); allfsrPhotons_iso.clear();

  // -------------------------
  // GEN level information
  // -------------------------

  //Event variables
  GENfinalState=-1;

  // lepton variables
  GENlep_pt.clear(); GENlep_eta.clear(); GENlep_phi.clear(); GENlep_mass.clear();
  GENlep_id.clear(); GENlep_status.clear(); GENlep_MomId.clear(); GENlep_MomMomId.clear();
  for (int i=0; i<4; ++i) {GENlep_Hindex[i]=-1;};//position of Higgs candidate leptons in lep_p4: 0 = Z1 lead, 1 = Z1 sub, 2 = Z2 lead, 3 = Z3 sub
  GENlep_isoCH.clear(); GENlep_isoNH.clear(); GENlep_isoPhot.clear(); GENlep_RelIso.clear();

  // Z candidate variables
  GENZ_pt.clear(); GENZ_eta.clear(); GENZ_phi.clear(); GENZ_mass.clear();
  GENZ_DaughtersId.clear(); GENZ_MomId.clear();

  //HadZ and quark
  GEN_Zq_pt.clear(); GEN_Zq_eta.clear(); GEN_Zq_phi.clear(); GEN_Zq_mass.clear(); GEN_Zq_id.clear(); GEN_Zq_Momid.clear(); GEN_Zq_MomMomid.clear();
  GEN_q_id.clear(); GEN_q_pt.clear(); GEN_q_eta.clear(); GEN_q_phi.clear(); GEN_q_mass.clear(); GEN_q_status.clear(); GEN_q_Momid.clear(); GEN_q_MomMomid.clear();
  GEN_qdau_id.clear(); GEN_qdau_pt.clear(); GEN_qdau_eta.clear(); GEN_qdau_phi.clear(); GEN_qdau_mass.clear(); GEN_qdau_status.clear(); GEN_q_nDaughters.clear();
  GEN_Zq_isHard.clear();

  // Higgs candidate variables (calculated using selected gen leptons)
  GENH_pt.clear(); GENH_eta.clear(); GENH_phi.clear(); GENH_mass.clear();
  GENMH=-1.0;
  GENH_Momid.clear(); GENH_MomMomid.clear();
  GENH_status.clear(); GENH_id.clear(); GENH_isHard.clear();
  GEN_id.clear(); GEN_status.clear();
  GENH_dau_id.clear(); GENH_nDaughters.clear(); GENH_dau_status.clear(); GENH_dau_pt.clear(); GENH_dau_eta.clear(); GENH_dau_phi.clear(); GENH_dau_mass.clear();

  //VBF quark
  GEN_VBF_pt.clear(); GEN_VBF_eta.clear(); GEN_VBF_phi.clear(); GEN_VBF_mass.clear(); GEN_VBF_id.clear(); GEN_VBF_Momid.clear(); GEN_VBF_MomMomid.clear(); GEN_VBF_status.clear();

  // Jets
  GENjet_pt.clear(); GENjet_eta.clear(); GENjet_phi.clear(); GENjet_mass.clear();
  GENnjets_pt30_eta4p7=0;
  GENnbjets_pt30_eta4p7=0;
  GENnjets_pt30_eta2p5=0;
  GENpt_leadingjet_pt30_eta4p7=-1.0; GENabsrapidity_leadingjet_pt30_eta4p7=-1.0; GENabsdeltarapidity_hleadingjet_pt30_eta4p7=-1.0;
  GENpt_leadingjet_pt30_eta2p5=-1.0;
  GENjet_id.clear();
  lheNb=0; lheNj=0; nGenStatus2bHad=0;

  if (verbose) {cout<<"clear other variables"<<endl; }
  // Resolution
  //massErrorUCSD=-1.0; massErrorUCSDCorr=-1.0; massErrorUF=-1.0; massErrorUFCorr=-1.0; massErrorUFADCorr=-1.0;

  // Event Category
  EventCat=-1;

  // Global variables not stored in tree
  lep_ptreco.clear(); lep_ptid.clear(); lep_ptindex.clear();
  recoMuons.clear(); recoElectrons.clear(); fsrPhotons.clear(); recoElectronsUnS.clear();
  HVec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  HVecNoFSR.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  Z1Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  Z2Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  GENZ1Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  GENZ2Vec.SetPtEtaPhiM(0.0,0.0,0.0,0.0);
  RecoFourMuEvent = false; RecoFourEEvent = false;
  RecoTwoETwoMuEvent = false; RecoTwoMuTwoEEvent = false;
  foundHiggsCandidate = false; foundZ1LCandidate = false;
  jet1pt=-1.0; jet2pt=-1.0;
  jet1pt2p5=-1.0; jet2pt2p5=-1.0;

  // ====================== Do Analysis ======================== //
  std::map<int, TLorentzVector> fsrmap;
  vector<reco::Candidate*> selectedLeptons;
  std::map<unsigned int, TLorentzVector> selectedFsrMap;

  fsrmap.clear(); selectedFsrMap.clear(); selectedLeptons.clear();

  if (verbose) cout<<"start pileup reweighting"<<endl;
  // PU information
  if(isMC && reweightForPU) {
      edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
      iEvent.getByToken(pileupSrc_, PupInfo);

      if (verbose) cout<<"got pileup info"<<endl;

      std::vector<PileupSummaryInfo>::const_iterator PVI;
      int npv = -1;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
          int BX = PVI->getBunchCrossing();
          if(BX == 0) { npv = PVI->getTrueNumInteractions(); continue;}
      }
      if (verbose) cout<<"N true interations = "<<npv<<endl;
      nInt = npv;
      //pileupWeight = pileUp.getPUWeight(npv,PUVersion);
      pileupWeight = pileUp.getPUWeight(h_pileup,npv);
      pileupWeightUp = pileUp.getPUWeight(h_pileupUp,npv);
      pileupWeightDn = pileUp.getPUWeight(h_pileupDn,npv);
      if (verbose) cout<<"pileup weight = "<<pileupWeight<<", filling histograms"<<endl;
      histContainer_["NINTERACT"]->Fill(npv);
      histContainer_["NINTERACT_RW"]->Fill(npv,pileupWeight);
  } else { pileupWeight = 1.0;}

  if (verbose) {cout<<"finished pileup reweighting"<<endl; }

  if(isMC){ //gen INFO
    float tmpWeight = genEventInfo->weight();
    genWeight = (tmpWeight > 0 ? 1.0 : -1.0);
    if (verbose) {cout<<"tmpWeight: "<<tmpWeight<<"; genWeight: "<<genWeight<<endl;}
    double rms = 0.0;

    if(lheInfo.isValid()){

      for(unsigned int i = 0; i < lheInfo->weights().size(); i++) {

          tmpWeight = genEventInfo->weight();
          tmpWeight *= lheInfo->weights()[i].wgt/lheInfo->originalXWGTUP();
          pdfWeights.push_back(tmpWeight);

          if (i<=8 or int(i)>=posNNPDF) {
              tmpWeight = genEventInfo->weight();
              tmpWeight *= lheInfo->weights()[i].wgt/lheInfo->originalXWGTUP();
              if (int(i)<posNNPDF) {qcdWeights.push_back(tmpWeight);}
          }
          else {
              tmpWeight = lheInfo->weights()[i].wgt;
              tmpWeight /= lheInfo->originalXWGTUP();
              //if (i==9) genWeight = tmpWeight;
              if (int(i)<posNNPDF) {nnloWeights.push_back(tmpWeight);}
          }
          // NNPDF30 variations
          if (int(i)>=posNNPDF && int(i)<=(posNNPDF+100)) {
              rms += tmpWeight*tmpWeight;
              if (tmpWeight>pdfENVup) pdfENVup=tmpWeight;
              if (tmpWeight<pdfENVdown) pdfENVdown=tmpWeight;
          }
      }

      pdfRMSup=sqrt(rms/100.0); pdfRMSdown=1.0/pdfRMSup;
      if (verbose) cout<<"pdfRMSup "<<pdfRMSup<<" pdfRMSdown "<<pdfRMSdown<<endl;

      const lhef::HEPEUP& lheEvent = lheInfo->hepeup();
      std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;

      for ( size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle ) {
          int id = std::abs(lheEvent.IDUP[idxParticle]);
          int status = lheEvent.ISTUP[idxParticle];
          if ( status == 1 && id==5 ) {
              lheNb += 1;
          }
          if ( status == 1 && ((id >= 1 && id <= 6) || id == 21) ) {
              lheNj += 1;
          }
      }
    }

    if (verbose) cout<<"setting gen variables"<<endl;
    setGENVariables(prunedgenParticles,packedgenParticles,genJets);
    if (verbose) {cout<<"finshed setting gen variables"<<endl;}

  }
  sumWeightsTotal += genWeight;
  sumWeightsTotalPU += pileupWeight*genWeight;
  eventWeight = pileupWeight*genWeight;

  unsigned int _tSize = trigger->size();
  // create a string with all passing trigger names
  for (unsigned int i=0; i<_tSize; ++i) {
      std::string triggerName = trigNames.triggerName(i);
      if (strstr(triggerName.c_str(),"_step")) continue;
      if (strstr(triggerName.c_str(),"MC_")) continue;
      if (strstr(triggerName.c_str(),"AlCa_")) continue;
      if (strstr(triggerName.c_str(),"DST_")) continue;
      if (strstr(triggerName.c_str(),"HLT_HI")) continue;
      if (strstr(triggerName.c_str(),"HLT_Physics")) continue;
      if (strstr(triggerName.c_str(),"HLT_Random")) continue;
      if (strstr(triggerName.c_str(),"HLT_ZeroBias")) continue;
      if (strstr(triggerName.c_str(),"HLT_IsoTrack")) continue;
      if (strstr(triggerName.c_str(),"Hcal")) continue;
      if (strstr(triggerName.c_str(),"Ecal")) continue;
      if (trigger->accept(i)) triggersPassed += triggerName;
  }
  if (firstEntry) cout<<"triggersPassed: "<<triggersPassed<<endl;
  firstEntry = false;
  // check if any of the triggers in the user list have passed
  bool passedSingleEl=false;
  bool passedSingleMu=false;
  bool passedAnyOther=false;
  for (unsigned int i=0; i<triggerList.size(); ++i) {
      if (strstr(triggersPassed.c_str(),triggerList.at(i).c_str())) {
          passedTrig=true;
          if (!isMC) {
              if (strstr(triggerList.at(i).c_str(),"_WP")) passedSingleEl=true;
              if (strstr(triggerList.at(i).c_str(),"HLT_Iso")) passedSingleMu=true;
              if (strstr(triggerList.at(i).c_str(),"CaloIdL")) passedAnyOther=true;
              if (strstr(triggerList.at(i).c_str(),"TrkIsoVVL")) passedAnyOther=true;
              if (strstr(triggerList.at(i).c_str(),"Triple")) passedAnyOther=true;
          }
      }
  }
  bool passedOnlySingle=((passedSingleEl && !passedAnyOther) || (passedSingleMu && !passedSingleEl && !passedAnyOther));
  bool trigConditionData = true;

  if (verbose) cout<<"checking PV"<<endl;
  const reco::Vertex *PV = 0;
  int theVertex = -1;
  for (unsigned int i=0; i<vertex->size(); i++) {
      PV = &(vertex->at(i));
      if (verbose) std::cout<<"isFake: "<<PV->isFake()<<" chi2 "<<PV->chi2()<<" ndof "<<PV->ndof()<<" rho "<<PV->position().Rho()<<" Z "<<PV->position().Z()<<endl;
      //if (PV->chi2()==0 && PV->ndof()==0) continue;
      if (PV->isFake()) continue;
      if (PV->ndof()<=4 || PV->position().Rho()>2.0 || fabs(PV->position().Z())>24.0) continue;
      theVertex=(int)i; break;
  }

  if (verbose) std::cout<<"vtx: "<<theVertex<<" trigConditionData "<<trigConditionData<<" passedTrig "<<passedTrig<<std::endl;

  if(theVertex>=0 && (isMC || (!isMC && trigConditionData))){
    if (verbose) cout<<"good PV "<<theVertex<<endl;
    PV_x =  PV->position().X();
    PV_y =  PV->position().Y();
    PV_z =  PV->position().Z();

    BS_x =  BS.position().X();
    BS_y =  BS.position().Y();
    BS_z =  BS.position().Z();
    BS_xErr =  BS.x0Error();
    BS_yErr =  BS.y0Error();
    BS_zErr =  BS.z0Error();

    BeamWidth_x = BS.BeamWidthX();
    BeamWidth_y = BS.BeamWidthY();
    BeamWidth_xErr = BS.BeamWidthXError();
    BeamWidth_yErr = BS.BeamWidthYError();

    //N Vertex
    if (verbose) {cout<<"fill nvtx histogram"<<endl;}
    nVtx = vertex->size();
    histContainer_["NVTX"]->Fill(nVtx);
    histContainer_["NVTX_RW"]->Fill(nVtx,pileupWeight);

    //MET
    if (verbose) {cout<<"get met value"<<endl;}
    if (!mets->empty()) {
        met = (*mets)[0].et();
        met_phi = (*mets)[0].phi();
        met_jesup = (*mets)[0].shiftedPt(pat::MET::JetEnUp);
        met_phi_jesup = (*mets)[0].shiftedPhi(pat::MET::JetEnUp);
        met_jesdn = (*mets)[0].shiftedPt(pat::MET::JetEnDown);
        met_phi_jesdn = (*mets)[0].shiftedPhi(pat::MET::JetEnDown);
        met_uncenup = (*mets)[0].shiftedPt(pat::MET::UnclusteredEnUp);
        met_phi_uncenup = (*mets)[0].shiftedPhi(pat::MET::UnclusteredEnUp);
        met_uncendn = (*mets)[0].shiftedPt(pat::MET::UnclusteredEnDown);
        met_phi_uncendn = (*mets)[0].shiftedPhi(pat::MET::UnclusteredEnDown);
    }

    if (verbose) cout<<"start lepton analysis"<<endl;
    vector<pat::Electron> AllElectrons; vector<pat::Muon> AllMuons;
    vector<pat::Electron> AllElectronsUnS;////uncorrected electron
    vector<pat::Photon> AllPhotons;

    AllElectrons = helper.goodLooseElectrons2012(electrons,_elecPtCut);
    AllElectronsUnS = helper.goodLooseElectrons2012(electrons,electronsUnS,_elecPtCut);
    AllMuons = helper.goodLooseMuons2012(muons,_muPtCut);
    AllPhotons = helper.goodLoosePhotons2015(photons,_phoPtCut);

    helper.cleanOverlappingLeptons(AllMuons,AllElectrons,PV);
    helper.cleanOverlappingLeptons(AllMuons,AllElectronsUnS,PV);
    recoMuons = helper.goodMuons2015_noIso_noPf(AllMuons,_muPtCut,PV,sip3dCut);
    recoElectrons = helper.goodElectrons2015_noIso_noBdt(AllElectrons,_elecPtCut,elecID,PV,iEvent,sip3dCut);
    recoElectronsUnS = helper.goodElectrons2015_noIso_noBdt(AllElectronsUnS,_elecPtCut,elecID,PV,iEvent,sip3dCut);
    recoPhotons = helper.goodPhotons2015(AllPhotons,_phoPtCut,year);

    if (verbose) cout<<AllMuons.size()<<" loose muons "<<AllElectrons.size()<<" loose electrons"<<endl;

    //sort electrons and muons by pt
    if (verbose) cout<<recoMuons.size()<<" good muons and "<<recoElectrons.size()<<" good electrons to be sorted"<<endl;
    if (verbose) cout<<"start pt-sorting leptons"<<endl;
    if (verbose) cout<<"adding muons to sorted list"<<endl;

    edm::ESHandle<TransientTrackBuilder> ttkb_recoLepton;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttkb_recoLepton);

    if( (recoMuons.size() + recoElectrons.size()) >= (uint)skimLooseLeptons ){
      if (verbose) cout<<"found >= " <<skimLooseLeptons << " leptons"<<endl;

      for(unsigned int i = 0; i < recoMuons.size(); i++){
        if (lep_ptreco.size()==0 || recoMuons[i].pt()<lep_ptreco[lep_ptreco.size()-1]){
          lep_ptreco.push_back(recoMuons[i].pt());
          lep_ptid.push_back(recoMuons[i].pdgId());
          lep_ptindex.push_back(i);
          continue;
        }
        for (unsigned int j=0; j<lep_ptreco.size(); j++){
          if (recoMuons[i].pt()>lep_ptreco[j]) {
              lep_ptreco.insert(lep_ptreco.begin()+j,recoMuons[i].pt());
              lep_ptid.insert(lep_ptid.begin()+j,recoMuons[i].pdgId());
              lep_ptindex.insert(lep_ptindex.begin()+j,i);
              break;
          }
        }
      }

      if (verbose) cout<<"adding electrons to sorted list"<<endl;
      for(unsigned int i = 0; i < recoElectrons.size(); i++) {
          if (lep_ptreco.size()==0 || recoElectrons[i].pt()<lep_ptreco[lep_ptreco.size()-1]) {
              lep_ptreco.push_back(recoElectrons[i].pt());
              lep_ptid.push_back(recoElectrons[i].pdgId());
              lep_ptindex.push_back(i);
              continue;
          }
          for (unsigned int j=0; j<lep_ptreco.size(); j++) {
              if (recoElectrons[i].pt()>lep_ptreco[j]) {
                  lep_ptreco.insert(lep_ptreco.begin()+j,recoElectrons[i].pt());
                  lep_ptid.insert(lep_ptid.begin()+j,recoElectrons[i].pdgId());
                  lep_ptindex.insert(lep_ptindex.begin()+j,i);
                  break;
              }
          }
      }

      std::vector<reco::TransientTrack> ttv_recoLepton;
      int n_lep_e = 0;
      int n_lep_mu = 0;

      for(unsigned int i = 0; i < lep_ptreco.size(); i++){
        if (verbose) cout<<"sorted lepton "<<i<<" pt "<<lep_ptreco[i]<<" id "<<lep_ptid[i]<<" index "<<lep_ptindex[i]<<endl;

        if (abs(lep_ptid[i])==11){

          if(n_lep_e < 2 && !isCode4l){
              if(helper.passTight_BDT_Id(recoElectronsUnS[lep_ptindex[i]],year)){
                  if(helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho) < isoCutEl){
                      ttv_recoLepton.push_back(ttkb_recoLepton->build(recoElectrons[lep_ptindex[i]].gsfTrack()));
                      n_lep_e++;
                  }
              }
          }

          lep_d0BS.push_back(recoElectrons[lep_ptindex[i]].gsfTrack()->dxy(beamSpot->position()));
          lep_d0PV.push_back(recoElectrons[lep_ptindex[i]].gsfTrack()->dxy(PV->position()));
          lep_numberOfValidPixelHits.push_back(recoElectrons[lep_ptindex[i]].gsfTrack()->hitPattern().numberOfValidPixelHits());
          lep_isEB.push_back(recoElectrons[lep_ptindex[i]].isEB());
          lep_isEE.push_back(recoElectrons[lep_ptindex[i]].isEE());
          lep_p.push_back(recoElectrons[lep_ptindex[i]].p());
          lep_ecalEnergy.push_back(recoElectrons[lep_ptindex[i]].correctedEcalEnergy());
          lep_id.push_back(recoElectrons[lep_ptindex[i]].pdgId());
          lep_pt.push_back(recoElectrons[lep_ptindex[i]].pt());
          lep_pterrold.push_back(recoElectrons[lep_ptindex[i]].p4Error(reco::GsfElectron::P4_COMBINATION));
          lep_pt_FromMuonBestTrack.push_back(-999);
          lep_eta_FromMuonBestTrack.push_back(-999);
          lep_phi_FromMuonBestTrack.push_back(-999);
          lep_position_x.push_back(-999);
          lep_position_y.push_back(-999);
          lep_position_z.push_back(-999);
          lep_pt_genFromReco.push_back(-999);

          double perr = 0.0;
          if (recoElectrons[lep_ptindex[i]].ecalDriven()) {
              perr = recoElectrons[lep_ptindex[i]].p4Error(reco::GsfElectron::P4_COMBINATION);
          } else{
            double ecalEnergy = recoElectrons[lep_ptindex[i]].correctedEcalEnergy();
            double err2 = 0.0;
            if (recoElectrons[lep_ptindex[i]].isEB()) {
                err2 += (5.24e-02*5.24e-02)/ecalEnergy;
                err2 += (2.01e-01*2.01e-01)/(ecalEnergy*ecalEnergy);
                err2 += 1.00e-02*1.00e-02;
            } else if (recoElectrons[lep_ptindex[i]].isEE()) {
                err2 += (1.46e-01*1.46e-01)/ecalEnergy;
                err2 += (9.21e-01*9.21e-01)/(ecalEnergy*ecalEnergy);
                err2 += 1.94e-03*1.94e-03;
            }
            perr = ecalEnergy * sqrt(err2);
          }

          double pterr = perr*recoElectrons[lep_ptindex[i]].pt()/recoElectrons[lep_ptindex[i]].p();
          lep_pterr.push_back(pterr);
          lep_eta.push_back(recoElectrons[lep_ptindex[i]].eta());
          lep_phi.push_back(recoElectrons[lep_ptindex[i]].phi());
          lep_mass.push_back(recoElectrons[lep_ptindex[i]].mass());
          lepFSR_pt.push_back(recoElectrons[lep_ptindex[i]].pt());
          lepFSR_eta.push_back(recoElectrons[lep_ptindex[i]].eta());
          lepFSR_phi.push_back(recoElectrons[lep_ptindex[i]].phi());
          lepFSR_mass.push_back(recoElectrons[lep_ptindex[i]].mass());
          lepFSR_ID.push_back(11);

          if (isoConeSizeEl==0.4) {
              lep_RelIso.push_back(helper.pfIso(recoElectrons[lep_ptindex[i]],elRho));
              lep_RelIsoNoFSR.push_back(helper.pfIso(recoElectrons[lep_ptindex[i]],elRho));
              lep_isoCH.push_back(recoElectrons[lep_ptindex[i]].chargedHadronIso());
              lep_isoNH.push_back(recoElectrons[lep_ptindex[i]].neutralHadronIso());
              lep_isoPhot.push_back(recoElectrons[lep_ptindex[i]].photonIso());
              lep_isoPU.push_back(recoElectrons[lep_ptindex[i]].puChargedHadronIso());
              lep_isoPUcorr.push_back(helper.getPUIso(recoElectrons[lep_ptindex[i]],elRho));
          } else if (isoConeSizeEl==0.3) {
              lep_RelIso.push_back(helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho));
              lep_RelIsoNoFSR.push_back(helper.pfIso03(recoElectrons[lep_ptindex[i]],elRho));
              lep_isoCH.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumChargedHadronPt);
              lep_isoNH.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumNeutralHadronEt);
              lep_isoPhot.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumPhotonEt);
              lep_isoPU.push_back(recoElectrons[lep_ptindex[i]].pfIsolationVariables().sumPUPt);
              lep_isoPUcorr.push_back(helper.getPUIso03(recoElectrons[lep_ptindex[i]],elRho));
          }

          lep_MiniIso.push_back(helper.miniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&recoElectrons[lep_ptindex[i]]), 0.05, 0.2, 10., rhoSUS, false));
          lep_Sip.push_back(helper.getSIP3D(recoElectrons[lep_ptindex[i]]));
          //lep_mva.push_back(recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Autumn18IdIsoValues"));
          //cout<<EleBDT_name_161718<<endl;
          lep_mva.push_back(recoElectrons[lep_ptindex[i]].userFloat(EleBDT_name_161718.c_str()));
          lep_ecalDriven.push_back(recoElectrons[lep_ptindex[i]].ecalDriven());
          //lep_tightId.push_back(helper.passTight_BDT_Id(recoElectrons[lep_ptindex[i]],recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Autumn18IdIsoValues"), year));
          lep_tightId.push_back(helper.passTight_BDT_Id(recoElectronsUnS[lep_ptindex[i]],year));
          //lep_tightId_old.push_back(helper.passTight_BDT_Id(recoElectronsUnS[lep_ptindex[i]],year));
          //cout<<"old "<<recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values") <<" new" <<recoElectrons[lep_ptindex[i]].userFloat("ElectronMVAEstimatorRun2Spring16HZZV1Values")<<endl;
          lep_tightIdSUS.push_back(helper.passTight_Id_SUS(recoElectrons[lep_ptindex[i]],elecID,PV,BS,*theConversions, year));
          //lep_tightIdHiPt.push_back(recoElectrons[lep_ptindex[i]].electronID("heepElectronID-HEEPV70"));
          lep_tightIdHiPt.push_back(recoElectrons[lep_ptindex[i]].electronID(heepID_name_161718.c_str()));
          lep_ptRatio.push_back(helper.ptRatio(recoElectrons[lep_ptindex[i]],jets,true)); // no L2L3 yet
          lep_ptRel.push_back(helper.ptRel(recoElectrons[lep_ptindex[i]],jets,true)); // no L2L3 yet
          lep_dataMC.push_back(helper.dataMC(recoElectrons[lep_ptindex[i]],hElecScaleFac,hElecScaleFac_Cracks,hElecScaleFacGsf,hElecScaleFacGsfLowET));
          lep_dataMCErr.push_back(helper.dataMCErr(recoElectrons[lep_ptindex[i]],hElecScaleFac,hElecScaleFac_Cracks));
          lep_genindex.push_back(-1.0);

        }

        if (abs(lep_ptid[i])==13){
          SingleTrackVertexConstraint stvc;
          SingleTrackVertexConstraint::BTFtuple a = stvc.constrain(ttkb_recoLepton->build(recoMuons[lep_ptindex[i]].muonBestTrack()), BS);
          reco::Track single_trk_BS = a.get<1>().track();
          TLorentzVector tmp;
          tmp.SetPxPyPzE(single_trk_BS.px(), single_trk_BS.py(), single_trk_BS.pz(), KalmanEnergy(single_trk_BS.px(), single_trk_BS.py(), single_trk_BS.pz(), recoMuons[lep_ptindex[i]].mass()));
          singleBS_RecoLep_pt.push_back(tmp.Pt());
          singleBS_RecoLep_ptError.push_back(single_trk_BS.ptError());
          singleBS_RecoLep_eta.push_back(tmp.Eta());
          singleBS_RecoLep_phi.push_back(tmp.Phi());
          singleBS_RecoLep_d0.push_back(single_trk_BS.dxy(BS.position()));
          singleBS_RecoLep_mass.push_back(recoMuons[lep_ptindex[i]].mass());

          if(n_lep_mu < 2 && !isCode4l){
              if(helper.passTight_Id(recoMuons[lep_ptindex[i]],PV)){
                  if(helper.pfIso03(recoMuons[lep_ptindex[i]],muRho) < isoCutMu){
                      ttv_recoLepton.push_back(ttkb_recoLepton->build(recoMuons[lep_ptindex[i]].muonBestTrack()));
                      n_lep_mu++;
                  }
              }
          }

          lep_d0BS.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->dxy(beamSpot->position()));
          lep_d0PV.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->dxy(PV->position()));
          lep_numberOfValidPixelHits.push_back(recoMuons[lep_ptindex[i]].innerTrack()->hitPattern().numberOfValidPixelHits());
          lep_trackerLayersWithMeasurement.push_back(recoMuons[lep_ptindex[i]].innerTrack()->hitPattern().trackerLayersWithMeasurement());
          lep_isEB.push_back(0);
          lep_isEE.push_back(0);
          lep_p.push_back(recoMuons[lep_ptindex[i]].p());
          lep_ecalEnergy.push_back(0);
          lep_id.push_back(recoMuons[lep_ptindex[i]].pdgId());
          lep_pt.push_back(recoMuons[lep_ptindex[i]].pt());
          lep_pterrold.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->ptError());

          lep_pt_FromMuonBestTrack.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->pt());
          lep_eta_FromMuonBestTrack.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->eta());
          lep_phi_FromMuonBestTrack.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->phi());
          lep_position_x.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->vx());
          lep_position_y.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->vy());
          lep_position_z.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->vz());

          auto gen_lep = recoMuons[lep_ptindex[i]].genParticle();
          if(gen_lep != 0)
              lep_pt_genFromReco.push_back(gen_lep->pt());
          else
              lep_pt_genFromReco.push_back(-999);

          if (recoMuons[lep_ptindex[i]].hasUserFloat("correctedPtError")){
            lep_pterr.push_back(recoMuons[lep_ptindex[i]].userFloat("correctedPtError"));
          } else{
            lep_pterr.push_back(recoMuons[lep_ptindex[i]].muonBestTrack()->ptError());
          }

          lep_eta.push_back(recoMuons[lep_ptindex[i]].eta());
          lep_phi.push_back(recoMuons[lep_ptindex[i]].phi());
          if (recoMuons[lep_ptindex[i]].mass()<0.105) cout<<"muon mass: "<<recoMuons[lep_ptindex[i]].mass()<<endl;
          lep_mass.push_back(recoMuons[lep_ptindex[i]].mass());
          lepFSR_pt.push_back(recoMuons[lep_ptindex[i]].pt());
          lepFSR_eta.push_back(recoMuons[lep_ptindex[i]].eta());
          lepFSR_phi.push_back(recoMuons[lep_ptindex[i]].phi());
          lepFSR_mass.push_back(recoMuons[lep_ptindex[i]].mass());
          lepFSR_ID.push_back(13);

          if (isoConeSizeMu==0.4) {
              lep_RelIso.push_back(helper.pfIso(recoMuons[lep_ptindex[i]],muRho));
              lep_RelIsoNoFSR.push_back(helper.pfIso(recoMuons[lep_ptindex[i]],muRho));
              lep_isoCH.push_back(recoMuons[lep_ptindex[i]].chargedHadronIso());
              lep_isoNH.push_back(recoMuons[lep_ptindex[i]].neutralHadronIso());
              lep_isoPhot.push_back(recoMuons[lep_ptindex[i]].photonIso());
              lep_isoPU.push_back(recoMuons[lep_ptindex[i]].puChargedHadronIso());
              lep_isoPUcorr.push_back(helper.getPUIso(recoMuons[lep_ptindex[i]],muRho));
          } else if (isoConeSizeMu==0.3) {
              lep_RelIso.push_back(helper.pfIso03(recoMuons[lep_ptindex[i]],muRho));
              lep_RelIsoNoFSR.push_back(helper.pfIso03(recoMuons[lep_ptindex[i]],muRho));
              lep_isoCH.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumChargedHadronPt);
              lep_isoNH.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumNeutralHadronEt);
              lep_isoPhot.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumPhotonEt);
              lep_isoPU.push_back(recoMuons[lep_ptindex[i]].pfIsolationR03().sumPUPt);
              lep_isoPUcorr.push_back(helper.getPUIso03(recoMuons[lep_ptindex[i]],muRho));
          }

          lep_MiniIso.push_back(helper.miniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&recoMuons[lep_ptindex[i]]), 0.05, 0.2, 10., rhoSUS, false));
          lep_Sip.push_back(helper.getSIP3D(recoMuons[lep_ptindex[i]]));
          //lep_mva.push_back(recoMuons[lep_ptindex[i]].isPFMuon());
          //lep_mva.push_back(helper.get_Muon_MVA_Value(recoMuons[lep_ptindex[i]],vertex,muRho,year,PV));
          lep_mva.push_back(recoMuons[lep_ptindex[i]].isPFMuon());
          lep_ecalDriven.push_back(0);
          lep_tightId.push_back(helper.passTight_Id(recoMuons[lep_ptindex[i]],PV));
          //lep_tightId.push_back(helper.passTight_BDT_Id(recoMuons[lep_ptindex[i]],vertex,muRho,year,PV));
          //lep_tightId_old.push_back(helper.passTight_Id(recoMuons[lep_ptindex[i]],PV));
          lep_tightIdSUS.push_back(helper.passTight_Id_SUS(recoMuons[lep_ptindex[i]],PV));
          lep_tightIdHiPt.push_back(recoMuons[lep_ptindex[i]].isHighPtMuon(*PV));
          lep_ptRatio.push_back(helper.ptRatio(recoMuons[lep_ptindex[i]],jets,true)); // no L2L3 yet
          lep_ptRel.push_back(helper.ptRel(recoMuons[lep_ptindex[i]],jets,true)); // no L2L3 yet
          lep_dataMC.push_back(helper.dataMC(recoMuons[lep_ptindex[i]],hMuScaleFac));
          lep_dataMCErr.push_back(helper.dataMCErr(recoMuons[lep_ptindex[i]],hMuScaleFacUnc));
          lep_genindex.push_back(-1.0);

        }

        if (verbose){
          cout<<" eta: "<<lep_eta[i]<<" phi: "<<lep_phi[i];
          if(abs(lep_ptid[i])==11)
              cout<<" eSuperClusterOverP: "<<recoElectrons[lep_ptindex[i]].eSuperClusterOverP()<<" ecalEnergy: "<<recoElectrons[lep_ptindex[i]].ecalEnergy()<<" p: "<<recoElectrons[lep_ptindex[i]].p();
          cout<<" RelIso: "<<lep_RelIso[i]<<" isoCH: "<<lep_isoCH[i]<<" isoNH: "<<lep_isoNH[i]
              <<" isoPhot: "<<lep_isoPhot[i]<<" lep_isoPU: "<<lep_isoPU[i]<<" isoPUcorr: "<<lep_isoPUcorr[i]<<" Sip: "<<lep_Sip[i]
              <<" MiniIso: "<<lep_MiniIso[i]<<" ptRatio: "<<lep_ptRatio[i]<<" ptRel: "<<lep_ptRel[i]<<" lep_mva: "<<lep_mva[i];
          if(abs(lep_ptid[i])==11)
              cout<<" SCeta: "<<recoElectrons[lep_ptindex[i]].superCluster()->eta()<<" dxy: "<<recoElectrons[lep_ptindex[i]].gsfTrack()->dxy(PV->position())<<" dz: "<<recoElectrons[lep_ptindex[i]].gsfTrack()->dz(PV->position());
          if(abs(lep_ptid[i])==11)
              cout<<" Rho: "<<elRho<<" EleBDT_name: "<<EleBDT_name_161718<<" Uncorrected electron pt: "<<recoElectronsUnS[lep_ptindex[i]].pt();
          if(abs(lep_ptid[i])==13)
              cout<<" Rho: "<<muRho;
          cout<<" dataMC: "<<lep_dataMC[i]<<" dataMCErr: "<<lep_dataMCErr[i];
          cout<<" lep_pterr: "<<lep_pterr[i]<<" lep_pterrold: "<<lep_pterrold[i]<<" lep_tightIdHiPt: "<<lep_tightIdHiPt[i]<<endl;
          if((abs(lep_ptid[i])==13)&&lep_pt[i]>200)
              cout<<"Muon pt over 200 isTrackerHighPtID? "<<helper.isTrackerHighPt(recoMuons[lep_ptindex[i]],PV)<<endl;
        }

      }

      if (verbose) cout<<"adding photons to sorted list"<<endl;
      for(int i = 0; i < (int)recoPhotons.size(); i++) {
          pho_pt.push_back(recoPhotons[i].pt());
          pho_eta.push_back(recoPhotons[i].eta());
          pho_phi.push_back(recoPhotons[i].phi());
          photonCutBasedIDLoose.push_back(recoPhotons[i].photonID("cutBasedPhotonID-Fall17-94X-V2-loose"));
      }

      if (doTriggerMatching){

        if (verbose) cout<<"start trigger matching"<<endl;
        // trigger Matching
        if (verbose) cout<<"start trigger matching"<<endl;
        for(unsigned int i = 0; i < lep_pt.size(); i++) {

            TLorentzVector reco;
            reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

            double reco_eta = reco.Eta();
            double reco_phi = reco.Phi();

            std::string filtersMatched = "";

            for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
                double hlt_eta = obj.eta();
                double hlt_phi = obj.phi();
                double dR =  deltaR(reco_eta,reco_phi,hlt_eta,hlt_phi);
                if (dR<0.5) {
                    obj.unpackFilterLabels(iEvent, *trigger);
                    for (unsigned h = 0; h < obj.filterLabels().size(); ++h) filtersMatched += obj.filterLabels()[h];
                }
            }

            if (verbose) cout<<"Trigger matching lep id: "<<lep_id[i]<<" pt: "<<reco.Pt()<<" filters: "<<filtersMatched<<endl;
            lep_filtersMatched.push_back(filtersMatched);

        }

      }

      // GEN matching
      if(isMC) {
          if (verbose) cout<<"begin gen matching"<<endl;
          // for each reco lepton find the nearest gen lepton with same ID
          for(unsigned int i = 0; i < lep_pt.size(); i++) {

              double minDr=9999.0;

              TLorentzVector reco, gen;
              reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

              for (unsigned int j = 0; j < GENlep_id.size(); j++) {

                  if (GENlep_id[j]!=lep_id[i]) continue;

                  gen.SetPtEtaPhiM(GENlep_pt[j],GENlep_eta[j],GENlep_phi[j],GENlep_mass[j]);
                  double thisDr = deltaR(reco.Eta(),reco.Phi(),gen.Eta(),gen.Phi());

                  if (thisDr<minDr && thisDr<0.5) {
                      lep_genindex[i]=j;
                      minDr=thisDr;
                  }

              } // all gen leptons

          } // all reco leptons

          // for each reco lepton find the nearest gen particle and save its ID and mothers
          for(unsigned int i = 0; i < lep_pt.size(); i++) {

              double minDr=9999.0;

              TLorentzVector reco;
              reco.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

              reco::GenParticleCollection::const_iterator genPart;
              int j = -1;
              int tmpPdgId = 0;
              int tmpMomId = 0;
              int tmpMomMomId = 0;
              for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
                  j++;
                  double thisDr = deltaR(reco.Eta(),reco.Phi(),genPart->eta(),genPart->phi());

                  if (thisDr<minDr && thisDr<0.3) {
                      tmpPdgId=genPart->pdgId();
                      tmpMomId=genAna.MotherID(&prunedgenParticles->at(j));
                      tmpMomMomId=genAna.MotherMotherID(&prunedgenParticles->at(j));

                      minDr=thisDr;
                  }

              } // all gen particles
              // storing the matches
              lep_matchedR03_PdgId.push_back(tmpPdgId);
              lep_matchedR03_MomId.push_back(tmpMomId);
              lep_matchedR03_MomMomId.push_back(tmpMomMomId);

          } // all reco leptons

          if (verbose) {cout<<"finished gen matching"<<endl;}
      } //isMC

      unsigned int Nleptons = lep_pt.size();
      // FSR Photons
      if(doFsrRecovery) {

          if (verbose) cout<<"checking "<<photonsForFsr->size()<<" fsr photon candidates"<<endl;

          // try to find an fsr photon for each lepton
          for (unsigned int i=0; i<Nleptons; i++) {

              TLorentzVector thisLep;
              thisLep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

              double minDrOEt2 = 999.0; double selectedPhotonIso=9999.0; double selectedPhotonDr=9999.0;
              bool selected=false; pat::PFParticle selectedPhoton;
              int index=-1;
              for(edm::View<pat::PFParticle>::const_iterator phot=photonsForFsr->begin(); phot!=photonsForFsr->end(); ++phot) {
                  index++;

                  // preselection
                  if (fabs(phot->eta()) > 2.4) continue;
                  if (phot->pt()<2.0) continue;
                  double fsrDr = deltaR(thisLep.Eta(), thisLep.Phi(), phot->eta(), phot->phi());
                  if ( verbose ) cout<<"fsr photon cand index: "<<index<<", pt: "<<phot->pt()<<" eta: "<<phot->eta()<<" phi: "<<phot->phi()<<"; fsrDr: "<<fsrDr<<endl;;
                  if (fsrDr>0.5) continue;

                  // check super cluster veto against all electrons for each photon
                  // at the same time check that this is the closest lepton for this photon
                  bool matched=false;
                  bool closest=true;

                  for (unsigned int j=0; j<Nleptons; j++) {

                      TLorentzVector otherLep;
                      otherLep.SetPtEtaPhiM(lep_pt[j],lep_eta[j],lep_phi[j],lep_mass[j]);

                      double fsrDrOther = deltaR(otherLep.Eta(), otherLep.Phi(), phot->eta(), phot->phi());
                      if (j!=i && fsrDrOther<fsrDr) {closest=false;}

                      if ( abs(lep_id[(int)j])==11) {

                          for(size_t ecand = 0; ecand < recoElectrons[lep_ptindex[j]].associatedPackedPFCandidates().size(); ecand++){


                              double ecandpt = recoElectrons[lep_ptindex[j]].associatedPackedPFCandidates()[ecand]->pt();
                              double ecandeta = recoElectrons[lep_ptindex[j]].associatedPackedPFCandidates()[ecand]->eta();
                              double ecandphi = recoElectrons[lep_ptindex[j]].associatedPackedPFCandidates()[ecand]->phi();
                              if (abs(ecandpt-phot->pt())<1e-10 && abs(ecandeta-phot->eta())<1e-10 && abs(ecandphi-phot->phi())<1e-10) matched=true;

                          }
                      }
                  }
                  if (matched) continue;
                  if (!closest) continue;

                  // comput iso, dR/ET^2
                  double photoniso = helper.photonPfIso03(*phot,pfCands)/phot->pt();
                  double fsrDrOEt2 = fsrDr/(phot->pt()*phot->pt());

                  // fill all fsr photons before iso, dR/Et^2 cuts
                  allfsrPhotons_dR.push_back(fsrDr);
                  allfsrPhotons_iso.push_back(photoniso);
                  allfsrPhotons_pt.push_back(phot->pt());

                  // require photon iso, dR/Et^2
                  if (photoniso>1.8) continue;
                  if (fsrDrOEt2>0.012) continue;

                  // this photon is now a good one, check if it is the best one

                  if ( verbose) cout<<"fsr photon cand, pt: "<<phot->pt()<<" eta: "<<phot->eta()<<" phi: "<<phot->phi()
                      //<<" isoCHPUNoPU: "<<phot->userFloat("fsrPhotonPFIsoChHadPUNoPU03pt02")
                      //<<" isoNHPhoton: "<<phot->userFloat("fsrPhotonPFIsoNHadPhoton03")
                      <<" photoniso: "<<photoniso<<" DrOEt2: "<< fsrDrOEt2 <<endl;

                  if( fsrDrOEt2 < minDrOEt2 ) {
                      selected = true;
                      selectedPhoton=(*phot);
                      selectedPhotonIso=photoniso;
                      selectedPhotonDr=fsrDr;
                      minDrOEt2 = fsrDrOEt2;
                      if (verbose) cout<<"****selected fsr: "<<i<<endl;
                      if (verbose) cout<<"photoniso: "<<photoniso<<" fsr dR over Et^2: "<<fsrDrOEt2<<" fsr dR: "<<fsrDr<<endl;
                  }

              } // all photons

              if (selected) {
                  nFSRPhotons++;
                  fsrPhotons.push_back(selectedPhoton);
                  fsrPhotons_pt.push_back(selectedPhoton.pt());
                  double perr = PFEnergyResolution().getEnergyResolutionEm(selectedPhoton.energy(), selectedPhoton.eta());
                  double pterr = perr*selectedPhoton.pt()/selectedPhoton.p();
                  fsrPhotons_pterr.push_back(pterr);
                  fsrPhotons_eta.push_back(selectedPhoton.eta());
                  fsrPhotons_phi.push_back(selectedPhoton.phi());
                  fsrPhotons_lepindex.push_back((int)i);
                  fsrPhotons_dR.push_back(selectedPhotonDr);
                  fsrPhotons_iso.push_back(selectedPhotonIso);
                  TLorentzVector phofsr;
                  phofsr.SetPtEtaPhiM(selectedPhoton.pt(),selectedPhoton.eta(),selectedPhoton.phi(),0.0);
                  TLorentzVector lepfsr;
                  lepfsr = thisLep+phofsr;
                  lepFSR_pt[i] = lepfsr.Pt();
                  lepFSR_eta[i] = lepfsr.Eta();
                  lepFSR_phi[i] = lepfsr.Phi();
                  lepFSR_mass[i] = lepfsr.M();

                  //                         std::cout<<"PHOTON = "<<selectedPhoton.eta()<<"\t"<<selectedPhoton.phi()<<"\t"<<selectedPhoton.pt()<<std::endl;
                  //                         std::cout<<"LEPTON = "<<lep_eta[i]<<"\t"<<lep_phi[i]<<"\t"<<lep_pt[i]<<"\t"<<fabs(lep_id[i])<<std::endl;
                  // 						std::cout<<"DeltaR = "<<minDrOEt2 * selectedPhoton.pt() * selectedPhoton.pt()<<std::endl;
                  fsrmap[i] = phofsr;
                  if (verbose) cout<<"****selected fsr: "<<i<<endl;
                  if (verbose) cout<<"phofsr pt: "<<phofsr.Pt()<<" eta: "<<phofsr.Eta()<<" phi: "<<phofsr.Phi()<<" mass: "<<phofsr.M()<<endl;
                  if (verbose) cout<<"lep pt: "<<thisLep.Pt()<<" eta: "<<thisLep.Eta()<<" phi: "<<thisLep.Phi()<<" mass: "<<thisLep.M()<<endl;
                  if (verbose) cout<<"lep+fsr pt: "<<lepFSR_pt[i]<<" eta: "<<lepFSR_eta[i]<<" phi: "<<lepFSR_phi[i]<<" mass: "<<lepfsr.M()<<endl;


              }

          } // all leptons

          //std::cout<<"Number = "<<fsrPhotons_eta.size()<<std::endl;

          // subtract selected photons from all leptons isolations
          for (unsigned int i=0; i<Nleptons; i++) {

              TLorentzVector lep_nofsr;
              lep_nofsr.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);

              double isoFSR=0.0;
              for (unsigned int j=0; j<fsrPhotons.size(); j++) {
                  double fsrDr = deltaR(lep_nofsr.Eta(),lep_nofsr.Phi(),fsrPhotons[j].eta(),fsrPhotons[j].phi());
                  bool isoVeto=true;
                  if (abs(lep_id[i])==13 && fsrDr>0.01) isoVeto=false;
                  if (abs(lep_id[i])==11 && (abs(recoElectrons[lep_ptindex[i]].superCluster()->eta())<1.479 || fsrDr>0.08)) isoVeto=false;
                  if (fsrDr<((abs(lep_id[i])==11)?isoConeSizeEl:isoConeSizeMu) && !isoVeto) isoFSR += fsrPhotons[j].pt();
              }

              double RelIsoNoFSR = (lep_isoCH[i]+std::max(lep_isoNH[i]+lep_isoPhot[i]-lep_isoPUcorr[i]-isoFSR,0.0))/lep_nofsr.Pt();
              lep_RelIsoNoFSR[i] = RelIsoNoFSR;
              if (verbose) cout<<"lep pt: "<<lep_nofsr.Pt()<<" eta: "<<lep_nofsr.Eta()<<" phi: "<<lep_nofsr.Phi()<<" RelIsoNoFSR: "<<RelIsoNoFSR<<" lep mva: "<<lep_mva[i]<<" tightId? "<<lep_tightId[i]<<endl;
          }
          if (verbose) {cout<<"finished filling fsr photon candidates"<<endl;}
      } // doFsrRecovery

      // count tight ID iso leptons
      uint ntight=0;
      for (unsigned int i=0; i<Nleptons; i++) {
          if (abs(lep_id[i])==11 && lep_RelIsoNoFSR[i]<isoCutEl && lep_tightId[i]==1) ntight+=1;
          if (abs(lep_id[i])==13 && lep_RelIsoNoFSR[i]<isoCutMu && lep_tightId[i]==1) ntight+=1;
      }

      if ( ntight >= (uint)skimTightLeptons ){

        // creat vectors for selected objects
        vector<pat::Muon> selectedMuons;
        vector<pat::Electron> selectedElectrons;


        // Jets
        if (verbose) cout<<"begin filling jet candidates"<<endl;

        vector<pat::Jet> goodJets;
        vector<float> patJetQGTagger, patJetaxis2, patJetptD;
        vector<float> goodJetQGTagger, goodJetaxis2, goodJetptD;
        vector<int> patJetmult, goodJetmult;

        for(auto jet = jets->begin();  jet != jets->end(); ++jet){
            edm::RefToBase<pat::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(jets, jet - jets->begin()));
            float qgLikelihood = (*qgHandle)[jetRef];
            float axis2 = (*axis2Handle)[jetRef];
            float ptD = (*ptDHandle)[jetRef];
            int mult = (*multHandle)[jetRef];
            patJetQGTagger.push_back(qgLikelihood);
            patJetaxis2.push_back(axis2);
            patJetmult.push_back(mult);
            patJetptD.push_back(ptD);
        }

        for(unsigned int i = 0; i < jets->size(); ++i){
          const pat::Jet & jet = jets->at(i);

          //JetID ID
          if (verbose) cout<<"checking jetid..."<<endl;
          float jpumva=0.;
          bool passPU;
          if (doJEC && (year==2017 || year==2018 )) {
              passPU = bool(jet.userInt("pileupJetIdUpdated:fullId") & (1 << 0));
              jpumva=jet.userFloat("pileupJetIdUpdated:fullDiscriminant");
          } else if (doJEC && year==2016) {
              passPU = bool(jet.userInt("pileupJetIdUpdated:fullId") & (1 << 2));
              jpumva=jet.userFloat("pileupJetIdUpdated:fullDiscriminant");
          } else {
              passPU = bool(jet.userInt("pileupJetId:fullId") & (1 << 0));
              jpumva=jet.userFloat("pileupJetId:fullDiscriminant");
          }
          if (verbose) cout<< " jet pu mva  "<<jpumva <<endl;

          if (verbose) cout<<"pt: "<<jet.pt()<<" eta: "<<jet.eta()<<" phi: "<<jet.phi()<<" passPU: "<<passPU<<" jetid: "<<jetHelper.patjetID(jet,year)<<endl;

          if( jetHelper.patjetID(jet,year)>=jetIDLevel ){
            if (verbose) cout<<"passed pf jet id and pu jet id"<<endl;
            // apply loose pt cut here (10 GeV cut is already applied in MINIAOD) since we are before JES/JER corrections
            if(fabs(jet.eta()) < jeteta_cut){ //move all pt cut after JES/JER corrections
              // apply scale factor for PU Jets by demoting 1-data/MC % of jets jets in certain pt/eta range
              // Configured now that SF is 1.0
              if (verbose) cout<<"adding pu jet scale factors..."<<endl;
              bool dropit=false;
              if (abs(jet.eta())>3.0 && isMC){
                TRandom3 rand;
                rand.SetSeed(abs(static_cast<int>(sin(jet.phi())*100000)));
                float coin = rand.Uniform(1.);
                if (jet.pt()>=20.0 && jet.pt()<36.0 && coin>1.0) dropit=true;
                if (jet.pt()>=36.0 && jet.pt()<50.0 && coin>1.0) dropit=true;
                if (jet.pt()>=50.0 && coin>1.0) dropit=true;
              }

              if (!dropit){
                if (verbose) cout<<"adding jet candidate, pt: "<<jet.pt()<<" eta: "<<jet.eta()<<endl;
                goodJets.push_back(jet);
                goodJetQGTagger.push_back(patJetQGTagger[i]);
                goodJetaxis2.push_back(patJetaxis2[i]);
                goodJetptD.push_back(patJetptD[i]);
                goodJetmult.push_back(patJetmult[i]);
              } // pu jet scale factor
            } // pass loose pt cut
          } // pass loose pf jet id and pu jet id
        } // all jets

        vector<pat::Jet> selectedMergedJets;

        for(unsigned int i = 0; i < mergedjets->size(); ++i) { // all merged jets

        const pat::Jet & mergedjet = mergedjets->at(i);
        double pt = double(mergedjet.pt());
        double eta = double(mergedjet.eta());

        if(pt>200 && abs(eta)<2.5) selectedMergedJets.push_back(mergedjet);

        } // all merged jets

        //Set All the Variables for Saved Trees (after finding higgs candidate)
        if (verbose) cout<<"begin setting tree variables"<<endl;
        setTreeVariables(iEvent, iSetup, selectedMuons, selectedElectrons, recoMuons, recoElectrons, goodJets, goodJetQGTagger,goodJetaxis2, goodJetptD, goodJetmult, selectedMergedJets, selectedFsrMap, pfCands, jets);
        if (verbose) cout<<"finshed setting tree variables"<<endl;

        eventWeight = genWeight*crossSection*pileupWeight*prefiringWeight;

        if (!isMC) passedEventsTree_All->Fill();

      } // 2 tight ID
      else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed  ntight ID"<<endl;}
    } //if 2 lepID
    else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed  nloose ID"<<endl;}
  } //primary vertex,notDuplicate
  else { if (verbose) cout<<Run<<":"<<LumiSect<<":"<<Event<<" failed primary vertex"<<endl;}

  if (isMC) passedEventsTree_All->Fill();
  if (nEventsTotal==1000.0) passedEventsTree_All->OptimizeBaskets();

}

// ------------ method called once each job just before starting event loop  ------------
void UFHZZ4LAna::beginJob(){
    using namespace edm;
    using namespace std;
    using namespace pat;

    bookPassedEventTree("passedEvents", passedEventsTree_All);

    firstEntry = true;
}
// ------------ method called once each job just after ending the event loop  ------------
void UFHZZ4LAna::endJob(){
    histContainer_["NEVENTS"]->SetBinContent(1,nEventsTotal);
    histContainer_["NEVENTS"]->GetXaxis()->SetBinLabel(1,"N Events in Sample");
    histContainer_["SUMWEIGHTS"]->SetBinContent(1,sumWeightsTotal);
    histContainer_["SUMWEIGHTSPU"]->SetBinContent(1,sumWeightsTotalPU);
    histContainer_["SUMWEIGHTS"]->GetXaxis()->SetBinLabel(1,"sum Weights in Sample");
    histContainer_["SUMWEIGHTSPU"]->GetXaxis()->SetBinLabel(1,"sum Weights PU in Sample");
}

void UFHZZ4LAna::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){

    //massErr.init(iSetup);
    if (isMC) {
        edm::Handle<LHERunInfoProduct> run;
        typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
        try {

            int pos=0;
            iRun.getByLabel( edm::InputTag("externalLHEProducer"), run );
            LHERunInfoProduct myLHERunInfoProduct = *(run.product());
            typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
            for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
                if (verbose) std::cout << iter->tag() << std::endl;
                std::vector<std::string> lines = iter->lines();
                for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
                    std::string pdfid=lines.at(iLine);
                    if (pdfid.substr(1,6)=="weight" && pdfid.substr(8,2)=="id") {
                        if (verbose) std::cout<<pdfid<<std::endl;
                        std::string pdf_weight_id = pdfid.substr(12,4);
                        int pdf_weightid=atoi(pdf_weight_id.c_str());
                        if (verbose) std::cout<<"parsed id: "<<pdf_weightid<<std::endl;
                        if (pdf_weightid==2001) {posNNPDF=int(pos);}
                        pos+=1;
                    }
                }
            }
        }
        catch(...) {
            std::cout<<"No LHERunInfoProduct"<<std::endl;
        }
    }

}

// ------------ method called when ending the processing of a run  ------------
void UFHZZ4LAna::endRun(const edm::Run& iRun, edm::EventSetup const&){
}


// ------------ method called when starting to processes a luminosity block  ------------
void UFHZZ4LAna::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}

// ------------ method called when ending the processing of a luminosity block  ------------
void UFHZZ4LAna::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg,edm::EventSetup const& eSetup){
    using namespace edm;
    using namespace std;
    // Keep track of all the events run over
    edm::Handle<MergeableCounter> numEventsCounter;
    lumiSeg.getByLabel("nEventsTotal", numEventsCounter);
    if(numEventsCounter.isValid()) {
        std::cout<<"numEventsCounter->value "<<numEventsCounter->value<<endl;
        nEventsTotal += numEventsCounter->value;
    }
}

// ============================ UF Functions ============================= //
void UFHZZ4LAna::setTreeVariables( const edm::Event& iEvent, const edm::EventSetup& iSetup,
                                  std::vector<pat::Muon> selectedMuons, std::vector<pat::Electron> selectedElectrons,
                                  std::vector<pat::Muon> recoMuons, std::vector<pat::Electron> recoElectrons,
                                  std::vector<pat::Jet> goodJets, std::vector<float> goodJetQGTagger,
                                  std::vector<float> goodJetaxis2, std::vector<float> goodJetptD, std::vector<int> goodJetmult,
                                  std::vector<pat::Jet> selectedMergedJets,
                                  std::map<unsigned int, TLorentzVector> selectedFsrMap,
                                  edm::Handle<pat::PackedCandidateCollection> &pfcands,
                                  edm::Handle<edm::View<pat::Jet> > &jets)
{
  using namespace edm;
  using namespace pat;
  using namespace std;

  // Jet Info
  double tempDeltaR = 999.0;
  // jets reco new declaration, TJ
  //int jet1index_jesup=-1, jet2index_jesup=-1;
  //float jet1pt_jesup=-1.0, jet2pt_jesup=-1.0;
  //int jet1index2p5_jesup=-1, jet2index2p5_jesup=-1;
  //float jet1pt2p5_jesup=-1.0, jet2pt2p5_jesup=-1.0;

  //int jet1index_jesdn=-1, jet2index_jesdn=-1;
  //float jet1pt_jesdn=-1.0, jet2pt_jesdn=-1.0;
  //int jet1index2p5_jesdn=-1, jet2index2p5_jesdn=-1;
  //float jet1pt2p5_jesdn=-1.0, jet2pt2p5_jesdn=-1.0;

  // jets reco new declaration end, TJ

  vector<pat::Jet> goodJets_JECJER_pt30_eta4p7;
  for( unsigned int k = 0; k < goodJets.size(); k++) {
      if (verbose) cout<<"jet pt: "<<goodJets[k].pt()<<" eta: "<<goodJets[k].eta()<<" phi: "<<goodJets[k].phi()<<endl;
      bool isclean_H4l = true;

      TLorentzVector thisJet;
      thisJet.SetPtEtaPhiM(goodJets[k].pt(),goodJets[k].eta(),goodJets[k].phi(),goodJets[k].mass());

      // check overlap with tight ID isolated leptons OR higgs candidate leptons
      unsigned int Nleptons = lep_pt.size();
      for (unsigned int i=0; i<Nleptons; i++) {
          bool passed_idiso=true;
          if (abs(lep_id[i])==13 && lep_RelIso[i]>isoCutMu) passed_idiso=false;
          if (abs(lep_id[i])==11 && lep_RelIso[i]>isoCutEl) passed_idiso=false;
          if (!(lep_tightId[i])) passed_idiso=false;
          //for (int l=0; l<4; l++) {
          //    if ((int)i==lep_Hindex[l]) cand_lep=true;
          //}
          //if (!(passed_idiso || cand_lep)) continue;
          if (!passed_idiso) continue;
          TLorentzVector thisLep;
          thisLep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
          tempDeltaR=999.0;
          tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),thisLep.Eta(),thisLep.Phi());
          if (verbose) cout<<" jet DeltaR between Lep"<<i<<" : "<<tempDeltaR;
          if (tempDeltaR<0.4) {
              isclean_H4l = false;
          }
      }

      // check overlap with fsr photons
      unsigned int N = fsrPhotons_pt.size();
      for(unsigned int i=0; i<N; i++) {
          // don't clean jet from fsr if the photon wasn't matched to tight Id and Isolated lepton
          if (!lep_tightId[fsrPhotons_lepindex[i]]) continue;
          double RelIsoNoFSR=lep_RelIsoNoFSR[fsrPhotons_lepindex[i]];
          if (RelIsoNoFSR>((abs(lep_id[fsrPhotons_lepindex[i]])==11) ? isoCutEl : isoCutMu)) continue;

          TLorentzVector pho;
          pho.SetPtEtaPhiM(fsrPhotons_pt[i],fsrPhotons_eta[i],fsrPhotons_phi[i],0.0);
          tempDeltaR=999.0;
          tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),pho.Eta(),pho.Phi());
          if (verbose) cout<<" jet DeltaR between pho"<<i<<" : "<<tempDeltaR;
          if (tempDeltaR<0.4) {
              isclean_H4l = false;
          }
      }
      if (N>0&&verbose) cout<<endl;

      if(isclean_H4l){
        //JER from database
        JME::JetParameters parameters;
        parameters.setJetPt(goodJets[k].pt());
        parameters.setJetEta(goodJets[k].eta());
        parameters.setRho(muRho);
        float relpterr = resolution_pt.getResolution(parameters);
        float phierr = resolution_phi.getResolution(parameters);

        double jercorr = 1.0; double jercorrup = 1.0; double jercorrdn = 1.0;

        if (isMC && doJER) {
            JME::JetParameters sf_parameters = {{JME::Binning::JetPt, goodJets[k].pt()}, {JME::Binning::JetEta, goodJets[k].eta()}, {JME::Binning::Rho, muRho}};
            float factor = resolution_sf.getScaleFactor(sf_parameters);
            float factorup = resolution_sf.getScaleFactor(sf_parameters, Variation::UP);
            float factordn = resolution_sf.getScaleFactor(sf_parameters, Variation::DOWN);

            double pt_jer, pt_jerup, pt_jerdn;
            const reco::GenJet * genJet = goodJets[k].genJet();
            if (genJet && deltaR(goodJets[k].eta(),goodJets[k].phi(),genJet->eta(),genJet->phi())<0.2
                && (abs(goodJets[k].pt()-genJet->pt())<3*relpterr*goodJets[k].pt())) {
                double gen_pt = genJet->pt();
                pt_jer = max(0.0,gen_pt+factor*(goodJets[k].pt()-gen_pt));
                pt_jerup = max(0.0,gen_pt+factorup*(goodJets[k].pt()-gen_pt));
                pt_jerdn = max(0.0,gen_pt+factordn*(goodJets[k].pt()-gen_pt));
            } else {
                TRandom3 rand;
                rand.SetSeed(abs(static_cast<int>(sin(goodJets[k].phi())*100000)));
                float smear = rand.Gaus(0,1.0);
                float sigma = sqrt(factor*factor-1.0)*relpterr*goodJets[k].pt();
                float sigmaup = sqrt(factorup*factorup-1.0)*relpterr*goodJets[k].pt();
                float sigmadn = sqrt(factordn*factordn-1.0)*relpterr*goodJets[k].pt();
                pt_jer = max(0.0,smear*sigma+goodJets[k].pt());
                pt_jerup = max(0.0,smear*sigmaup+goodJets[k].pt());
                pt_jerdn = max(0.0,smear*sigmadn+goodJets[k].pt());
            }

            jercorr = pt_jer/goodJets[k].pt();
            jercorrup = pt_jerup/goodJets[k].pt();
            jercorrdn = pt_jerdn/goodJets[k].pt();
        }

        TLorentzVector *jet_jer = new TLorentzVector(jercorr*goodJets[k].px(),jercorr*goodJets[k].py(),jercorr*goodJets[k].pz(),jercorr*goodJets[k].energy());
        TLorentzVector *jet_jerup = new TLorentzVector(jercorrup*goodJets[k].px(),jercorrup*goodJets[k].py(),jercorrup*goodJets[k].pz(),jercorrup*goodJets[k].energy());
        TLorentzVector *jet_jerdn = new TLorentzVector(jercorrdn*goodJets[k].px(),jercorrdn*goodJets[k].py(),jercorrdn*goodJets[k].pz(),jercorrdn*goodJets[k].energy());

        bool passPU_;
        if (doJEC && (year==2017 || year==2018 )) {
            passPU_ = bool(goodJets[k].userInt("pileupJetIdUpdated:fullId") & (1 << 0));
        } else if (doJEC && year==2016) {
            passPU_ = bool(goodJets[k].userInt("pileupJetIdUpdated:fullId") & (1 << 2));
        } else {
            passPU_ = bool(goodJets[k].userInt("pileupJetId:fullId") & (1 << 0));
        }
        //if(!(passPU_ || !doPUJetID || jet_jer->Pt()>50)) continue;
        if(!(passPU_ || !doPUJetID || goodJets[k].pt()>50)) continue;

        if(jet_jer->Pt()<10) continue;
        if (verbose) cout<<"Jet nominal: "<<goodJets[k].pt()<<" JER corrected: "<<jet_jer->Pt()<<" JER up: "<<jet_jerup->Pt()<<" JER dn: "<<jet_jerdn->Pt()<<" check Delta between jet and lep / pho: "<<isclean_H4l<<std::endl;

        jecunc->setJetPt(jet_jer->Pt());
        jecunc->setJetEta(goodJets[k].eta());
        double jecunc_up = 1.0+jecunc->getUncertainty(true);
        jecunc->setJetPt(jet_jer->Pt());
        jecunc->setJetEta(goodJets[k].eta());
        double jecunc_dn = 1.0-jecunc->getUncertainty(false);

        if (jet_jer->Pt() > 30.0 && fabs(goodJets[k].eta())<4.7) {
            if (isclean_H4l) {
                pat::Jet goodJets_JECJER_pt30_eta4p7_tmp=goodJets[k];
                goodJets_JECJER_pt30_eta4p7_tmp.setP4(reco::Particle::PolarLorentzVector(jet_jer->Pt(), jet_jer->Eta(), jet_jer->Phi(), jet_jer->M()));
                if (verbose)    std::cout<<"goodJets "<<k<<" : pt: "<<goodJets[k].pt()<<"; jer_pt: "<<jet_jer->Pt()<<"; after setP4 pt: "<<goodJets_JECJER_pt30_eta4p7_tmp.pt()<<std::endl;
                goodJets_JECJER_pt30_eta4p7.push_back(goodJets_JECJER_pt30_eta4p7_tmp);
                njets_pt30_eta4p7++;
                jet_iscleanH4l.push_back((int)jet_pt.size());
            if (jet_jer->Pt() > pt_leadingjet_pt30_eta4p7) {
                    pt_leadingjet_pt30_eta4p7 = jet_jer->Pt();
                    absrapidity_leadingjet_pt30_eta4p7 = jet_jer->Rapidity(); //take abs later
            }
                if (jet_jer->Pt() > jet1pt )  {
                jet2pt=jet1pt; jet2index=jet1index;
                  jet1pt=jet_jer->Pt(); jet1index=(int)jet_pt.size();
                } else if (jet_jer->Pt()>jet2pt) {
                    jet2pt=jet_jer->Pt(); jet2index=(int)jet_pt.size();
                }
                if (fabs(goodJets[k].eta())<2.5) {
                    njets_pt30_eta2p5++;
                if (jet_jer->Pt() > pt_leadingjet_pt30_eta2p5) {
                        pt_leadingjet_pt30_eta2p5 = jet_jer->Pt();
                    }
                    if (jet_jer->Pt() > jet1pt2p5) {
                  jet2pt2p5 = jet1pt2p5; jet2index_2p5 = jet1index_2p5;
                        jet1pt2p5 = jet_jer->Pt(); jet1index_2p5 = (int)jet_pt.size();
                    } else if (jet_jer->Pt() > jet2pt2p5) {
                  jet2pt2p5 = jet_jer->Pt(); jet2index_2p5 = (int)jet_pt.size();
                }
                }
            }  // isclean_H4l
            jet_pt.push_back(jet_jer->Pt());
            jet_pt_raw.push_back(goodJets[k].correctedJet("Uncorrected").pt());///jet Pt without JEC applied
            jet_eta.push_back(jet_jer->Eta());
            jet_phi.push_back(jet_jer->Phi());
            jet_mass.push_back(jet_jer->M());
            if (doJEC && (year==2017 || year==2018 || year==2016)) {    //FIXME for UL
                jet_pumva.push_back(goodJets[k].userFloat("pileupJetIdUpdated:fullDiscriminant"));
            } else {
                jet_pumva.push_back(goodJets[k].userFloat("pileupJetId:fullDiscriminant"));
            }
            //jet_csvv2.push_back(goodJets[k].bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll"));
            jet_csvv2.push_back(goodJets[k].bDiscriminator("pfDeepCSVJetTags:probb")+goodJets[k].bDiscriminator("pfDeepCSVJetTags:probbb"));
            jet_csvv2_.push_back(goodJets[k].bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll"));

            TRandom3 rand;
            rand.SetSeed(abs(static_cast<int>(sin(goodJets[k].phi())*100000)));
            float coin = rand.Uniform(1.);

            double jet_scalefactor = 1.0;

            jet_scalefactor = reader->eval_auto_bounds(
                "central",
                BTagEntry::FLAV_B,
                goodJets[k].eta(),
                goodJets[k].pt()
                );

            if ((goodJets[k].bDiscriminator("pfDeepCSVJetTags:probb")+goodJets[k].bDiscriminator("pfDeepCSVJetTags:probbb"))>BTagCut && coin>(1.0-jet_scalefactor))
            //if (goodJets[k].bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll")>BTagCut && coin>(1.0-jet_scalefactor))
            {
                jet_isbtag.push_back(1);
            } else {
                jet_isbtag.push_back(0);
            }

            if ((goodJets[k].bDiscriminator("pfDeepCSVJetTags:probb")+goodJets[k].bDiscriminator("pfDeepCSVJetTags:probbb"))>BTagCut && isclean_H4l) nbjets_pt30_eta4p7++;
            //if (goodJets[k].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")>BTagCut && isclean_H4l) nbjets_pt30_eta4p7++;
            if (isMC) {
                jet_hadronFlavour.push_back(goodJets[k].hadronFlavour());
                jet_partonFlavour.push_back(goodJets[k].partonFlavour());
            } else {
                jet_hadronFlavour.push_back(-1);
                jet_partonFlavour.push_back(-1);
            }
            jet_QGTagger.push_back(goodJetQGTagger[k]);
            jet_axis2.push_back(goodJetaxis2[k]);
            jet_ptD.push_back(goodJetptD[k]);
            jet_mult.push_back(goodJetmult[k]);
            jet_relpterr.push_back(relpterr);
            jet_phierr.push_back(phierr);
            jet_bTagEffi.push_back(helper.get_bTagEffi(jet_jer->Pt(), jet_jer->Eta(), hbTagEffi));
            jet_cTagEffi.push_back(helper.get_bTagEffi(jet_jer->Pt(), jet_jer->Eta(), hcTagEffi));
            jet_udsgTagEffi.push_back(helper.get_bTagEffi(jet_jer->Pt(), jet_jer->Eta(), hudsgTagEffi));
        }   // if (jet_jer->Pt() > 30.0 && fabs(goodJets[k].eta())<4.7)

        // JER up
        if (jet_jerup->Pt() > 30.0 && fabs(jet_jerup->Eta())<4.7) {
            if (isclean_H4l) {
                njets_pt30_eta4p7_jerup++;
                jet_jerup_iscleanH4l.push_back((int)jet_jerup_pt.size());
                if (jet_jerup->Pt() > pt_leadingjet_pt30_eta4p7_jerup) {
                    pt_leadingjet_pt30_eta4p7_jerup = jet_jerup->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jerup = jet_jerup->Rapidity(); //take abs later
                    //jet2pt_jesup=jet1pt_jesup; jet2index_jesup=jet1index_jesup;
                    //jet1pt_jesup=jet_jerup->Pt(); jet1index_jesup=(int)jet_jerup_pt.size();
                } //else if (jet_jerup->Pt() > jet2pt_jesup) {
                    //jet2pt_jesup=jet_jerup->Pt(); jet2index_jesup=(int)jet_jerup_pt.size();
                //}
                if (fabs(jet_jerup->Eta())<2.5) {
                    njets_pt30_eta2p5_jerup++;
                    if (jet_jerup->Pt() > pt_leadingjet_pt30_eta2p5_jerup) {
                        pt_leadingjet_pt30_eta2p5_jerup = jet_jerup->Pt();
                        //jet2pt2p5_jesup=jet1pt2p5_jesup; jet2index2p5_jesup=jet1index2p5_jesup;
                        //jet1pt2p5_jesup=jet_jerup->Pt(); jet1index2p5_jesup=(int)jet_jerup_pt.size();
                    } //else if (jet_jerup->Pt()>jet2pt2p5_jesup) {
                        //jet2pt2p5_jesup=jet_jerup->Pt(); jet2index2p5_jesup=(int)jet_jerup_pt.size();
                    //}
                }
            }
            jet_jerup_pt.push_back(jet_jerup->Pt());
            jet_jerup_eta.push_back(jet_jerup->Eta());
            jet_jerup_phi.push_back(jet_jerup->Phi());
            jet_jerup_mass.push_back(jet_jerup->M());
        }

        // JER dn
        if (jet_jerdn->Pt() > 30.0 && fabs(jet_jerdn->Eta())<4.7) {
            if (isclean_H4l) {
                njets_pt30_eta4p7_jerdn++;
                jet_jerdn_iscleanH4l.push_back((int)jet_jerdn_pt.size());
                if (jet_jerdn->Pt() > pt_leadingjet_pt30_eta4p7_jerdn) {
                    pt_leadingjet_pt30_eta4p7_jerdn = jet_jerdn->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jerdn = jet_jerdn->Rapidity(); //take abs later
                    //jet2pt_jesdn=jet1pt_jesdn; jet2index_jesdn=jet1index_jesdn;
                    //jet1pt_jesdn=jet_jerdn->Pt(); jet1index_jesdn=(int)jet_jerdn_pt.size();
                } //else if (jet_jerdn->Pt()>jet2pt_jesdn) {
                    //jet2pt_jesdn=jet_jerdn->Pt(); jet2index_jesdn=(int)jet_jerdn_pt.size();
                //}
                if (fabs(jet_jerdn->Eta())<2.5) {
                    njets_pt30_eta2p5_jerdn++;
                    if (jet_jerdn->Pt() > pt_leadingjet_pt30_eta2p5_jerdn) {
                        pt_leadingjet_pt30_eta2p5_jerdn = jet_jerdn->Pt();
                        //jet2pt2p5_jesdn=jet1pt2p5_jesdn; jet2index2p5_jesdn=jet1index2p5_jesdn;
                        //jet1pt2p5_jesdn=jet_jerdn->Pt(); jet1index2p5_jesdn=(int)jet_jerdn_pt.size();
                    } //else if (jet_jerdn->Pt()>jet2pt2p5_jesdn) {
                        //jet2pt2p5_jesdn=jet_jerdn->Pt(); jet2index2p5_jesdn=(int)jet_jerdn_pt.size();
                    //}
                }
            }
            jet_jerdn_pt.push_back(jet_jerdn->Pt());
            jet_jerdn_eta.push_back(jet_jerdn->Eta());
            jet_jerdn_phi.push_back(jet_jerdn->Phi());
            jet_jerdn_mass.push_back(jet_jerdn->M());
        }

        double jetPx_jesup = jecunc_up * jet_jer->Px();
        double jetPy_jesup = jecunc_up * jet_jer->Py();
        double jetPz_jesup = jecunc_up * jet_jer->Pz();
        double jetE_jesup = sqrt(jetPx_jesup*jetPx_jesup + jetPy_jesup*jetPy_jesup + jetPz_jesup*jetPz_jesup + jet_jer->M()*jet_jer->M());
        TLorentzVector *jet_jesup = new TLorentzVector(jetPx_jesup,jetPy_jesup,jetPz_jesup,jetE_jesup);

        if (jet_jesup->Pt() > 30.0 && fabs(jet_jesup->Eta())<4.7) {
            if (isclean_H4l) {
                njets_pt30_eta4p7_jesup++;
                jet_jesup_iscleanH4l.push_back((int)jet_jesup_pt.size());
                if (jet_jesup->Pt() > pt_leadingjet_pt30_eta4p7_jesup) {
                    pt_leadingjet_pt30_eta4p7_jesup = jet_jesup->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jesup = jet_jesup->Rapidity(); //take abs later
                    //jet2pt_jesup=jet1pt_jesup; jet2index_jesup=jet1index_jesup;
                    //jet1pt_jesup=jet_jesup->Pt(); jet1index_jesup=(int)jet_jesup_pt.size();
                }  //else if (jet_jesup->Pt()>jet2pt_jesup) {
                    //jet2pt_jesup=jet_jesup->Pt(); jet2index_jesup=(int)jet_jesup_pt.size();
                //}
                if (fabs(jet_jesup->Eta())<2.5) {
                    njets_pt30_eta2p5_jesup++;
                    if (jet_jesup->Pt() > pt_leadingjet_pt30_eta2p5_jesup) {
                        pt_leadingjet_pt30_eta2p5_jesup = jet_jesup->Pt();
                        //jet2pt2p5_jesup=jet1pt2p5_jesup; jet2index2p5_jesup=jet1index2p5_jesup;
                        //jet1pt2p5_jesup=jet_jesup->Pt(); jet1index2p5_jesup=(int)jet_jesup_pt.size();
                    } //else if (jet_jesup->Pt() > jet2pt2p5_jesup) {
                        //jet2pt2p5_jesup=jet_jesup->Pt(); jet2index2p5_jesup=(int)jet_jesup_pt.size();
                    //}
                }
            }
            TLorentzVector jet_jesup(jetPx_jesup,jetPy_jesup,jetPz_jesup,jetE_jesup);
            jet_jesup_pt.push_back(jet_jesup.Pt());
            jet_jesup_eta.push_back(jet_jesup.Eta());
            jet_jesup_phi.push_back(jet_jesup.Phi());
            jet_jesup_mass.push_back(jet_jesup.M());
            jet_QGTagger_jesup.push_back(goodJetQGTagger[k]);
        }

        double jetPx_jesdn = jecunc_dn * jet_jer->Px();
        double jetPy_jesdn = jecunc_dn * jet_jer->Py();
        double jetPz_jesdn = jecunc_dn * jet_jer->Pz();
        double jetE_jesdn = sqrt(jetPx_jesdn*jetPx_jesdn + jetPy_jesdn*jetPy_jesdn + jetPz_jesdn*jetPz_jesdn + jet_jer->M()*jet_jer->M());
        TLorentzVector *jet_jesdn = new TLorentzVector(jetPx_jesdn,jetPy_jesdn,jetPz_jesdn,jetE_jesdn);

        if (jet_jesdn->Pt() > 30.0 && fabs(jet_jesdn->Eta())<4.7) {
            if (isclean_H4l) {
                njets_pt30_eta4p7_jesdn++;
                jet_jesdn_iscleanH4l.push_back(jet_jesdn_pt.size());
                if (jet_jesdn->Pt() > pt_leadingjet_pt30_eta4p7_jesdn) {
                    pt_leadingjet_pt30_eta4p7_jesdn = jet_jesdn->Pt();
                    absrapidity_leadingjet_pt30_eta4p7_jesdn = jet_jesdn->Rapidity(); //take abs later
                    //jet2pt_jesdn=jet1pt_jesdn; jet2index_jesdn=jet1index_jesdn;
                    //jet1pt_jesdn=jet_jesdn->Pt(); jet1index_jesdn=jet_jesdn_pt.size();
                } //else if (jet_jesdn->Pt()>jet2pt_jesdn) {
                    //jet2pt_jesdn=jet_jesdn->Pt(); jet2index_jesdn=jet_jesdn_pt.size();
                //}
                if (fabs(jet_jesdn->Eta())<2.5) {
                    njets_pt30_eta2p5_jesdn++;
                    if (jet_jesdn->Pt() > pt_leadingjet_pt30_eta2p5_jesdn) {
                        pt_leadingjet_pt30_eta2p5_jesdn = jet_jesdn->Pt();
                        //jet2pt2p5_jesdn=jet1pt2p5_jesdn; jet2index2p5_jesdn=jet1index2p5_jesdn;
                        //jet1pt2p5_jesdn=jet_jesdn->Pt(); jet1index2p5_jesdn=jet_jesdn_pt.size();
                    } //else if (jet_jesdn->Pt()>jet2pt2p5_jesdn) {
                        //jet2pt2p5_jesdn=jet_jesdn->Pt(); jet2index2p5_jesdn=jet_jesdn_pt.size();
                    //}
                }
            }
            TLorentzVector jet_jesdn(jetPx_jesdn,jetPy_jesdn,jetPz_jesdn,jetE_jesdn);
            jet_jesdn_pt.push_back(jet_jesdn.Pt());
            jet_jesdn_eta.push_back(jet_jesdn.Eta());
            jet_jesdn_phi.push_back(jet_jesdn.Phi());
            jet_jesdn_mass.push_back(jet_jesdn.M());
            jet_QGTagger_jesdn.push_back(goodJetQGTagger[k]);

        }

      }//jet are clean from tight leptons and fsr photons

  } // loop over jets

  // merged jet
  for( unsigned int k = 0; k < selectedMergedJets.size(); k++) {

    double tempDeltaR = 999.0;
    bool isclean_H4l = true;

    unsigned int Nleptons = lep_pt.size();
    for (unsigned int i=0; i<Nleptons; i++) {

      if (abs(lep_id[i])==13 && lep_RelIsoNoFSR[i]>isoCutMu) continue;
      if (abs(lep_id[i])==11 && lep_RelIsoNoFSR[i]>isoCutEl) continue;
      if (!(lep_tightId[i])) continue;
      TLorentzVector thisLep;
      thisLep.SetPtEtaPhiM(lep_pt[i],lep_eta[i],lep_phi[i],lep_mass[i]);
      tempDeltaR=999.0;
      tempDeltaR=deltaR(selectedMergedJets[k].eta(),selectedMergedJets[k].phi(),thisLep.Eta(),thisLep.Phi());
      if (tempDeltaR<0.8) {
        isclean_H4l = false;
      }
    }

    // check overlap with fsr photons
    unsigned int N = fsrPhotons_pt.size();
    for(unsigned int i=0; i<N; i++) {

        // don't clean jet from fsr if the photon wasn't matched to tight Id and Isolated lepton
        if (!lep_tightId[fsrPhotons_lepindex[i]]) continue;
        double RelIsoNoFSR=lep_RelIsoNoFSR[fsrPhotons_lepindex[i]];
        if (RelIsoNoFSR>((abs(lep_id[fsrPhotons_lepindex[i]])==11) ? isoCutEl : isoCutMu)) continue;

        TLorentzVector pho;
        pho.SetPtEtaPhiM(fsrPhotons_pt[i],fsrPhotons_eta[i],fsrPhotons_phi[i],0.0);
        tempDeltaR=999.0;
        tempDeltaR=deltaR(goodJets[k].eta(),goodJets[k].phi(),pho.Eta(),pho.Phi());
        if (verbose) cout<<" jet DeltaR between pho"<<i<<" : "<<tempDeltaR;
        if (tempDeltaR<0.8) {
            isclean_H4l = false;
        }
    }

    if (isclean_H4l){

      //JER from database
      JME::JetParameters ak8_parameters;
      ak8_parameters.setJetPt(selectedMergedJets[k].pt());
      ak8_parameters.setJetEta(selectedMergedJets[k].eta());
      ak8_parameters.setRho(muRho);
      float reak8pterr = ak8_resolution_pt.getResolution(ak8_parameters);
      float phiak8err = ak8_resolution_phi.getResolution(ak8_parameters);
      double ak8_jercorr = 1.0; double ak8_jercorrup = 1.0; double ak8_jercorrdn = 1.0;

      if(isMC && doJER){

        JME::JetParameters sf_ak8parameters = {{JME::Binning::JetPt, selectedMergedJets[k].pt()}, {JME::Binning::JetEta, selectedMergedJets[k].eta()}, {JME::Binning::Rho, muRho}};
        float ak8_factor = ak8_resolution_sf.getScaleFactor(sf_ak8parameters);
        float ak8_factorup = ak8_resolution_sf.getScaleFactor(sf_ak8parameters, Variation::UP);
        float ak8_factordonw = ak8_resolution_sf.getScaleFactor(sf_ak8parameters, Variation::DOWN);

        double ak8_ptjer, ak8_ptjer_up, ak8_ptjer_donw;
        const reco::GenJet * ak8_genJet = selectedMergedJets[k].genJet();
        if(ak8_genJet && deltaR(selectedMergedJets[k].eta(),selectedMergedJets[k].phi(),ak8_genJet->eta(),ak8_genJet->phi())<0.4){//match to gen ak8 jet
          double ak8_gen_pt = ak8_genJet->pt();
          ak8_ptjer = max(0.0,ak8_gen_pt+ak8_factor*(selectedMergedJets[k].pt()-ak8_gen_pt));
          ak8_ptjer_up = max(0.0,ak8_gen_pt+ak8_factorup*(selectedMergedJets[k].pt()-ak8_gen_pt));
          ak8_ptjer_donw = max(0.0,ak8_gen_pt+ak8_factordonw*(selectedMergedJets[k].pt()-ak8_gen_pt));
        } else{
          TRandom3 rand;
          rand.SetSeed(abs(static_cast<int>(sin(selectedMergedJets[k].phi())*100000)));
          float ak8_smear = rand.Gaus(0,1.0);
          float ak8_sigma = sqrt(ak8_factor*ak8_factor-1.0)*reak8pterr*selectedMergedJets[k].pt();
          float ak8_sigmaup = sqrt(ak8_factorup*ak8_factorup-1.0)*reak8pterr*selectedMergedJets[k].pt();
          float ak8_sigmadn = sqrt(ak8_factordonw*ak8_factordonw-1.0)*reak8pterr*selectedMergedJets[k].pt();
          ak8_ptjer = max(0.0,ak8_smear*ak8_sigma+selectedMergedJets[k].pt());
          ak8_ptjer_up = max(0.0,ak8_smear*ak8_sigmaup+selectedMergedJets[k].pt());
          ak8_ptjer_donw = max(0.0,ak8_smear*ak8_sigmadn+selectedMergedJets[k].pt());
        }

        mergedjet_jer.push_back(ak8_ptjer/selectedMergedJets[k].pt());
        ak8_jercorr = ak8_ptjer/selectedMergedJets[k].pt();
        ak8_jercorrup = ak8_ptjer_up/selectedMergedJets[k].pt();
        ak8_jercorrdn = ak8_ptjer_donw/selectedMergedJets[k].pt();

      }

      TLorentzVector *ak8_jer = new TLorentzVector(ak8_jercorr*selectedMergedJets[k].px(),ak8_jercorr*selectedMergedJets[k].py(),ak8_jercorr*selectedMergedJets[k].pz(),ak8_jercorr*selectedMergedJets[k].energy());
      TLorentzVector *ak8_jerup = new TLorentzVector(ak8_jercorrup*selectedMergedJets[k].px(),ak8_jercorrup*selectedMergedJets[k].py(),ak8_jercorrup*selectedMergedJets[k].pz(),ak8_jercorrup*selectedMergedJets[k].energy());
      TLorentzVector *ak8_jerdonw = new TLorentzVector(ak8_jercorrdn*selectedMergedJets[k].px(),ak8_jercorrdn*selectedMergedJets[k].py(),ak8_jercorrdn*selectedMergedJets[k].pz(),ak8_jercorrdn*selectedMergedJets[k].energy());

      //jer center
      mergedjet_jer_pt.push_back(ak8_jer->Pt());
      mergedjet_jer_eta.push_back(ak8_jer->Eta());
      mergedjet_jer_phi.push_back(ak8_jer->Phi());
      mergedjet_jer_mass.push_back(ak8_jer->M());
      //jer up
      mergedjet_jerup_pt.push_back(ak8_jerup->Pt());
      mergedjet_jerup_eta.push_back(ak8_jerup->Eta());
      mergedjet_jerup_phi.push_back(ak8_jerup->Phi());
      mergedjet_jerup_mass.push_back(ak8_jerup->M());
      //jer down
      mergedjet_jerdn_pt.push_back(ak8_jerdonw->Pt());
      mergedjet_jerdn_eta.push_back(ak8_jerdonw->Eta());
      mergedjet_jerdn_phi.push_back(ak8_jerdonw->Phi());
      mergedjet_jerdn_mass.push_back(ak8_jerdonw->M());

      //jer error
      mergedjet_jer_pterr.push_back(reak8pterr);
      mergedjet_jer_phierr.push_back(phiak8err);



      mergedjet_iscleanH4l.push_back((int)mergedjet_pt.size());
      mergedjet_pt.push_back((float)selectedMergedJets[k].pt());
      mergedjet_ptrow.push_back((float)selectedMergedJets[k].correctedJet("Uncorrected").pt());
      mergedjet_eta.push_back((float)selectedMergedJets[k].eta());
      mergedjet_phi.push_back((float)selectedMergedJets[k].phi());
      mergedjet_mass.push_back((float)selectedMergedJets[k].mass());
      //mergedjet_L1.push_back((float)selectedMergedJets[k].jecFactor("L1FastJet")); // current JEC to L1

      //mergedjet_softdropmass.push_back((float)selectedMergedJets[k].userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass"));
      mergedjet_softdropmass.push_back((float)selectedMergedJets[k].userFloat("ak8PFJetsPuppiSoftDropMass"));
      //mergedjet_prunedmass.push_back((float)selectedMergedJets[k].userFloat("ak8PFJetsCHSCorrPrunedMass"));
      mergedjet_tau1.push_back((float)selectedMergedJets[k].userFloat("NjettinessAK8Puppi:tau1") );
      mergedjet_tau2.push_back((float)selectedMergedJets[k].userFloat("NjettinessAK8Puppi:tau2") );
      mergedjet_btag.push_back((float)selectedMergedJets[k].bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") );

      mergedjet_ZvsQCD.push_back((float)selectedMergedJets[k].bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD"));
      //std::cout<<"(float)selectedMergedJets[k].bDiscriminator(pfDeepBoostedDiscriminatorsJetTags:ZvsQCD) = " <<(float)selectedMergedJets[k].bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZvsQCD")<<std::endl;
      mergedjet_ZbbvsQCD.push_back((float)selectedMergedJets[k].bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"));
      mergedjet_WvsQCD.push_back((float)selectedMergedJets[k].bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:WvsQCD"));
      mergedjet_ZHbbvsQCD.push_back((float)selectedMergedJets[k].bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD"));
      mergedjet_HbbvsQCD.push_back((float)selectedMergedJets[k].bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:HbbvsQCD"));
      mergedjet_H4qvsQCD.push_back((float)selectedMergedJets[k].bDiscriminator("pfDeepBoostedDiscriminatorsJetTags:H4qvsQCD"));
      mergedjet_ZvsQCD_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD"));
      mergedjet_ZbbvsQCD_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZbbvsQCD"));
      mergedjet_WvsQCD_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"));
      mergedjet_ZHbbvsQCD_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD"));
      mergedjet_HbbvsQCD_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:HbbvsQCD"));
      mergedjet_H4qvsQCD_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:H4qvsQCD"));

      mergedjet_Net_Xbb_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probXbb"));
      mergedjet_Net_Xcc_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probXcc"));
      mergedjet_Net_Xqq_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probXqq"));
      mergedjet_Net_QCDbb_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDbb"));
      mergedjet_Net_QCDcc_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDcc"));
      mergedjet_Net_QCDother_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDothers"));
      mergedjet_Net_QCDb_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDb"));
      mergedjet_Net_QCDc_de.push_back((float)selectedMergedJets[k].bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDc"));

      mergedjet_nbHadrons.push_back((float)selectedMergedJets[k].jetFlavourInfo().getbHadrons().size());
      mergedjet_ncHadrons.push_back((float)selectedMergedJets[k].jetFlavourInfo().getcHadrons().size());
      mergedjet_partonFlavour.push_back((float)selectedMergedJets[k].partonFlavour());
      mergedjet_hadronFlavour.push_back((float)selectedMergedJets[k].hadronFlavour());

      if (verbose) cout<<"double btag: "<<selectedMergedJets[k].bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags")<<endl;

      auto wSubjets = selectedMergedJets[k].subjets("SoftDropPuppi");
      int nsub = 0;
      math::XYZTLorentzVector fatJet;
      TLorentzVector fatjetUncorr, fatjetUncorr_subjet;
      vector<float> subjets_pt, subjets_eta, subjets_phi, subjets_mass, subjets_btag, subjetsUncorr_pt, subjetsUncorr_eta, subjetsUncorr_phi, subjetsUncorr_mass;
      vector<int> subjets_hadronFlavour, subjets_partonFlavour;
      subjets_pt.clear(); subjets_eta.clear(); subjets_phi.clear(); subjets_mass.clear(); subjets_btag.clear();
      subjets_hadronFlavour.clear(); subjets_partonFlavour.clear();
      for ( auto const & iw : wSubjets ) {
        nsub = nsub + 1;
        fatJet = fatJet + iw->p4();
        fatjetUncorr.SetPtEtaPhiM(iw->correctedP4(0).pt(),iw->correctedP4(0).eta(),iw->correctedP4(0).phi(),iw->correctedP4(0).mass());
        fatjetUncorr_subjet += fatjetUncorr;
        subjets_pt.push_back((float)iw->pt());
        subjets_eta.push_back((float)iw->eta());
        subjets_phi.push_back((float)iw->phi());
        subjets_mass.push_back((float)iw->mass());
        subjetsUncorr_pt.push_back((float)iw->correctedP4(0).pt());
        subjetsUncorr_eta.push_back((float)iw->correctedP4(0).eta());
        subjetsUncorr_phi.push_back((float)iw->correctedP4(0).phi());
        subjetsUncorr_mass.push_back((float)iw->correctedP4(0).mass());
        subjets_btag.push_back((float)iw->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
        subjets_hadronFlavour.push_back(iw->hadronFlavour());
        subjets_partonFlavour.push_back(iw->partonFlavour());
        if (verbose) cout<<"subjet parton "<<iw->partonFlavour()<<" hadron "<<iw->hadronFlavour()<<endl;

      }

      mergedjet_nsubjet.push_back(nsub);
      mergedjet_subjet_pt.push_back(subjets_pt);
      mergedjet_subjet_eta.push_back(subjets_eta);
      mergedjet_subjet_phi.push_back(subjets_phi);
      mergedjet_subjet_mass.push_back(subjets_mass);
      mergedjet_subjet_btag.push_back(subjets_btag);
      mergedjet_subjetUncorr_pt.push_back(subjetsUncorr_pt);
      mergedjet_subjetUncorr_eta.push_back(subjetsUncorr_eta);
      mergedjet_subjetUncorr_phi.push_back(subjetsUncorr_phi);
      mergedjet_subjetUncorr_mass.push_back(subjetsUncorr_mass);
      mergedjet_subjet_softDropMass.push_back(fatJet.M());
      mergedjet_subjet_softDropMassUncorr.push_back(fatjetUncorr_subjet.M());
      mergedjet_subjet_partonFlavour.push_back(subjets_partonFlavour);
      mergedjet_subjet_hadronFlavour.push_back(subjets_hadronFlavour);
      //cout<<"sdmass with CMSSW = "<<mergedjet_softdropmass[k]<<endl;
      //cout<<"sdmass with older jec = "<<mergedjet_subjet_softDropMass[k]<<endl;
      //cout<<"sdmass without jer  = "<<mergedjet_subjet_softDropMassUncorr[k]<<endl;
      //cout<<"sdmass with jer = "<<mergedjet_subjet_softDropMassUncorr[k]*mergedjet_jer[k]<<endl;

    }

  }

}

void UFHZZ4LAna::setGENVariables(edm::Handle<reco::GenParticleCollection> prunedgenParticles,
                                 edm::Handle<edm::View<pat::PackedGenParticle> > packedgenParticles,
                                 edm::Handle<edm::View<reco::GenJet> > genJets)
{
  reco::GenParticleCollection::const_iterator genPart;
  int j = -1;
  int nGENLeptons=0;
  TLorentzVector GENmom1, GENmom2;
  TLorentzVector LS3_Z1_1, LS3_Z1_2, LS3_Z2_1, LS3_Z2_2, GEN_HVec;
  //int GENmom1_id=-999, GENmom2_id=-999;
  int counter_initParticle=0;
  if (verbose) cout<<"begin looping on gen particles"<<endl;
  for(genPart = prunedgenParticles->begin(); genPart != prunedgenParticles->end(); genPart++) {
       j++;

       if (genPart->status() == 21) {
           counter_initParticle++;
           if (counter_initParticle==1){
               GENmom1.SetPxPyPzE(genPart->px(),genPart->py(),genPart->pz(),genPart->energy());
               //GENmom1_id=genPart->pdgId();
           }
           if (counter_initParticle==2){
               GENmom2.SetPxPyPzE(genPart->px(),genPart->py(),genPart->pz(),genPart->energy());
               //GENmom2_id=genPart->pdgId();
           }
       }  //end if status, incoming particles
       if (counter_initParticle > 2)
       {
           cout << "initial particle can't be more than 2... please check for issue" << endl;
           // exit(0);
       }

      if (abs(genPart->pdgId())==11  || abs(genPart->pdgId())==13 || abs(genPart->pdgId())==15) {

          if (!(genPart->status()==1 || abs(genPart->pdgId())==15)) continue;
          if (!(genAna.MotherID(&prunedgenParticles->at(j))==23 || genAna.MotherID(&prunedgenParticles->at(j))==443 || genAna.MotherID(&prunedgenParticles->at(j))==553 || abs(genAna.MotherID(&prunedgenParticles->at(j)))==24) ) continue;

          nGENLeptons++;
          if (verbose) cout<<"found a gen lepton: id "<<genPart->pdgId()<<" pt: "<<genPart->pt()<<" eta: "<<genPart->eta()<<" status: "<<genPart->status()<<endl;

          // Collect FSR photons
          TLorentzVector lep_dressed;
          lep_dressed.SetPtEtaPhiE(genPart->pt(),genPart->eta(),genPart->phi(),genPart->energy());
          set<int> gen_fsrset;
          for(size_t k=0; k<packedgenParticles->size();k++){
              if( (*packedgenParticles)[k].status() != 1) continue; // stable particles only
              if( (*packedgenParticles)[k].pdgId() != 22) continue; // only photons
              double this_dR_lgamma = deltaR(genPart->eta(), genPart->phi(), (*packedgenParticles)[k].eta(), (*packedgenParticles)[k].phi());
              bool idmatch=false;
              if ((*packedgenParticles)[k].mother(0)->pdgId()==genPart->pdgId() ) idmatch=true;
              const reco::Candidate * mother = (*packedgenParticles)[k].mother(0);
              for(size_t m=0;m<mother->numberOfMothers();m++) {
                  if ( (*packedgenParticles)[k].mother(m)->pdgId() == genPart->pdgId() ) idmatch=true;
              }
              if (!idmatch) continue;
              if(this_dR_lgamma<((abs(genPart->pdgId())==11)?genIsoConeSizeEl:genIsoConeSizeMu)) {
                  gen_fsrset.insert(k);
                  TLorentzVector gamma;
                  gamma.SetPtEtaPhiE((*packedgenParticles)[k].pt(),(*packedgenParticles)[k].eta(),(*packedgenParticles)[k].phi(),(*packedgenParticles)[k].energy());
                  lep_dressed = lep_dressed+gamma;
              }
          } // Dressed leptons loop
          if (verbose) cout<<"gen lep pt "<<genPart->pt()<< " dressed pt: " << lep_dressed.Pt()<<endl;

          GENlep_id.push_back( genPart->pdgId() );
          GENlep_status.push_back(genPart->status());
          GENlep_pt.push_back( lep_dressed.Pt() );
          GENlep_eta.push_back( lep_dressed.Eta() );
          GENlep_phi.push_back( lep_dressed.Phi() );
          GENlep_mass.push_back( lep_dressed.M() );
          GENlep_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
          GENlep_MomMomId.push_back(genAna.MotherMotherID(&prunedgenParticles->at(j)));

          TLorentzVector thisLep;
          thisLep.SetPtEtaPhiM(lep_dressed.Pt(),lep_dressed.Eta(),lep_dressed.Phi(),lep_dressed.M());
          // GEN iso calculation
          double this_GENiso=0.0;
          double this_GENneutraliso=0.0;
          double this_GENchargediso=0.0;
          if (verbose) cout<<"gen iso calculation"<<endl;
          for(size_t k=0; k<packedgenParticles->size();k++){
              if( (*packedgenParticles)[k].status() != 1 ) continue; // stable particles only
              if (abs((*packedgenParticles)[k].pdgId())==12 || abs((*packedgenParticles)[k].pdgId())==14 || abs((*packedgenParticles)[k].pdgId())==16) continue; // exclude neutrinos
              if ((abs((*packedgenParticles)[k].pdgId())==11 || abs((*packedgenParticles)[k].pdgId())==13)) continue; // exclude leptons
              if (gen_fsrset.find(k)!=gen_fsrset.end()) continue; // exclude particles which were selected as fsr photons
              double this_dRvL = deltaR(thisLep.Eta(), thisLep.Phi(), (*packedgenParticles)[k].eta(), (*packedgenParticles)[k].phi());
              if(this_dRvL<((abs(genPart->pdgId())==11)?genIsoConeSizeEl:genIsoConeSizeMu)) {
                  if (verbose) cout<<"adding to geniso id: "<<(*packedgenParticles)[k].pdgId()<<" status: "<<(*packedgenParticles)[k].status()<<" pt: "<<(*packedgenParticles)[k].pt()<<" dR: "<<this_dRvL<<endl;
                  this_GENiso = this_GENiso + (*packedgenParticles)[k].pt();
                  if ((*packedgenParticles)[k].charge()==0) this_GENneutraliso = this_GENneutraliso + (*packedgenParticles)[k].pt();
                  if ((*packedgenParticles)[k].charge()!=0) this_GENchargediso = this_GENchargediso + (*packedgenParticles)[k].pt();
              }
          } // GEN iso loop
          this_GENiso = this_GENiso/thisLep.Pt();
          if (verbose) cout<<"gen lep pt: "<<thisLep.Pt()<<" rel iso: "<<this_GENiso<<endl;
          GENlep_RelIso.push_back(this_GENiso);
          // END GEN iso calculation

      } // leptons

      // Higgs
      vector<int> this_id, this_status;
      vector<double> this_pt, this_eta, this_phi, this_mass;
      //if (genPart->pdgId()==25 && genPart->status()==22) {
      if (genPart->pdgId()==25) {

          GENMH=genPart->mass();
          GENH_pt.push_back(genPart->pt());
          GENH_eta.push_back(genPart->eta());
          GENH_phi.push_back(genPart->phi());
          GENH_mass.push_back(genPart->mass());
          GENH_Momid.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
          GENH_MomMomid.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
          GENH_id.push_back(genPart->pdgId());
          GENH_status.push_back(genPart->status());
          GENH_isHard.push_back((int)genPart->isHardProcess());
          GENH_nDaughters.push_back(genPart->numberOfDaughters());

          int ndau = genPart->numberOfDaughters();
          for(int d=0; d<ndau; d++){
            const reco::Candidate *Hdau=genPart->daughter(d);
            this_id.push_back(Hdau->pdgId());
            this_status.push_back(Hdau->status());
            this_pt.push_back(Hdau->pt());
            this_eta.push_back(Hdau->eta());
            this_phi.push_back(Hdau->phi());
            this_mass.push_back(Hdau->mass());
          }

          GENH_dau_id.push_back(this_id);
          GENH_dau_status.push_back(this_status);
          GENH_dau_pt.push_back(this_pt);
          GENH_dau_eta.push_back(this_eta);
          GENH_dau_phi.push_back(this_phi);
          GENH_dau_mass.push_back(this_mass);
      }

      if(genPart->isHardProcess()){
        GEN_id.push_back(genPart->pdgId());
        GEN_status.push_back(genPart->status());
      }

      //find quark come from Z decay
      if(abs(genPart->pdgId())<=6 && abs(genAna.MotherID(&prunedgenParticles->at(j)))==23 && (genPart->status()==1 || genPart->status()==23)){
        GEN_Zq_pt.push_back(genPart->pt());
        GEN_Zq_eta.push_back(genPart->eta());
        GEN_Zq_phi.push_back(genPart->phi());
        GEN_Zq_mass.push_back(genPart->mass());
        GEN_Zq_id.push_back(genPart->pdgId());
        GEN_Zq_Momid.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
        GEN_Zq_MomMomid.push_back(genAna.MotherMotherID(&prunedgenParticles->at(j)));
        GEN_Zq_isHard.push_back((int)genPart->isHardProcess());
      }

      //VBF quark
      if(abs(genPart->pdgId())<=6 && abs(genAna.MotherID(&prunedgenParticles->at(j)))!=23 && genPart->status()==23 && genPart->isHardProcess()){
      //if(abs(genPart->pdgId())<=6 && abs(genAna.MotherID(&prunedgenParticles->at(j)))!=23 && genPart->isHardProcess()){
        GEN_VBF_pt.push_back(genPart->pt());
        GEN_VBF_eta.push_back(genPart->eta());
        GEN_VBF_phi.push_back(genPart->phi());
        GEN_VBF_mass.push_back(genPart->mass());
        GEN_VBF_id.push_back(genPart->pdgId());
        GEN_VBF_status.push_back(genPart->status());
        GEN_VBF_Momid.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
        GEN_VBF_MomMomid.push_back(genAna.MotherMotherID(&prunedgenParticles->at(j)));
      }

      //find quark not come from Z decay, associated quark candidate
      vector<double> dau_pt, dau_eta, dau_phi, dau_mass;
      vector<int> dau_id, dau_status;
      if(abs(genPart->pdgId())<=6 && abs(genAna.MotherID(&prunedgenParticles->at(j)))!=23){
        GEN_q_id.push_back(genPart->pdgId());
        GEN_q_pt.push_back(genPart->pt());
        GEN_q_eta.push_back(genPart->eta());
        GEN_q_phi.push_back(genPart->phi());
        GEN_q_mass.push_back(genPart->mass());
        GEN_q_status.push_back(genPart->status());
        GEN_q_Momid.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
        GEN_q_MomMomid.push_back(genAna.MotherMotherID(&prunedgenParticles->at(j)));

        GEN_q_nDaughters.push_back(genPart->numberOfDaughters());
        int ndau = genPart->numberOfDaughters();
        for (int d=0; d<ndau; d++){
          const reco::Candidate *qdau=genPart->daughter(d);
          dau_id.push_back(qdau->pdgId());
          dau_pt.push_back(qdau->pt());
          dau_eta.push_back(qdau->eta());
          dau_phi.push_back(qdau->phi());
          dau_mass.push_back(qdau->mass());
          dau_status.push_back(qdau->status());
        }

        GEN_qdau_id.push_back(dau_id);
        GEN_qdau_pt.push_back(dau_pt);
        GEN_qdau_eta.push_back(dau_eta);
        GEN_qdau_phi.push_back(dau_phi);
        GEN_qdau_mass.push_back(dau_mass);
        GEN_qdau_status.push_back(dau_status);

      }


      if ((genPart->pdgId()==23 || genPart->pdgId()==443 || genPart->pdgId()==553) && (genPart->status()>=20 && genPart->status()<30) ) {
          const reco::Candidate *Zdau0=genPart->daughter(0);
          int ZdauId = fabs(Zdau0->pdgId());
          if (fabs(Zdau0->pdgId())==23) {
              int ndau = genPart->numberOfDaughters();
              for (int d=0; d<ndau; d++) {
                  const reco::Candidate *Zdau=genPart->daughter(d);
                  if (verbose) cout<<"ZDau "<<d<<" id "<<fabs(Zdau->pdgId())<<endl;
                  if (fabs(Zdau->pdgId())<17) {
                      ZdauId = fabs(Zdau->pdgId());
                      break;
                  }
              }
          }
          if (verbose) cout<<"GENZ status "<<genPart->status()<<" MomId: "<<genAna.MotherID(&prunedgenParticles->at(j))<< "DauId: "<<ZdauId<<endl;

          if (Zdau0) GENZ_DaughtersId.push_back(ZdauId);
          GENZ_MomId.push_back(genAna.MotherID(&prunedgenParticles->at(j)));
          GENZ_pt.push_back(genPart->pt());
          GENZ_eta.push_back(genPart->eta());
          GENZ_phi.push_back(genPart->phi());
          GENZ_mass.push_back(genPart->mass());
      }

      if (abs(genPart->pdgId())>500 && abs(genPart->pdgId())<600 && genPart->status()==2) {
          nGenStatus2bHad+=1;
      }

    }

  // DO GEN JETS
  if (verbose) cout<<"begin filling gen jets"<<endl;
  edm::View<reco::GenJet>::const_iterator genjet;
  for(genjet = genJets->begin(); genjet != genJets->end(); genjet++) {

    double pt = genjet->pt();
    //double eta = genjet->eta();
    //if (pt<30.0 || abs(eta)>4.7) continue;
    if (pt<30.0) continue;

    bool inDR_pt30_eta4p7 = false;
    unsigned int N=GENlep_pt.size();
    for(unsigned int i = 0; i<N; i++) {
      //if (GENlep_status[i]!=1) continue;
      if (!(abs(GENlep_id[i])==11 || abs(GENlep_id[i])==13)) continue;
      TLorentzVector genlep;
      genlep.SetPtEtaPhiM(GENlep_pt[i],GENlep_eta[i],GENlep_phi[i],GENlep_mass[i]);
      double dR = deltaR(genlep.Eta(), genlep.Phi(), genjet->eta(),genjet->phi());
      if(dR<0.4) {
        inDR_pt30_eta4p7=true;
      }
    }

    if (verbose) cout<<"check overlap of gen jet with gen leptons"<<endl;
    // count number of gen jets which no gen leptons are inside its cone
    if (!inDR_pt30_eta4p7) {
      GENnjets_pt30_eta4p7++;
      GENjet_pt.push_back(genjet->pt());
      GENjet_eta.push_back(genjet->eta());
      GENjet_phi.push_back(genjet->phi());
      GENjet_mass.push_back(genjet->mass());
      if (pt>GENpt_leadingjet_pt30_eta4p7) {
        GENpt_leadingjet_pt30_eta4p7=pt;
        GENabsrapidity_leadingjet_pt30_eta4p7=genjet->rapidity(); //take abs later
      }
      if (abs(genjet->eta())<2.5) {
        GENnjets_pt30_eta2p5++;
        if (pt>GENpt_leadingjet_pt30_eta2p5) {
          GENpt_leadingjet_pt30_eta2p5=pt;
        }
      }
    }

  }// loop over gen jets



}

void UFHZZ4LAna::bookPassedEventTree(TString treeName, TTree *tree){

  using namespace edm;
  using namespace pat;
  using namespace std;

  // -------------------------
  // RECO level information
  // -------------------------
  // Event variables
  tree->Branch("Run",&Run,"Run/l");
  tree->Branch("Event",&Event,"Event/l");
  tree->Branch("LumiSect",&LumiSect,"LumiSect/l");
  tree->Branch("nVtx",&nVtx,"nVtx/I");
  tree->Branch("nInt",&nInt,"nInt/I");
  tree->Branch("PV_x", &PV_x, "PV_x/F");
  tree->Branch("PV_y", &PV_y, "PV_y/F");
  tree->Branch("PV_z", &PV_z, "PV_z/F");
  tree->Branch("BS_x", &BS_x, "BS_x/F");
  tree->Branch("BS_y", &BS_y, "BS_y/F");
  tree->Branch("BS_z", &BS_z, "BS_z/F");
  tree->Branch("BS_xErr", &BS_xErr, "BS_xErr/F");
  tree->Branch("BS_yErr", &BS_yErr, "BS_yErr/F");
  tree->Branch("BS_zErr", &BS_zErr, "BS_zErr/F");
  tree->Branch("BeamWidth_x", &BeamWidth_x, "BeamWidth_x/F");
  tree->Branch("BeamWidth_y", &BeamWidth_y, "BeamWidth_y/F");
  tree->Branch("BeamWidth_xErr", &BeamWidth_xErr, "BeamWidth_xErr/F");
  tree->Branch("BeamWidth_yErr", &BeamWidth_yErr, "BeamWidth_yErr/F");
  tree->Branch("finalState",&finalState,"finalState/I");
  tree->Branch("triggersPassed",&triggersPassed);
  tree->Branch("passedTrig",&passedTrig,"passedTrig/O");
  tree->Branch("genWeight",&genWeight,"genWeight/F");
  tree->Branch("qcdWeights",&qcdWeights);
  tree->Branch("nnloWeights",&nnloWeights);
  tree->Branch("pdfWeights",&pdfWeights);
  tree->Branch("pdfRMSup",&pdfRMSup,"pdfRMSup/F");
  tree->Branch("pdfRMSdown",&pdfRMSdown,"pdfRMSdown/F");
  tree->Branch("pdfENVup",&pdfENVup,"pdfENVup/F");
  tree->Branch("pdfENVdown",&pdfENVdown,"pdfENVdown/F");
  tree->Branch("pileupWeight",&pileupWeight,"pileupWeight/F");
  tree->Branch("pileupWeightUp",&pileupWeightUp,"pileupWeightUp/F");
  tree->Branch("pileupWeightDn",&pileupWeightDn,"pileupWeightDn/F");
  tree->Branch("dataMCWeight",&dataMCWeight,"dataMCWeight/F");
  tree->Branch("eventWeight",&eventWeight,"eventWeight/F");
  tree->Branch("prefiringWeight",&prefiringWeight,"prefiringWeight/F");
  tree->Branch("prefiringWeightECAL",&prefiringWeightECAL,"prefiringWeightECAL/F");
  tree->Branch("prefiringWeightMuon",&prefiringWeightMuon,"prefiringWeightMuon/F");
  tree->Branch("crossSection",&crossSection,"crossSection/F");

  // Lepton variables
  tree->Branch("lep_d0BS",&lep_d0BS);
  tree->Branch("lep_d0PV",&lep_d0PV);
  tree->Branch("lep_numberOfValidPixelHits",&lep_numberOfValidPixelHits);
  tree->Branch("lep_trackerLayersWithMeasurement",&lep_trackerLayersWithMeasurement);

  tree->Branch("lep_id",&lep_id);
  tree->Branch("lep_pt",&lep_pt);
  tree->Branch("lep_pterr",&lep_pterr);
  tree->Branch("lep_pterrold",&lep_pterrold);
  tree->Branch("lep_eta",&lep_eta);
  tree->Branch("lep_phi",&lep_phi);
  tree->Branch("lep_mass",&lep_mass);
  tree->Branch("lepFSR_pt",&lepFSR_pt);
  tree->Branch("lepFSR_eta",&lepFSR_eta);
  tree->Branch("lepFSR_phi",&lepFSR_phi);
  tree->Branch("lepFSR_mass",&lepFSR_mass);
  //tree->Branch("lep_Hindex",&lep_Hindex,"lep_Hindex[4]/I");
  tree->Branch("lep_genindex",&lep_genindex);
  tree->Branch("lep_matchedR03_PdgId",&lep_matchedR03_PdgId);
  tree->Branch("lep_matchedR03_MomId",&lep_matchedR03_MomId);
  tree->Branch("lep_matchedR03_MomMomId",&lep_matchedR03_MomMomId);
  tree->Branch("lep_missingHits",&lep_missingHits);
  tree->Branch("lep_mva",&lep_mva);
  tree->Branch("lep_ecalDriven",&lep_ecalDriven);
  tree->Branch("lep_tightId",&lep_tightId);
  //tree->Branch("lep_tightId_old",&lep_tightId_old);
  tree->Branch("lep_tightIdSUS",&lep_tightIdSUS);
  tree->Branch("lep_tightIdHiPt",&lep_tightIdHiPt);
  tree->Branch("lep_Sip",&lep_Sip);
  tree->Branch("lep_IP",&lep_IP);
  tree->Branch("lep_isoNH",&lep_isoNH);
  tree->Branch("lep_isoCH",&lep_isoCH);
  tree->Branch("lep_isoPhot",&lep_isoPhot);
  tree->Branch("lep_isoPU",&lep_isoPU);
  tree->Branch("lep_isoPUcorr",&lep_isoPUcorr);
  tree->Branch("lep_RelIso",&lep_RelIso);
  tree->Branch("lep_RelIsoNoFSR",&lep_RelIsoNoFSR);
  tree->Branch("lep_MiniIso",&lep_MiniIso);
  tree->Branch("lep_ptRatio",&lep_ptRatio);
  tree->Branch("lep_ptRel",&lep_ptRel);
  tree->Branch("lep_filtersMatched",&lep_filtersMatched);
  tree->Branch("lep_dataMC",&lep_dataMC);
  tree->Branch("lep_dataMCErr",&lep_dataMCErr);
  tree->Branch("dataMC_VxBS",&dataMC_VxBS);
  tree->Branch("dataMCErr_VxBS",&dataMCErr_VxBS);
  tree->Branch("nisoleptons",&nisoleptons,"nisoleptons/I");

  tree->Branch("singleBS_RecoLep_pt",&singleBS_RecoLep_pt);
  tree->Branch("singleBS_RecoLep_ptError",&singleBS_RecoLep_ptError);
  tree->Branch("singleBS_RecoLep_eta",&singleBS_RecoLep_eta);
  tree->Branch("singleBS_RecoLep_phi",&singleBS_RecoLep_phi);
  tree->Branch("singleBS_RecoLep_mass",&singleBS_RecoLep_mass);
  tree->Branch("singleBS_RecoLep_d0",&singleBS_RecoLep_d0);

  tree->Branch("pho_pt",&pho_pt);
  tree->Branch("pho_eta",&pho_eta);
  tree->Branch("pho_phi",&pho_phi);
  tree->Branch("photonCutBasedIDLoose",&photonCutBasedIDLoose);

  // MET
  tree->Branch("met",&met,"met/F");
  tree->Branch("met_phi",&met_phi,"met_phi/F");
  tree->Branch("met_jesup",&met_jesup,"met_jesup/F");
  tree->Branch("met_phi_jesup",&met_phi_jesup,"met_phi_jesup/F");
  tree->Branch("met_jesdn",&met_jesdn,"met_jesdn/F");
  tree->Branch("met_phi_jesdn",&met_phi_jesdn,"met_phi_jesdn/F");
  tree->Branch("met_uncenup",&met_uncenup,"met_uncenup/F");
  tree->Branch("met_phi_uncenup",&met_phi_uncenup,"met_phi_uncenup/F");
  tree->Branch("met_uncendn",&met_uncendn,"met_uncendn/F");
  tree->Branch("met_phi_uncendn",&met_phi_uncendn,"met_phi_uncendn/F");

  // Jets
  tree->Branch("jet_iscleanH4l",&jet_iscleanH4l);
  tree->Branch("jet1index",&jet1index,"jet1index/I");
  tree->Branch("jet2index",&jet2index,"jet2index/I");
  tree->Branch("jet_pt",&jet_pt);
  tree->Branch("jet_pt_raw",&jet_pt_raw);
  tree->Branch("jet_relpterr",&jet_relpterr);
  tree->Branch("jet_eta",&jet_eta);
  tree->Branch("jet_phi",&jet_phi);
  tree->Branch("jet_phierr",&jet_phierr);
  tree->Branch("jet_bTagEffi",&jet_bTagEffi);
  tree->Branch("jet_cTagEffi",&jet_cTagEffi);
  tree->Branch("jet_udsgTagEffi",&jet_udsgTagEffi);
  tree->Branch("jet_mass",&jet_mass);
  tree->Branch("jet_jesup_iscleanH4l",&jet_jesup_iscleanH4l);
  tree->Branch("jet_jesup_pt",&jet_jesup_pt);
  tree->Branch("jet_jesup_eta",&jet_jesup_eta);
  tree->Branch("jet_jesup_phi",&jet_jesup_phi);
  tree->Branch("jet_jesup_mass",&jet_jesup_mass);
  tree->Branch("jet_jesdn_iscleanH4l",&jet_jesdn_iscleanH4l);
  tree->Branch("jet_jesdn_pt",&jet_jesdn_pt);
  tree->Branch("jet_jesdn_eta",&jet_jesdn_eta);
  tree->Branch("jet_jesdn_phi",&jet_jesdn_phi);
  tree->Branch("jet_jesdn_mass",&jet_jesdn_mass);
  tree->Branch("jet_jerup_iscleanH4l",&jet_jerup_iscleanH4l);
  tree->Branch("jet_jerup_pt",&jet_jerup_pt);
  tree->Branch("jet_jerup_eta",&jet_jerup_eta);
  tree->Branch("jet_jerup_phi",&jet_jerup_phi);
  tree->Branch("jet_jerup_mass",&jet_jerup_mass);
  tree->Branch("jet_jerdn_iscleanH4l",&jet_jerdn_iscleanH4l);
  tree->Branch("jet_jerdn_pt",&jet_jerdn_pt);
  tree->Branch("jet_jerdn_eta",&jet_jerdn_eta);
  tree->Branch("jet_jerdn_phi",&jet_jerdn_phi);
  tree->Branch("jet_jerdn_mass",&jet_jerdn_mass);
  tree->Branch("jet_pumva",&jet_pumva);
  tree->Branch("jet_csvv2",&jet_csvv2);
  tree->Branch("jet_csvv2_",&jet_csvv2_);
  tree->Branch("jet_isbtag",&jet_isbtag);
  tree->Branch("jet_hadronFlavour",&jet_hadronFlavour);
  tree->Branch("jet_partonFlavour",&jet_partonFlavour);
  tree->Branch("jet_QGTagger",&jet_QGTagger);
  tree->Branch("jet_QGTagger_jesup",&jet_QGTagger_jesup);
  tree->Branch("jet_QGTagger_jesdn",&jet_QGTagger_jesdn);
  tree->Branch("jet_axis2",&jet_axis2);
  tree->Branch("jet_ptD",&jet_ptD);
  tree->Branch("jet_mult",&jet_mult);
  tree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7,"njets_pt30_eta4p7/I");
  tree->Branch("njets_pt30_eta4p7_jesup",&njets_pt30_eta4p7_jesup,"njets_pt30_eta4p7_jesup/I");
  tree->Branch("njets_pt30_eta4p7_jesdn",&njets_pt30_eta4p7_jesdn,"njets_pt30_eta4p7_jesdn/I");
  tree->Branch("njets_pt30_eta4p7_jerup",&njets_pt30_eta4p7_jerup,"njets_pt30_eta4p7_jerup/I");
  tree->Branch("njets_pt30_eta4p7_jerdn",&njets_pt30_eta4p7_jerdn,"njets_pt30_eta4p7_jerdn/I");
  tree->Branch("pt_leadingjet_pt30_eta4p7",&pt_leadingjet_pt30_eta4p7,"pt_leadingjet_pt30_eta4p7/F");
  tree->Branch("pt_leadingjet_pt30_eta4p7_jesup",&pt_leadingjet_pt30_eta4p7_jesup,"pt_leadingjet_pt30_eta4p7_jesup/F");
  tree->Branch("pt_leadingjet_pt30_eta4p7_jesdn",&pt_leadingjet_pt30_eta4p7_jesdn,"pt_leadingjet_pt30_eta4p7_jesdn/F");
  tree->Branch("pt_leadingjet_pt30_eta4p7_jerup",&pt_leadingjet_pt30_eta4p7_jerup,"pt_leadingjet_pt30_eta4p7_jerup/F");
  tree->Branch("pt_leadingjet_pt30_eta4p7_jerdn",&pt_leadingjet_pt30_eta4p7_jerdn,"pt_leadingjet_pt30_eta4p7_jerdn/F");
  tree->Branch("absrapidity_leadingjet_pt30_eta4p7",&absrapidity_leadingjet_pt30_eta4p7,"absrapidity_leadingjet_pt30_eta4p7/F");
  tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jesup",&absrapidity_leadingjet_pt30_eta4p7_jesup,"absrapidity_leadingjet_pt30_eta4p7_jesup/F");
  tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jesdn",&absrapidity_leadingjet_pt30_eta4p7_jesdn,"absrapidity_leadingjet_pt30_eta4p7_jesdn/F");
  tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jerup",&absrapidity_leadingjet_pt30_eta4p7_jerup,"absrapidity_leadingjet_pt30_eta4p7_jerup/F");
  tree->Branch("absrapidity_leadingjet_pt30_eta4p7_jerdn",&absrapidity_leadingjet_pt30_eta4p7_jerdn,"absrapidity_leadingjet_pt30_eta4p7_jerdn/F");
  tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7",&absdeltarapidity_hleadingjet_pt30_eta4p7,"absdeltarapidity_hleadingjet_pt30_eta4p7/F");
  tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jesup",&absdeltarapidity_hleadingjet_pt30_eta4p7_jesup,"absdeltarapidity_hleadingjet_pt30_eta4p7_jesup/F");
  tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn",&absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn,"absdeltarapidity_hleadingjet_pt30_eta4p7_jesdn/F");
  tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jerup",&absdeltarapidity_hleadingjet_pt30_eta4p7_jerup,"absdeltarapidity_hleadingjet_pt30_eta4p7_jerup/F");
  tree->Branch("absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn",&absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn,"absdeltarapidity_hleadingjet_pt30_eta4p7_jerdn/F");
  tree->Branch("nbjets_pt30_eta4p7",&nbjets_pt30_eta4p7,"nbjets_pt30_eta4p7/I");
  tree->Branch("nvjets_pt40_eta2p4",&nvjets_pt40_eta2p4,"nvjets_pt40_eta2p4/I");
  tree->Branch("DijetMass",&DijetMass,"DijetMass/F");
  tree->Branch("DijetDEta",&DijetDEta,"DijetDEta/F");
  tree->Branch("DijetFisher",&DijetFisher,"DijetFisher/F");
  tree->Branch("njets_pt30_eta2p5",&njets_pt30_eta2p5,"njets_pt30_eta2p5/I");
  tree->Branch("njets_pt30_eta2p5_jesup",&njets_pt30_eta2p5_jesup,"njets_pt30_eta2p5_jesup/I");
  tree->Branch("njets_pt30_eta2p5_jesdn",&njets_pt30_eta2p5_jesdn,"njets_pt30_eta2p5_jesdn/I");
  tree->Branch("njets_pt30_eta2p5_jerup",&njets_pt30_eta2p5_jerup,"njets_pt30_eta2p5_jerup/I");
  tree->Branch("njets_pt30_eta2p5_jerdn",&njets_pt30_eta2p5_jerdn,"njets_pt30_eta2p5_jerdn/I");
  tree->Branch("pt_leadingjet_pt30_eta2p5",&pt_leadingjet_pt30_eta2p5,"pt_leadingjet_pt30_eta2p5/F");
  tree->Branch("pt_leadingjet_pt30_eta2p5_jesup",&pt_leadingjet_pt30_eta2p5_jesup,"pt_leadingjet_pt30_eta2p5_jesup/F");
  tree->Branch("pt_leadingjet_pt30_eta2p5_jesdn",&pt_leadingjet_pt30_eta2p5_jesdn,"pt_leadingjet_pt30_eta2p5_jesdn/F");
  tree->Branch("pt_leadingjet_pt30_eta2p5_jerup",&pt_leadingjet_pt30_eta2p5_jerup,"pt_leadingjet_pt30_eta2p5_jerup/F");
  tree->Branch("pt_leadingjet_pt30_eta2p5_jerdn",&pt_leadingjet_pt30_eta2p5_jerdn,"pt_leadingjet_pt30_eta2p5_jerdn/F");

  // merged jets
  tree->Branch("mergedjet_iscleanH4l",&mergedjet_iscleanH4l);
  tree->Branch("mergedjet_pt",&mergedjet_pt);
  tree->Branch("mergedjet_ptrow",&mergedjet_ptrow);
  tree->Branch("mergedjet_jer",&mergedjet_jer);
  tree->Branch("mergedjet_eta",&mergedjet_eta);
  tree->Branch("mergedjet_phi",&mergedjet_phi);
  tree->Branch("mergedjet_mass",&mergedjet_mass);

  tree->Branch("mergedjet_jer_pt",&mergedjet_jer_pt);
  tree->Branch("mergedjet_jer_eta",&mergedjet_jer_eta);
  tree->Branch("mergedjet_jer_phi",&mergedjet_jer_phi);
  tree->Branch("mergedjet_jer_mass",&mergedjet_jer_mass);
  tree->Branch("mergedjet_jerup_pt",&mergedjet_jerup_pt);
  tree->Branch("mergedjet_jerup_eta",&mergedjet_jerup_eta);
  tree->Branch("mergedjet_jerup_phi",&mergedjet_jerup_phi);
  tree->Branch("mergedjet_jerup_mass",&mergedjet_jerup_mass);
  tree->Branch("mergedjet_jerdn_pt",&mergedjet_jerdn_pt);
  tree->Branch("mergedjet_jerdn_eta",&mergedjet_jerdn_eta);
  tree->Branch("mergedjet_jerdn_phi",&mergedjet_jerdn_phi);
  tree->Branch("mergedjet_jerdn_mass",&mergedjet_jerdn_mass);
  tree->Branch("mergedjet_jer_pterr",&mergedjet_jer_pterr);
  tree->Branch("mergedjet_jer_phierr",&mergedjet_jer_phierr);


  tree->Branch("mergedjet_tau1",&mergedjet_tau1);
  tree->Branch("mergedjet_tau2",&mergedjet_tau2);
  tree->Branch("mergedjet_btag",&mergedjet_btag);

  tree->Branch("mergedjet_ZvsQCD",&mergedjet_ZvsQCD);
  tree->Branch("mergedjet_ZbbvsQCD",&mergedjet_ZbbvsQCD);
  tree->Branch("mergedjet_WvsQCD",&mergedjet_WvsQCD);
  tree->Branch("mergedjet_ZHbbvsQCD",&mergedjet_ZHbbvsQCD);
  tree->Branch("mergedjet_HbbvsQCD",&mergedjet_H4qvsQCD);
  tree->Branch("mergedjet_ZvsQCD_de",&mergedjet_ZvsQCD_de);
  tree->Branch("mergedjet_ZbbvsQCD_de",&mergedjet_ZbbvsQCD_de);
  tree->Branch("mergedjet_WvsQCD_de",&mergedjet_WvsQCD_de);
  tree->Branch("mergedjet_ZHbbvsQCD_de",&mergedjet_ZHbbvsQCD_de);
  tree->Branch("mergedjet_HbbvsQCD_de",&mergedjet_H4qvsQCD_de);

  tree->Branch("mergedjet_Net_Xbb_de",&mergedjet_Net_Xbb_de);
  tree->Branch("mergedjet_Net_Xcc_de",&mergedjet_Net_Xcc_de);
  tree->Branch("mergedjet_Net_Xqq_de",&mergedjet_Net_Xqq_de);
  tree->Branch("mergedjet_Net_QCDbb_de",&mergedjet_Net_QCDbb_de);
  tree->Branch("mergedjet_Net_QCDcc_de",&mergedjet_Net_QCDcc_de);
  tree->Branch("mergedjet_Net_QCDother_de",&mergedjet_Net_QCDother_de);
  tree->Branch("mergedjet_Net_QCDb_de",&mergedjet_Net_QCDb_de);
  tree->Branch("mergedjet_Net_QCDc_de",&mergedjet_Net_QCDc_de);

  //tree->Branch("mergedjet_L1",&mergedjet_L1);
  tree->Branch("mergedjet_softdropmass",&mergedjet_softdropmass);
  //tree->Branch("mergedjet_prunedmass",&mergedjet_prunedmass);

  tree->Branch("mergedjet_nsubjet",&mergedjet_nsubjet);
  tree->Branch("mergedjet_subjet_pt",&mergedjet_subjet_pt);
  tree->Branch("mergedjet_subjet_eta",&mergedjet_subjet_eta);
  tree->Branch("mergedjet_subjet_phi",&mergedjet_subjet_phi);
  tree->Branch("mergedjet_subjet_mass",&mergedjet_subjet_mass);
  tree->Branch("mergedjet_subjet_softDropMass",&mergedjet_subjet_softDropMass);
  tree->Branch("mergedjet_subjet_btag",&mergedjet_subjet_btag);
  tree->Branch("mergedjet_subjetUncorr_pt",&mergedjet_subjetUncorr_pt);
  tree->Branch("mergedjet_subjetUncorr_eta",&mergedjet_subjetUncorr_eta);
  tree->Branch("mergedjet_subjetUncorr_phi",&mergedjet_subjetUncorr_phi);
  tree->Branch("mergedjet_subjetUncorr_mass",&mergedjet_subjetUncorr_mass);
  tree->Branch("mergedjet_subjet_softDropMassUncorr",&mergedjet_subjet_softDropMassUncorr);
  tree->Branch("mergedjet_nbHadrons",&mergedjet_nbHadrons);
  tree->Branch("mergedjet_ncHadrons",&mergedjet_ncHadrons);
  tree->Branch("mergedjet_subjet_partonFlavour",&mergedjet_subjet_partonFlavour);
  tree->Branch("mergedjet_subjet_hadronFlavour",&mergedjet_subjet_hadronFlavour);


  // FSR Photons
  tree->Branch("nFSRPhotons",&nFSRPhotons,"nFSRPhotons/I");
  tree->Branch("allfsrPhotons_dR",&allfsrPhotons_dR);
  tree->Branch("allfsrPhotons_iso",&allfsrPhotons_iso);
  tree->Branch("allfsrPhotons_pt",&allfsrPhotons_pt);
  tree->Branch("fsrPhotons_lepindex",&fsrPhotons_lepindex);
  tree->Branch("fsrPhotons_pt",&fsrPhotons_pt);
  tree->Branch("fsrPhotons_pterr",&fsrPhotons_pterr);
  tree->Branch("fsrPhotons_eta",&fsrPhotons_eta);
  tree->Branch("fsrPhotons_phi",&fsrPhotons_phi);
  tree->Branch("fsrPhotons_dR",&fsrPhotons_dR);
  tree->Branch("fsrPhotons_iso",&fsrPhotons_iso);

  // Event Category
  tree->Branch("EventCat",&EventCat,"EventCat/I");

  // -------------------------
  // GEN level information
  // -------------------------
  //Event variables
  // lepton variables
  tree->Branch("GENlep_pt",&GENlep_pt);
  tree->Branch("GENlep_eta",&GENlep_eta);
  tree->Branch("GENlep_phi",&GENlep_phi);
  tree->Branch("GENlep_mass",&GENlep_mass);
  tree->Branch("GENlep_id",&GENlep_id);
  tree->Branch("GENlep_status",&GENlep_status);
  tree->Branch("GENlep_MomId",&GENlep_MomId);
  tree->Branch("GENlep_MomMomId",&GENlep_MomMomId);
  tree->Branch("GENlep_Hindex",&GENlep_Hindex,"GENlep_Hindex[4]/I");
  tree->Branch("GENlep_isoCH",&GENlep_isoCH);
  tree->Branch("GENlep_isoNH",&GENlep_isoNH);
  tree->Branch("GENlep_isoPhot",&GENlep_isoPhot);
  tree->Branch("GENlep_RelIso",&GENlep_RelIso);

  // Higgs candidate variables (calculated using selected gen leptons)
  tree->Branch("GENH_pt",&GENH_pt);
  tree->Branch("GENH_eta",&GENH_eta);
  tree->Branch("GENH_phi",&GENH_phi);
  tree->Branch("GENH_mass",&GENH_mass);
  tree->Branch("GENMH",&GENMH,"GENMH/F");
  tree->Branch("GENH_Momid",&GENH_Momid);
  tree->Branch("GENH_MomMomid",&GENH_MomMomid);
  tree->Branch("GENH_status",&GENH_status);
  tree->Branch("GENH_id",&GENH_id);
  tree->Branch("GENH_isHard",&GENH_isHard);
  tree->Branch("GENH_nDaughters",&GENH_nDaughters);
  tree->Branch("GEN_id",&GEN_id);
  tree->Branch("GEN_status",&GEN_status);
  tree->Branch("GENH_dau_id",&GENH_dau_id);
  tree->Branch("GENH_dau_status",&GENH_dau_status);
  tree->Branch("GENH_dau_pt",&GENH_dau_pt);
  tree->Branch("GENH_dau_eta",&GENH_dau_eta);
  tree->Branch("GENH_dau_phi",&GENH_dau_phi);
  tree->Branch("GENH_dau_mass",&GENH_dau_mass);


  // Z candidate variables
  tree->Branch("GENZ_pt",&GENZ_pt);
  tree->Branch("GENZ_eta",&GENZ_eta);
  tree->Branch("GENZ_phi",&GENZ_phi);
  tree->Branch("GENZ_mass",&GENZ_mass);
  tree->Branch("GENZ_DaughtersId",&GENZ_DaughtersId);
  tree->Branch("GENZ_MomId",&GENZ_MomId);

  //hadZ and quark
  tree->Branch("GEN_Zq_pt",&GEN_Zq_pt);
  tree->Branch("GEN_Zq_eta",&GEN_Zq_eta);
  tree->Branch("GEN_Zq_phi",&GEN_Zq_phi);
  tree->Branch("GEN_Zq_mass",&GEN_Zq_mass);
  tree->Branch("GEN_Zq_id",&GEN_Zq_id);
  tree->Branch("GEN_Zq_Momid",&GEN_Zq_Momid);
  tree->Branch("GEN_Zq_MomMomid",&GEN_Zq_MomMomid);
  tree->Branch("GEN_Zq_isHard",&GEN_Zq_isHard);

  tree->Branch("GEN_q_id",&GEN_q_id);
  tree->Branch("GEN_q_pt",&GEN_q_pt);
  tree->Branch("GEN_q_eta",&GEN_q_eta);
  tree->Branch("GEN_q_phi",&GEN_q_phi);
  tree->Branch("GEN_q_mass",&GEN_q_mass);
  tree->Branch("GEN_q_status",&GEN_q_status);
  tree->Branch("GEN_q_Momid",&GEN_q_Momid);
  tree->Branch("GEN_q_MomMomid",&GEN_q_MomMomid);
  tree->Branch("GEN_q_nDaughters",&GEN_q_nDaughters);

  tree->Branch("GEN_qdau_id",&GEN_qdau_id);
  tree->Branch("GEN_qdau_pt",&GEN_qdau_pt);
  tree->Branch("GEN_qdau_eta",&GEN_qdau_eta);
  tree->Branch("GEN_qdau_phi",&GEN_qdau_phi);
  tree->Branch("GEN_qdau_mass",&GEN_qdau_mass);
  tree->Branch("GEN_qdau_status",&GEN_qdau_status);

  tree->Branch("GEN_VBF_pt",&GEN_VBF_pt);
  tree->Branch("GEN_VBF_eta",&GEN_VBF_eta);
  tree->Branch("GEN_VBF_phi",&GEN_VBF_phi);
  tree->Branch("GEN_VBF_mass",&GEN_VBF_mass);
  tree->Branch("GEN_VBF_id",&GEN_VBF_id);
  tree->Branch("GEN_VBF_Momid",&GEN_VBF_Momid);
  tree->Branch("GEN_VBF_MomMomid",&GEN_VBF_MomMomid);
  tree->Branch("GEN_VBF_status",&GEN_VBF_status);




  // Jets
  tree->Branch("GENjet_pt",&GENjet_pt);
  tree->Branch("GENjet_eta",&GENjet_eta);
  tree->Branch("GENjet_phi",&GENjet_phi);
  //tree->Branch("GENjet_id",&GENjet_id);
  tree->Branch("GENjet_id",&GENjet_id);
  tree->Branch("GENjet_mass",&GENjet_mass);
  tree->Branch("GENnjets_pt30_eta4p7",&GENnjets_pt30_eta4p7,"GENnjets_pt30_eta4p7/I");
  tree->Branch("GENpt_leadingjet_pt30_eta4p7",&GENpt_leadingjet_pt30_eta4p7,"GENpt_leadingjet_pt30_eta4p7/F");
  tree->Branch("GENnbjets_pt30_eta4p7",&GENnbjets_pt30_eta4p7,"GENnbjets_pt30_eta4p7/I");
  tree->Branch("GENabsrapidity_leadingjet_pt30_eta4p7",&GENabsrapidity_leadingjet_pt30_eta4p7,"GENabsrapidity_leadingjet_pt30_eta4p7/F");
  tree->Branch("GENabsdeltarapidity_hleadingjet_pt30_eta4p7",&GENabsdeltarapidity_hleadingjet_pt30_eta4p7,"GENabsdeltarapidity_hleadingjet_pt30_eta4p7/F");
  tree->Branch("GENnjets_pt30_eta2p5",&GENnjets_pt30_eta2p5,"GENnjets_pt30_eta2p5/I");
  tree->Branch("GENpt_leadingjet_pt30_eta2p5",&GENpt_leadingjet_pt30_eta2p5,"GENpt_leadingjet_pt30_eta2p5/F");
  tree->Branch("lheNj",&lheNj,"lheNj/I");
  tree->Branch("lheNb",&lheNb,"lheNb/I");
  tree->Branch("nGenStatus2bHad",&nGenStatus2bHad,"nGenStatus2bHad/I");
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
UFHZZ4LAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(UFHZZ4LAna);

//  LocalWords:  ecalDriven
