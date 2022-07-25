#ifndef UFHZZ4LAna_H
#define UFHZZ4LAna_H


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
//#include "IvyBase.h"
//#include <IvyFramework/IvyAutoMELA/interface/IvyMELAHelpers.h>
//#include <IvyFramework/IvyAutoMELA/interface/IvyMELAOutputStreamerExt.h>
//#include <IvyFramework/IvyDataTools/interface/HostHelpersCore.h>


//
// class declaration
//
//using namespace IvyStreamHelpers;
using namespace EwkCorrections;

class UFHZZ4LAna : public edm::EDAnalyzer{
public:
  explicit UFHZZ4LAna(const edm::ParameterSet&);
  ~UFHZZ4LAna();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static bool sortByPt( const reco::GenParticle &p1, const reco::GenParticle &p2 ){ return (p1.pt() > p2.pt()); };
  void Reset();
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
  int nlooseleps, ntightleps, nGenV; //##am
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

  //#######am 

  vector<double>    GenV_pt; 
  vector<int>       GenV_status; 
  vector<double>    GenV_eta; 
  vector<double>    GenV_phi; 
  vector<double>    GenV_mass;  
  vector<int>       GenV_pdgId; 
  vector<int>       GenV_ndau;
  std::vector<bool> GenV_hadronic;
  vector<int>       GenVdau_pdgId;
  vector<int>       GenVdau_MompdgId;
  vector<double>    GenVdau_pt;
  vector<double>    GenVdau_eta;
  vector<double>    GenVdau_phi;
  vector<double>    GenVdau_mass;
  vector<int>       GenVdau_status;
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

  std::vector<float> lhepart_pt;
  std::vector<float> lhepart_eta;
  std::vector<float> lhepart_phi;
  std::vector<float> lhepart_mass;
  std::vector<int>   lhepart_pdgId;
  std::vector<int>   lhepart_status;


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
  //vector<string> lheMElist;
  //vector<string> recoMElist;
  //IvyMELAHelpers::GMECBlock MEblock;

};

void UFHZZ4LAna::Reset(){


  nVtx = -1.0; nInt = -1.0;
  nlooseleps=-9999; ntightleps=-9999,nGenV=-999;//##am
  finalState = -1;
  triggersPassed="";
  passedTrig=false;

  // Event Weights
  genWeight=1.0; pileupWeight=1.0; pileupWeightUp=1.0; pileupWeightDn=1.0; dataMCWeight=1.0; eventWeight=1.0;

  qcdWeights.clear(); nnloWeights.clear(); pdfWeights.clear();
  pdfRMSup=1.0; pdfRMSdown=1.0; pdfENVup=1.0; pdfENVdown=1.0;
  //######am 
  lhepart_pt.clear();
  lhepart_eta.clear();
  lhepart_phi.clear();
  lhepart_mass.clear();
  lhepart_pdgId.clear();
  lhepart_status.clear();

  GenV_pt.clear(); 
  GenV_status.clear(); 
  GenV_eta.clear(); 
  GenV_phi.clear(); 
  GenV_mass.clear();  
  GenV_pdgId.clear(); 
  GenV_ndau.clear();
  GenV_hadronic.clear();
  GenVdau_pdgId.clear();
  GenVdau_MompdgId.clear();
  GenVdau_pt.clear();
  GenVdau_eta.clear();
  GenVdau_phi.clear();
  GenVdau_mass.clear();
  GenVdau_status.clear();
 
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
  GENZ_pt.clear(); GENZ_eta.clear(); GENZ_phi.clear(); GENZ_mass.clear();  GENZ_DaughtersId.clear(); GENZ_MomId.clear();


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
}//reset

#endif
