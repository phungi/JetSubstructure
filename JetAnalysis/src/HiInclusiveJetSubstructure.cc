/*
  Based on the jet response analyzer
  Modified by Matt Nguyen, November 2010
  Modified by Leticia Cunqueiro, September 2021
*/
#include "HeavyIonsAnalysis/JetAnalysis/interface/HiInclusiveJetSubstructure.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
// #include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "TLorentzVector.h"
#include <random>
#include "TRandom.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include <math.h>

using namespace std;
using namespace edm;
using namespace reco;

// float delta_phi(float phi1, float phi2){
//     float result = phi1 - phi2;
//     while ( result < -M_PI )
//       {
//         result += 2.*M_PI;
//       }
//     while ( result > M_PI )
//       {
//         result -= 2.*M_PI;
//       }
//     return result;
// }

HiInclusiveJetSubstructure::HiInclusiveJetSubstructure(const edm::ParameterSet& iConfig)
{
  doMatch_ = iConfig.getUntrackedParameter<bool>("matchJets",false);
  jetTag_ = consumes<pat::JetCollection> (iConfig.getParameter<InputTag>("jetTag"));
  matchTag_ = consumes<pat::JetCollection> (iConfig.getUntrackedParameter<InputTag>("matchTag"));

  vtxTag_ = consumes<vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vtxTag"));
  //vtxTag_ = consumes<vector<reco::Vertex>>(iConfig.getUntrackedParameter<edm::InputTag>("vtxTag",edm::InputTag("offlinePrimaryVertices")));
  trackTag_ = consumes<reco::TrackCollection> (iConfig.getParameter<InputTag>("trackTag"));
  useQuality_ = iConfig.getUntrackedParameter<bool>("useQuality",1);
  trackQuality_ = iConfig.getUntrackedParameter<string>("trackQuality","highPurity");

  jetName_ = iConfig.getUntrackedParameter<string>("jetName");
  doGenTaus_ = iConfig.getUntrackedParameter<bool>("doGenTaus",0);
  doGenSym_ = iConfig.getUntrackedParameter<bool>("doGenSym",0);
  doSubJets_ = iConfig.getUntrackedParameter<bool>("doSubJets",0);
  doJetConstituents_ = iConfig.getUntrackedParameter<bool>("doJetConstituents", false);
  doGenSubJets_ = iConfig.getUntrackedParameter<bool>("doGenSubJets", false);
  if (doGenSubJets_)
    subjetGenTag_ = consumes<reco::JetView> (iConfig.getUntrackedParameter<InputTag>("subjetGenTag"));
  // subjetGenTag_ = consumes<reco::JetView> (iConfig.getUntrackedParameter<InputTag>("subjetGenTag"));

  //reWTA reclustering
  doWTARecluster_ = iConfig.getUntrackedParameter<bool>("doWTARecluster", false);
/*
  if(doGenTaus_){
    tokenGenTau1_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("genTau1"));
    tokenGenTau2_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("genTau2"));
    tokenGenTau3_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("genTau3"));
  }
*/
  if (doGenSym_){
    tokenGenSym_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("genSym"));
    tokenGenDroppedBranches_ = consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("genDroppedBranches"));
  }

  isMC_ = iConfig.getUntrackedParameter<bool>("isMC",false);
  useHepMC_ = iConfig.getUntrackedParameter<bool> ("useHepMC",false);
  fillGenJets_ = iConfig.getUntrackedParameter<bool>("fillGenJets",false);

  doHiJetID_ = iConfig.getUntrackedParameter<bool>("doHiJetID",false);
  doStandardJetID_ = iConfig.getUntrackedParameter<bool>("doStandardJetID",false);

  rParam = iConfig.getParameter<double>("rParam");
  hardPtMin_ = iConfig.getUntrackedParameter<double>("hardPtMin",4);
  jetPtMin_ = iConfig.getParameter<double>("jetPtMin");
  mysdcut1 = iConfig.getParameter<double>("mysdcut1");
  mysdcut2 = iConfig.getParameter<double>("mysdcut2");
  mydynktcut = iConfig.getParameter<double>("mydynktcut");
  groom_type = iConfig.getParameter<double>("groom_type");
  groom_combine = iConfig.getParameter<double>("groom_combine");
  jetAbsEtaMax_ = iConfig.getUntrackedParameter<double>("jetAbsEtaMax", 2.5);

  if(isMC_){
    genjetTag_ = consumes<edm::View<reco::GenJet>>(iConfig.getParameter<InputTag>("genjetTag"));
    if(useHepMC_) eventInfoTag_ = consumes<HepMCProduct> (iConfig.getParameter<InputTag>("eventInfoTag"));
    eventGenInfoTag_ = consumes<GenEventInfoProduct> (iConfig.getParameter<InputTag>("eventInfoTag"));
  }
  verbose_ = iConfig.getUntrackedParameter<bool>("verbose",false);
  useVtx_ = iConfig.getUntrackedParameter<bool>("useVtx",false);
  useRawPt_ = iConfig.getUntrackedParameter<bool>("useRawPt",true);

  doLifeTimeTagging_ = iConfig.getUntrackedParameter<bool>("doLifeTimeTagging",false);
  doLifeTimeCandidateTagging_ = iConfig.getUntrackedParameter<bool>("doLifeTimeCandidateTagging",false);
  doLifeTimeTaggingExtras_ = iConfig.getUntrackedParameter<bool>("doLifeTimeTaggingExtras",true);
  saveBfragments_  = iConfig.getUntrackedParameter<bool>("saveBfragments",false);

  // pfCandidateLabel_ = consumes<reco::PFCandidateCollection>(iConfig.getUntrackedParameter<edm::InputTag>("pfCandidateLabel")); 
  // pfCandidateToken_ = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandSource"));
  doTower = iConfig.getUntrackedParameter<bool>("doTower",false);
  if(doTower){
    TowerSrc_ = consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("towersSrc"));
  }

  doExtraCTagging_ = iConfig.getUntrackedParameter<bool>("doExtraCTagging",false);

  pfCandidateToken_ = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandSource"));


  if(isMC_){
    genParticleSrc_ = consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles"));
  }

  doSubEvent_ = 0;
  doChargedConstOnly_ = iConfig.getUntrackedParameter<bool>("doChargedConstOnly",0);
  TrackVariation_ = -1;
  pfChargedCandidateEnergyScale_ = -1;
  pfNeutralCandidateEnergyScale_ = -1;
  pfGammaCandidateEnergyScale_ = -1;
  pfNeutralSmear_ = false;
  doNaiveNeuPFScaling_ = false;
  doRatioNeuPFScaling_ = false;
  doPeripheralNeuPFScaling_ = false;
  doCompensatoryNeuPFScaling_ = false;

  doPrimaryLJPReco_ = iConfig.getUntrackedParameter<bool>("doPrimaryLJPReco", false);
  if(isMC_){
    doPrimaryLJPTruth_ = iConfig.getUntrackedParameter<bool>("doPrimaryLJPTruth", false);
    pfChargedCandidateEnergyScale_ = iConfig.getUntrackedParameter<double>("pfChargedEnergyScaleVar",1.);
    pfNeutralCandidateEnergyScale_ = iConfig.getUntrackedParameter<double>("pfNeutralEnergyScaleVar",1.);
    pfGammaCandidateEnergyScale_ = iConfig.getUntrackedParameter<double>("pfGammaEnergyScaleVar",1.);
    pfNeutralSmear_ = iConfig.getUntrackedParameter<bool>("pfNeutralSmear", false);
    doNaiveNeuPFScaling_ = iConfig.getUntrackedParameter<bool>("doNaiveNeuPFScaling", false);
    doRatioNeuPFScaling_ = iConfig.getUntrackedParameter<bool>("doRatioNeuPFScaling", false);
    doPeripheralNeuPFScaling_ = iConfig.getUntrackedParameter<bool>("doPeripheralNeuPFScaling", false);
    doCompensatoryNeuPFScaling_ = iConfig.getUntrackedParameter<bool>("doCompensatoryNeuPFScaling", false);
    TrackVariation_ = iConfig.getUntrackedParameter<double>("TrackVariation",0.);

    genPtMin_ = iConfig.getUntrackedParameter<double>("genPtMin",10);
    doSubEvent_ = iConfig.getUntrackedParameter<bool>("doSubEvent",0);
    doSubjetPurity = iConfig.getUntrackedParameter<bool>("doSubjetPurity",0);
    dopthatcut = iConfig.getUntrackedParameter<bool>("dopthatcut",0);
    doHardestSplitMatching_ = iConfig.getUntrackedParameter<bool>("doHardestSplitMatching",0);
  }
}

HiInclusiveJetSubstructure::~HiInclusiveJetSubstructure() { 
}

void HiInclusiveJetSubstructure::beginRun(const edm::Run& run, const edm::EventSetup & es) {
}

void HiInclusiveJetSubstructure::beginJob() {
  std::cout << "Running job with systematics" << std::endl;
  std::cout << "Doing Charged only: " << doChargedConstOnly_ << std::endl;
  std::cout << "Track efficiency var: " << TrackVariation_ << std::endl;
  std::cout << "Charged pfCand 4-mom var: " << pfChargedCandidateEnergyScale_ << std::endl;
  std::cout << "Neutral pfCand 4-mom var: " << pfNeutralCandidateEnergyScale_ << std::endl;
  std::cout << "Photon pfCand 4-mom var: " << pfGammaCandidateEnergyScale_ << std::endl;
  std::cout << "Neutral phi and rapidity smearing: " << pfNeutralSmear_ << std::endl;
  string jetTagTitle = jetTagLabel_.label()+" Jet Analysis Tree";
  t = fs1->make<TTree>("t",jetTagTitle.c_str());
  t->Branch("run",&jets_.run,"run/I");
  t->Branch("evt",&jets_.evt,"evt/I");
  t->Branch("lumi",&jets_.lumi,"lumi/I");
  if (useVtx_) {
    t->Branch("vx",&jets_.vx,"vx/F");
    t->Branch("vy",&jets_.vy,"vy/F");
    t->Branch("vz",&jets_.vz,"vz/F");
  }

  t->Branch("nref",&jets_.nref,"nref/I");
  t->Branch("jtptUncorrected",jets_.rawpt,"jtptUncorrected[nref]/F");
  t->Branch("jtEUncorrected",jets_.jtrawE,"jtEUncorrected[nref]/F");
  t->Branch("jtpt",jets_.jtpt,"jtpt[nref]/F");
  t->Branch("jteta",jets_.jteta,"jteta[nref]/F");
  t->Branch("jtphi",jets_.jtphi,"jtphi[nref]/F");
  if (doHiJetID_) {
    t->Branch("trackMax", jets_.trackMax, "trackMax[nref]/F");
    t->Branch("trackSum", jets_.trackSum, "trackSum[nref]/F");
    t->Branch("trackN", jets_.trackN, "trackN[nref]/I");
    t->Branch("trackHardSum", jets_.trackHardSum, "trackHardSum[nref]/F");
    t->Branch("trackHardN", jets_.trackHardN, "trackHardN[nref]/I");

    t->Branch("chargedMax", jets_.chargedMax, "chargedMax[nref]/F");
    t->Branch("chargedSum", jets_.chargedSum, "chargedSum[nref]/F");
    t->Branch("chargedN", jets_.chargedN, "chargedN[nref]/I");
    t->Branch("chargedHardSum", jets_.chargedHardSum, "chargedHardSum[nref]/F");
    t->Branch("chargedHardN", jets_.chargedHardN, "chargedHardN[nref]/I");

    t->Branch("photonMax", jets_.photonMax, "photonMax[nref]/F");
    t->Branch("photonSum", jets_.photonSum, "photonSum[nref]/F");
    t->Branch("photonN", jets_.photonN, "photonN[nref]/I");
    t->Branch("photonHardSum", jets_.photonHardSum, "photonHardSum[nref]/F");
    t->Branch("photonHardN", jets_.photonHardN, "photonHardN[nref]/I");

    t->Branch("neutralMax", jets_.neutralMax, "neutralMax[nref]/F");
    t->Branch("neutralSum", jets_.neutralSum, "neutralSum[nref]/F");
    t->Branch("neutralN", jets_.neutralN, "neutralN[nref]/I");

    t->Branch("eMax", jets_.eMax, "eMax[nref]/F");
    t->Branch("eSum", jets_.eSum, "eSum[nref]/F");
    t->Branch("eN", jets_.eN, "eN[nref]/I");

    t->Branch("eg_HFMax", jets_.eg_HFMax, "eg_HFMax[nref]/F");
    t->Branch("eg_HFSum", jets_.eg_HFSum, "eg_HFSum[nref]/F");
    t->Branch("eg_HFN", jets_.eg_HFN, "eg_HFN[nref]/I");

    t->Branch("h_HFMax", jets_.h_HFMax, "h_HFMax[nref]/F");
    t->Branch("h_HFSum", jets_.h_HFSum, "h_HFSum[nref]/F");
    t->Branch("h_HFN", jets_.h_HFN, "h_HFN[nref]/I");

    t->Branch("muMax", jets_.muMax, "muMax[nref]/F");
    t->Branch("muSum", jets_.muSum, "muSum[nref]/F");
    t->Branch("muN", jets_.muN, "muN[nref]/I");
  }
  if (doStandardJetID_) {
    t->Branch("fHPD", jets_.fHPD, "fHPD[nref]/F");
    t->Branch("fRBX", jets_.fRBX, "fRBX[nref]/F");
    t->Branch("n90", jets_.n90, "n90[nref]/I");
    t->Branch("fSubDet1", jets_.fSubDet1, "fSubDet1[nref]/F");
    t->Branch("fSubDet2", jets_.fSubDet2, "fSubDet2[nref]/F");
    t->Branch("fSubDet3", jets_.fSubDet3, "fSubDet3[nref]/F");
    t->Branch("fSubDet4", jets_.fSubDet4, "fSubDet4[nref]/F");
    t->Branch("restrictedEMF", jets_.restrictedEMF, "restrictedEMF[nref]/F");
    t->Branch("nHCAL", jets_.nHCAL, "nHCAL[nref]/I");
    t->Branch("nECAL", jets_.nECAL, "nECAL[nref]/I");
    t->Branch("apprHPD", jets_.apprHPD, "apprHPD[nref]/F");
    t->Branch("apprRBX", jets_.apprRBX, "apprRBX[nref]/F");
    t->Branch("n2RPC", jets_.n2RPC, "n2RPC[nref]/I");
    t->Branch("n3RPC", jets_.n3RPC, "n3RPC[nref]/I");
    t->Branch("nRPC", jets_.nRPC, "nRPC[nref]/I");

    t->Branch("fEB", jets_.fEB, "fEB[nref]/F");
    t->Branch("fEE", jets_.fEE, "fEE[nref]/F");
    t->Branch("fHB", jets_.fHB, "fHB[nref]/F");
    t->Branch("fHE", jets_.fHE, "fHE[nref]/F");
    t->Branch("fHO", jets_.fHO, "fHO[nref]/F");
    t->Branch("fLong", jets_.fLong, "fLong[nref]/F");
    t->Branch("fShort", jets_.fShort, "fShort[nref]/F");
    t->Branch("fLS", jets_.fLS, "fLS[nref]/F");
    t->Branch("fHFOOT", jets_.fHFOOT, "fHFOOT[nref]/F");
  }



  // t->Branch("jtsym",jets_.jtsym,"jtsym[nref]/F");
  // t->Branch("jtrg",jets_.jtrg,"jtrg[nref]/F");
  // t->Branch("jtdynkt",jets_.jtdynkt,"jtdynkt[nref]/F");
  // t->Branch("jtdyn_pt1",jets_.jtdyn_pt1,"jtdyn_pt1[nref]/F");
  // t->Branch("jtdyn_var",jets_.jtdyn_var,"jtdyn_var[nref]/F");
  t->Branch("jtdyn_split",jets_.jtdyn_split,"jtdyn_split[nref]/I");
  t->Branch("jtdyn_eta",jets_.jtdyn_eta,"jtdyn_eta[nref]/F");
  t->Branch("jtdyn_phi",jets_.jtdyn_phi,"jtdyn_phi[nref]/F");
  // t->Branch("jtdyn_theta",jets_.jtdyn_theta,"jtdyn_theta[nref]/F");
  t->Branch("jtdyn_deltaR",jets_.jtdyn_deltaR,"jtdyn_deltaR[nref]/F");
  t->Branch("jtdyn_kt",jets_.jtdyn_kt,"jtdyn_kt[nref]/F");
  t->Branch("jtdyn_z",jets_.jtdyn_z,"jtdyn_z[nref]/F");

  t->Branch("jt_intjet_multi", jets_.jt_intjet_multi,"jt_intjet_multi[nref]/I");
  t->Branch("jt_girth", jets_.jt_girth, "jt_girth[nref]/F");
  t->Branch("jt_girth_new", jets_.jt_girth_new, "jt_girth_new[nref]/F");
  t->Branch("jt_thrust", jets_.jt_thrust,"jt_thrust[nref]/I");
  t->Branch("jt_LHA", jets_.jt_LHA,"jt_LHA[nref]/I");
  t->Branch("jt_pTD", jets_.jt_pTD,"jt_pTD[nref]/I");
  if(doPrimaryLJPReco_){
    t->Branch("jt_PLJPkT",&jets_.jt_PLJPkT);
    t->Branch("jt_PLJPdR",&jets_.jt_PLJPdR);
    t->Branch("jt_PLJPeta",&jets_.jt_PLJPeta);
    t->Branch("jt_PLJPphi",&jets_.jt_PLJPphi);
  }
  t->Branch("triggerJetInAcceptance", &jets_.triggerJetInAcceptance, "triggerJetInAcceptance/O");
  // t->Branch("jtangu",jets_.jtangu,"jtangu[nref]/F");

  if(isMC_){
    if (useHepMC_) {
      t->Branch("beamId1",&jets_.beamId1,"beamId1/I");
      t->Branch("beamId2",&jets_.beamId2,"beamId2/I");
    }
    // Only matched gen jets
    t->Branch("refpt",jets_.refpt,"refpt[nref]/F");
    t->Branch("refeta",jets_.refeta,"refeta[nref]/F");
    t->Branch("refphi",jets_.refphi,"refphi[nref]/F");
    // t->Branch("refsym",jets_.refsym,"refsym[nref]/F");
    // t->Branch("refrg",jets_.refrg,"rg[nref]/F");
    // t->Branch("refdynkt",jets_.refdynkt,"refdynkt[nref]/F");
    // t->Branch("refangu",jets_.refangu,"refangu[nref]/F");
    // t->Branch("refdyn_pt1",jets_.refdyn_pt1,"refdyn_pt1[nref]/F");
    // t->Branch("refdyn_var",jets_.refdyn_var,"refdyn_var[nref]/F");
    t->Branch("refdyn_split",jets_.refdyn_split,"refdyn_split[nref]/I");
    t->Branch("refdyn_eta",jets_.refdyn_eta,"refdyn_eta[nref]/F");
    t->Branch("refdyn_phi",jets_.refdyn_phi,"refdyn_phi[nref]/F");
    // t->Branch("refdyn_theta",jets_.refdyn_theta,"refdyn_theta[nref]/F");
    t->Branch("refdyn_deltaR",jets_.refdyn_deltaR,"refdyn_deltaR[nref]/F");
    t->Branch("refdyn_kt",jets_.refdyn_kt,"refdyn_kt[nref]/F");
    t->Branch("refdyn_z",jets_.refdyn_z,"refdyn_z[nref]/F");
    
    t->Branch("jtdyn_isClosestToTruth", jets_.jtdyn_isClosestToTruth,"jtdyn_isClosestToTruth[nref]/O");
    t->Branch("refdyn_isClosestToReco", jets_.refdyn_isClosestToReco,"refdyn_isClosestToReco[nref]/O");
    t->Branch("jtdyn_refdyn_dR", jets_.jtdyn_refdyn_dR,"jtdyn_refdyn_dR[nref]/F");

    t->Branch("ref_intjet_multi", jets_.ref_intjet_multi,"ref_intjet_multi[nref]/I");
    t->Branch("ref_girth", jets_.ref_girth, "ref_girth[nref]/F");
    t->Branch("ref_girth_new", jets_.ref_girth_new, "ref_girth_new[nref]/F");
    t->Branch("ref_thrust", jets_.ref_thrust,"ref_thrust[nref]/I");
    t->Branch("ref_LHA", jets_.ref_LHA,"ref_LHA[nref]/I");
    t->Branch("ref_pTD", jets_.ref_pTD,"ref_pTD[nref]/I");
    if(doPrimaryLJPTruth_){
      t->Branch("ref_PLJPkT",&jets_.ref_PLJPkT);
      t->Branch("ref_PLJPdR",&jets_.ref_PLJPdR);
      t->Branch("ref_PLJPeta",&jets_.ref_PLJPeta);
      t->Branch("ref_PLJPphi",&jets_.ref_PLJPphi);
    }
    

    if(doSubjetPurity){
      t->Branch("refsub11",jets_.refsub11,"sub11[nref]/F");
      t->Branch("refsub12",jets_.refsub12,"sub12[nref]/F");
      t->Branch("refsub21",jets_.refsub21,"sub21[nref]/F");
      t->Branch("refsub22",jets_.refsub22,"sub22[nref]/F");
    }
    t->Branch("refparton_pt",jets_.refparton_pt,"refparton_pt[nref]/F");
    t->Branch("refparton_flavor",jets_.refparton_flavor,"refparton_flavor[nref]/I");
  }    

  if(doSubEvent_){
    t->Branch("subid",jets_.subid,"subid[nref]/I");
  }
  //declare asymmetry map here so we don't have to generate it every time the function is called
  float x[83] = {-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191};
  float y[73] = {-3.14159, -3.05433, -2.96706, -2.87979, -2.79253, -2.70526, -2.61799, -2.53073, -2.44346, -2.35619, -2.26893, -2.18166, -2.0944, -2.00713, -1.91986, -1.8326, -1.74533, -1.65806, -1.5708, -1.48353, -1.39626, -1.309, -1.22173, -1.13446, -1.0472, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.0872665, 0, 0.0872665, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, 0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.0472, 1.13446, 1.22173, 1.309, 1.39626, 1.48353, 1.5708, 1.65806, 1.74533, 1.8326, 1.91986, 2.00713, 2.0944, 2.18166, 2.26893, 2.35619, 2.44346, 2.53073, 2.61799, 2.70526, 2.79253, 2.87979, 2.96706, 3.05433, 3.14159};
  //jetasymmetry
  std::vector<float> TH2_contents = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.692481, 0.0889995, 0.0488732, 0.0271711, 0.591431, -0.0220031, 0.0499767, -0.0317, -0.222624, -0.00447521, -0.859672, -0.0507967, 0.139846, 0.0418028, 0.0564749, 0.0223388, -0.0239424, -0.0332173, -0.012542, -0.0251912, -0.0140944, -0.0053001, -0.0110252, 0.0119864, -0.130413, 0.0174637, 0.0130642, -0.0210857, 0.00365597, -0.00414784, -0.0109237, -0.00183268, -0.0145658, -0.0148098, 0.00420747, -0.0113168, -0.0100645, -0.00173105, -0.0063189, -0.0164269, -0.0155545, -0.0128474, -0.00517136, -0.0146287, -0.00406517, -0.00172483, -0.00630371, 0.00957061, 0.00255217, -0.00152076, 0.00632337, -0.00673325, -0.0116147, -0.0897954, -0.00527975, 0.000436002, 0.00347066, 0.0142754, 0.00451353, -0.0257207, -0.0096567, 0.0135829, 0.00858068, -0.0112624, -0.0306688, -0.100311, -0.00124851, 0.0615901, 0.0219068, -0.933002, 0.293214, -0.0318881, -0.0463794, -0.0171273, -1.09762, -0.755443, -0.419978, -0.860532, 0.105163, 0.346869, -0.201116, -0.283902, 0, 0, 0.57727, -0.00185421, 0.493765, 0.0441692, 0.0583022, -0.185069, 0.0485006, -0.0428545, -0.0599695, -0.752321, -0.0252789, 0.14263, 0.424159, 0.0654519, 0.0433135, -0.0161727, -0.0309603, -0.027013, -0.0139957, -0.0069596, -0.00872675, 0.0190906, 0.0243925, 0.0227801, -0.0199543, -0.00532035, 0.0134203, -0.0155735, -0.00815319, -0.0178682, -0.0149579, -0.00688743, 0.0105501, -0.011147, -0.00675572, -0.00443346, -0.0168502, -0.0126955, -0.00523446, -0.0153959, -0.00511822, 0.0043422, 0.00119243, -0.00569269, -0.0223246, -0.00946905, 9.70352e-05, -0.0042707, -0.00933542, 0.00574674, -0.00181296, -0.0168258, 0.00856457, -0.031274, 0.0106197, 0.00295752, 0.00272576, 0.00557588, 0.0164022, -0.00783731, 0.00392558, 0.0208606, 0.00963894, -8.32081e-05, -0.0165951, -0.0901057, -0.00162511, 0.0285461, -0.10459, -0.0967904, -0.303101, -0.0678063, -0.0952985, -0.0621147, -0.709681, 0.148665, -0.149715, 0.358851, -0.0288054, 0.0805165, -0.0912166, -0.041935, 0, 0, -0.094413, 0.0797872, -0.915948, -0.126334, -0.3275, -0.0586907, 0.297263, -0.245401, -0.0372844, -0.0314002, 0.246943, 0.120641, 0.17085, 0.12131, 0.0303635, -0.032442, -0.0239907, -0.0298847, -0.0507977, -0.0226592, -0.00945497, 0.00250452, -0.00235422, -0.0229798, 0.00236698, 0.00236396, 0.0227597, -0.016443, -0.00279459, 0.00517963, -0.00276472, -0.00696321, 0.00231039, -0.0063227, -0.0159179, -0.0125469, -0.0230273, -0.0142568, 0.00265222, 0.00188145, -0.0171193, -0.0264957, -0.0229323, -0.0168065, -0.0207187, -0.0139671, -0.0118489, -0.00241866, 0.0120056, -0.00702765, 0.00111993, -0.00851971, -0.0018665, -0.00926114, 0.0112867, 0.000936917, 0.00230619, 0.0131138, 0.0224393, -0.00111217, -0.0121739, -0.00270997, 0.0253676, 0.0047566, 0.00032666, -0.0358684, 0.0108563, 0.101476, 0.0762351, 0.242003, -0.335354, 0.18152, -0.0746725, 0.151479, -0.610246, 0.215392, -0.0743208, -0.149682, -0.380448, 0.0645457, 0.0353304, -0.0355584, 0, 0, -0.179057, -0.651616, -0.317413, -0.214704, 0.00569736, -0.454127, 0.0496812, -0.465714, -0.00488801, -0.0165198, -0.0149351, -0.0790763, -0.0483222, 0.0675668, -0.0177272, -0.0141974, -0.0199854, -0.0358977, -0.0349543, -0.0189815, -0.0141925, -0.000495171, -0.0151553, -0.0308546, -0.00959729, 0.0199475, 0.014882, 0.00237883, -0.00540407, -0.00418991, -0.00285692, 0.00167154, 0.0137326, 0.00405793, -0.0114689, -0.0075544, -0.00875609, -0.0118922, 0.00739902, 0.000903536, -0.00283981, 0.0101093, -0.0142075, -0.0168805, -0.0153334, -0.02036, 0.00760125, -0.00193258, 0.000315754, -0.00936756, -0.0125014, -0.0110745, 0.00472172, 0.00793685, 0.000756488, 0.0072499, 0.0200901, -0.00663225, -0.0061667, -0.0325506, 0.00165758, -0.0151671, 0.000412318, 0.00534114, -0.0350076, -0.0561563, -0.035415, 0.0298955, 0.0746828, -0.685617, -0.381398, 0.0586197, -0.0451635, -0.204068, -0.0998712, -0.79404, 0.256899, -0.315723, 0.0973528, -0.345976, -0.210872, 0.0928538, 0, 0, 0.241217, 0.400394, -0.0666695, -0.158005, 0.26225, -0.108582, -0.26941, -0.75131, -0.185876, -0.00700733, -0.0781575, -0.627763, 0.282856, 0.0821524, 0.026343, -0.00802542, -0.0306676, -0.0101172, -0.0149613, -0.0218905, -0.00243048, -0.00079285, 0.0156385, -0.0123283, -0.0177813, 0.0079354, 0.00546082, -0.000626886, 0.00683124, -0.0163611, -0.00351077, -0.00460977, -0.00092901, -0.00928307, -0.020251, -0.00884757, -0.01474, -0.022146, -0.000346573, -0.000936699, -0.0094029, -0.00333975, -0.0155691, -0.0165001, -0.0176222, -0.0213735, 0.00524076, -0.00549542, -0.00777862, -0.0107096, -0.0192286, 0.0131772, -0.00109536, 0.0114477, 0.0194206, -0.00260564, 0.0232004, 0.00314938, -0.0447919, 0.00589243, -0.00308165, 0.000139447, -0.00543521, -0.0103415, -0.0104856, -0.0535124, -0.00836283, 0.0448769, 0.0104794, -0.0388582, 0.0513463, 0.513094, 0.1955, 0.408535, -0.190217, -0.158773, -0.0106476, -0.05763, 0.137491, -0.447991, -0.174172, 0.640554, 0, 0, 0.635523, 0.383194, 0.0471924, -0.448492, -0.532818, 0.105994, -0.412574, -0.0209031, -0.18324, -0.464226, -0.696381, 0.353804, 0.401172, 0.108224, 0.0253444, -0.0163406, -0.0297756, -0.018598, -0.0192483, -0.0165515, 0.000276613, 0.00964332, 0.00763258, 0.0304398, 0.00667342, -0.036511, -0.0518301, -0.0128837, -0.0150542, -0.0276558, -0.00805501, -0.0121524, -0.019746, -0.0174639, -0.0264037, -0.0221524, -0.0223379, -0.0157773, -0.013965, -0.0295115, -0.0111453, -0.0080501, -0.0135597, -0.0201641, -0.0133529, -0.0272214, -0.0386694, -0.00988204, -0.0152411, -0.0113014, -0.00176733, -0.0174914, 0.00358105, -0.0106377, 0.00704942, 0.00696026, 0.00326909, -0.00515381, 0.0231363, -0.00652144, -0.00827459, 0.00428358, 0.00921997, -0.0235977, -0.0153305, -0.0609575, -0.0160713, 0.0401529, 0.090507, 0.503379, -0.00164319, 0.320254, -0.270456, -0.122933, -0.113325, -0.338322, -0.0922963, 0.178692, -0.0578076, 0.603073, -0.0767917, 0.0357519, 0, 0, 0.223967, 0.0820101, 0.0594998, -0.00696237, 0.38287, -0.31273, -0.155371, -0.0834628, -1.08342, -0.0309528, -0.0127956, 0.218202, 0.353339, 0.0952157, 0.015665, 0.0148236, -0.012928, -0.0929718, -0.0719128, -0.0237461, -0.000848476, -0.00272153, 0.0128599, 0.0102214, 0.00475587, 0.0079596, 0.00761983, -0.0160167, 0.0017728, -0.00526156, -0.00577721, -0.0134193, -0.0077349, -0.020991, -0.0123413, -0.0120822, -0.00713296, -0.00651997, -0.0163721, -0.125554, -0.0371654, -0.0043485, -0.019857, -0.0320028, -0.0176095, -0.0101254, -0.0256106, -0.0260723, -0.0182457, -0.0234095, -0.0243768, -0.0131801, -0.0175333, -0.0141, 0.0274373, 0.0258866, -0.0215948, -0.0591452, 0.0049807, 0.00476878, 0.0097065, 0.0245877, 0.020531, 0.0277221, -0.0167622, -0.0820686, 0.00408396, 0.0225291, 0.0840927, 0.302244, -0.144083, -0.0356616, -0.0688779, -0.0898273, -0.499005, -0.121706, -0.832438, -0.156388, -0.0635448, -0.0772618, 0.445295, -1.08198, 0, 0, 0.402111, -0.897993, -0.213741, -0.397141, -0.404678, -0.0333172, -0.361323, -0.122302, -0.543129, -0.0175922, 0.211651, 0.137641, 0.302261, 0.0733409, 0.0805957, 0.0132403, -0.00820685, -0.232306, -0.179226, -0.010014, -0.0101136, 0.00329426, 0.0110309, 0.0130043, 0.00826398, 0.00738476, 0.0173034, 0.00912987, 0.017987, 0.00950827, 0.00447975, -0.00186073, -0.0140161, -0.0276679, -0.00826424, -0.00777828, 0.00873003, 0.00114451, 0.00922456, -0.000657646, -0.00766935, 0.00144061, -0.0145528, -0.0146207, -0.0239155, -0.0169029, -0.00142693, -0.00331613, -0.00424719, -0.00569148, -0.000558325, 0.00542426, -0.0013854, 0.0104042, 0.0126233, 0.0025047, 0.021519, 0.00525557, 0.0445078, 0.00460212, -0.00499751, -0.00391799, 0.0265155, 0.00989075, -0.0145478, -0.0743591, -0.00388425, 0.0855553, 0.031261, 0.135778, 0.0136394, 0.0511522, -0.0444034, -0.051919, -0.215204, 0.246392, -0.101076, 0.158482, -0.817024, -0.202225, -0.101926, 0.421931, 0, 0, 0.778521, -0.081385, -0.0982709, 0.0892587, 0.443332, -0.00848676, -0.0465547, -0.600336, -0.0185422, -0.515981, -0.227618, 0.201082, 0.0266532, 0.0784823, 0.0434776, 0.0331402, -0.0246829, -0.100254, -0.101474, -0.0135935, 0.0223962, -0.000428967, -0.0090278, -0.00287544, 0.0107303, 0.0234042, 0.0173564, 0.0010533, 0.00892887, 0.0211536, 0.00535458, -0.00327053, -0.0172822, -0.00903424, -0.00486848, -0.00451604, -0.011407, 0.00635238, -0.00758598, -0.00337, -0.0174488, 0.00429635, -0.0112207, -0.0366266, -0.0187475, -0.00265132, -0.0103889, -0.0117039, -0.00999488, -0.00971019, 0.00435461, -0.00484959, 0.00390806, -0.0127897, 0.0228542, 0.00654182, 0.0166379, 0.0252971, 0.0231911, 0.0143937, 0.0186902, 0.0331305, 0.0028324, -0.00201826, -0.00620561, -0.0411134, 0.000333978, 0.0923436, 0.0213859, 0.178626, 0.752333, 0.000110625, -0.0989621, -0.190667, -0.0783878, -0.13253, -0.313936, 0.519067, -0.0122425, 0.0350773, -0.188402, -0.152124, 0, 0, -0.352215, -0.393448, 0.0488434, 0.125899, -0.0958107, 0.537927, 0.084583, -0.0967815, -0.0425838, -0.0400473, -0.0386545, 0.107311, 0.236004, 0.039269, 0.0436118, 0.00227775, -0.0177748, -0.0196556, -0.0223534, -0.00234434, -0.031381, -0.0217671, 0.0247793, 0.0205904, 0.0124014, 0.0085037, 0.0173181, -0.0133052, -0.0113499, -0.00106419, -0.00464989, -0.013887, -0.0403803, -0.0204153, -0.00936807, -0.015015, -0.0154599, -0.00624564, -0.0107153, -0.00675703, -0.0190415, -0.0177749, -0.021626, -0.0335676, -0.0196363, -0.0196789, -0.0297371, -0.0195992, -0.02468, -0.0311443, -0.00651764, -0.000699608, -0.00635447, -0.00806952, 0.0105818, 0.0118741, 0.00831456, 0.0081103, 0.0137638, 0.00919886, -0.00500086, 0.0120812, 0.0079366, 0.00645288, 0.00632571, -0.0596419, -0.00561277, 0.0640793, 0.0683026, 0.349955, -0.32298, -0.052906, -0.0659777, -0.328101, 0.10134, 0.250765, -0.313818, -0.545086, -0.831482, -0.0142488, -0.0794105, -0.0538377, 0, 0, 0.630796, 0.00407919, -0.0999201, -0.166125, -0.153519, 0.166715, 0.113169, -0.0753073, -0.0339575, -0.0102734, -0.661451, -0.242633, 0.229896, 0.0770016, 0.0374344, -0.00478962, -0.010084, -0.0570212, -0.0272859, -0.00482065, -0.0260493, -0.00294451, 0.00361729, 0.00393428, 0.0257346, -0.016381, 0.0118486, -0.0145691, 0.0028475, 0.00157045, -0.0109091, -0.0268666, -0.101676, -0.0502863, -0.0390606, -0.0205377, -0.0166666, -0.0110481, -0.00744403, -0.0110262, -0.0222995, -0.0282915, -0.0151853, -0.0286117, -0.0149096, -0.0302809, -0.031216, -0.0320598, -0.0246323, -0.0364602, -0.0192392, -0.00508147, -0.00189403, -0.0132494, 0.00858424, 0.000777411, -0.00515244, -0.00147004, 0.0352409, 0.0173434, -0.0116848, 0.00748003, 0.0176932, -0.00746976, -0.0165967, -0.0481361, 0.00436496, 0.0334366, 0.0491357, -0.0016776, 0.428081, -0.245933, -0.332742, -0.0389898, -0.0754877, -0.116659, -0.155858, 0.0526248, -0.0902852, -0.174906, 0.0986154, 0.0780115, 0, 0, 0.0294969, 0.710542, 0.246853, -0.181819, -1.10553, -0.108341, -0.0968914, -0.741393, -0.00804364, 0.00788821, 0.0199012, -0.0523333, 0.22182, 0.0722078, 0.0623054, -0.0407065, -0.0305878, -0.0124447, -0.0123877, -0.00569012, -0.0109109, -0.0133494, -0.00524988, -0.00894854, 0.00680346, 0.0234693, 0.0140284, -0.0174579, 0.00277419, 0.00426957, 0.010228, -0.0129847, -0.0392953, -0.000996321, -0.0180593, -0.0136793, -0.0113985, -0.00940063, 0.00572052, 0.00401037, -0.00748454, -0.0107947, -0.0192912, -0.0248065, -0.00352487, -0.00995329, -0.00697424, -0.0166744, -0.0256582, -0.0124119, -0.0166892, -0.00108947, 0.00739528, -0.00169881, 0.0161027, -0.00244309, 0.0229118, 0.00622133, -0.0203428, -0.0410362, -0.0104823, -0.0780274, -0.0177096, -0.00674507, -0.0169333, -0.0690751, -0.0111946, 0.0281637, 0.0726271, 0.105913, 0.167731, 0.0744736, -0.641828, -0.479662, 0.60511, -0.586241, -0.0654343, -0.0715008, -0.22106, -0.0760932, -0.031357, 0.467903, 0, 0, 0.249053, 0.989868, -0.00371472, -0.141392, -0.430668, -0.766347, -0.00336666, -0.0249267, -0.00970645, 0.262195, 0.00874484, 0.386975, 0.646159, 0.116328, 0.0315476, -0.00634763, -0.0238191, -0.0204939, -0.0276844, 0.00292402, 0.0157932, 0.0106556, -0.00328335, 0.0243152, -0.00583176, 0.00542975, -0.00972421, -0.0312602, 0.00356797, 0.011631, 0.00734189, -0.00254677, -0.0143934, -0.0191418, -0.00790896, -0.0234786, -0.00466583, -0.0178492, -0.00306653, -0.0137909, -0.0167064, -0.00790835, -0.0175689, -0.0134265, -0.0137238, 0.00076683, -0.0117192, -0.0128731, -0.000482551, -0.0136895, -0.00113186, 0.00894425, 0.000101945, -0.00995268, 0.00725163, 0.0177828, 0.0183539, 0.0136346, -0.0143643, -0.0339278, -0.0302801, -0.198482, -0.0757668, 0.0147628, -0.00472809, -0.0804901, -0.021419, 0.0231373, 0.0175283, 0.181176, -0.0499659, 0.441314, -0.456779, 0.00645369, -0.0479331, -0.997304, -0.480699, 0.278233, -0.0736059, -0.0587284, 0.08207, 0.181002, 0, 0, -0.389443, -0.148789, -0.858688, -0.370014, -0.370851, -0.541955, -0.0772077, -0.900191, -0.00302828, -0.0308162, -0.0369823, 0.349356, 0.233338, -0.0221938, 0.0443064, 0.00304124, -0.051543, -0.0194929, 0.00242344, -0.0127757, 0.00364768, 0.0140199, -0.00524501, -0.00768281, -0.016381, -0.0447044, -0.0229204, -0.0231399, -0.0443288, -0.0251841, -0.0135384, -0.0219644, -0.022706, -0.0240075, -0.0274266, -0.0528264, -0.0346044, -0.0127472, -0.0117596, -0.00941735, -0.0250471, -0.013843, -0.0338644, -0.0292246, -0.0217411, -0.0210789, -0.00930685, -0.0170774, -0.0120392, -0.0363605, -0.0253981, -0.0293458, -0.00649897, -0.00357375, 0.00474007, -0.0147378, -0.00763752, -0.016706, -0.0137936, -0.00500838, -0.000555554, -0.0299957, -0.0171872, -0.00967026, -0.0105729, -0.0802368, -0.0125737, 0.0395725, 0.0215275, 0.449801, -6.92833e-05, -0.437318, -0.0661687, -0.211653, -0.0842417, -0.63763, -0.280629, -0.503877, -0.230137, 0.166285, 0.00543485, 0.087686, 0, 0, 0.28185, -0.459855, 0.0796597, 0.238767, 0.659407, -0.0178574, 0.182376, -0.685928, 0.018133, -1.11309, -0.328435, 0.178425, -0.576355, -0.0394662, -0.01033, -0.0137273, -0.0532932, -0.0289949, -0.0416981, -0.0218728, -0.0361008, -0.0118195, -0.00716551, -0.0128887, -0.0180558, -0.0269275, -0.0125904, -0.0153188, -0.0272923, -0.0251575, -0.0152636, -0.0183252, -0.0348238, -0.0318612, -0.0424587, -0.0230369, -0.0147893, -0.0163467, -0.0300189, -0.022036, -0.0249379, -0.0227961, -0.0202033, -0.0104154, -0.0205427, -0.00752462, -0.0188932, -0.0216169, -0.020606, -0.0209418, -0.00528157, -0.0153525, -0.00268057, -0.000874482, 0.00244179, -0.0198601, -0.0350665, 0.0170332, -0.00583035, 0.0231644, 0.00584754, 0.0177185, 0.0223859, 0.000331274, 0.012227, -0.0636185, -0.0206553, 0.0292794, 0.0552782, -0.411339, -0.515819, -0.0259648, -0.0654042, -0.0376427, -0.996853, -0.715456, -1.00314, -0.390182, 0.0687867, 0.0101767, 0.060564, -0.0635338, 0, 0, 0.300993, -0.961547, -0.148008, -0.0498296, 0.166536, -0.0591029, -0.126643, -0.878997, -0.461792, -0.0250288, -0.905822, 0.0892704, -0.227084, 0.0277642, -0.0254591, -0.00735138, -0.0354956, -0.0288948, -0.0295874, -0.0302729, -0.0146654, -0.0105238, -0.0148108, -0.0154247, -0.00702647, 0.00511626, 0.00903628, -0.009179, -0.0284544, -0.0246433, -0.0449161, -0.0341853, -0.0163645, -0.0279609, -0.0346951, -0.0256351, -0.0212052, -0.0154316, -0.020518, -0.0164559, -0.00630044, -0.00260369, -0.00158553, -0.0157032, -0.00306127, -0.00253658, -0.0158594, -0.00674932, -0.0129851, -0.011132, -0.00540862, 0.0184149, 0.018082, 0.016687, 0.0159309, 0.0128742, 0.0139188, 0.0172888, 0.00404881, 0.0155121, 0.00747395, 0.00730698, 0.0212248, 0.00776203, -0.00183232, -0.0652221, 0.00213611, 0.00929812, 0.0198699, 0.216866, 0.131415, 0.0128178, -0.0662269, -0.0391224, -0.0461694, -0.0504201, -0.563149, -0.0175883, -0.0796148, -0.761598, 0.764126, -0.286025, 0, 0, -0.041431, -0.158913, -0.427602, 0.103955, 0.0942781, -0.0950267, 0.00520918, -0.0546388, -0.0297707, -0.0369032, -0.0455329, -0.0845437, -0.203323, 0.0608365, 0.0224457, -0.0261639, -0.0174784, -0.0481156, -0.0453423, -0.0274567, -0.00280218, 0.012375, 0.0140053, 0.00821749, -0.0030262, -0.00434924, -0.00107374, 0.00301484, -0.0215459, -0.0381095, -0.0558011, -0.121491, -0.0249246, -0.0229143, -0.0368204, -0.0375643, -0.0274013, -0.0273612, -0.0181671, -0.0212474, -0.0265247, 0.00697195, -0.0219561, 0.0061888, -0.0207724, -0.00725199, -0.00749185, -0.00316676, -0.00453047, -0.00599191, 0.00352235, 0.0129495, 0.0201342, 0.00365264, 0.0220816, 0.0045046, 0.0185548, 0.0403589, 0.0349068, 0.0068003, 0.00928776, 0.0019501, 0.00997909, 0.00420352, 0.00169043, -0.0874413, 0.00215645, 0.0409604, 0.0939347, 0.0201855, 0.240988, -0.04507, -0.293279, -0.0325009, 0.0234813, -0.270725, -0.0593942, -0.0852383, 0.167657, 0.226025, 0.0477906, 0.68768, 0, 0, 0.212925, 0.00167504, -0.268522, -0.556555, -0.277276, 0.130338, 0.0221382, 0.405686, -0.0102554, -0.0164467, -0.0237906, -0.694015, 0.56998, 0.0565794, -0.0109281, -0.0585175, -0.0645067, -0.038687, 0.000111035, -0.0100362, -0.0195029, 0.00387161, 0.0342483, 0.00986556, -0.0302034, -0.0278349, -0.0166448, -0.0153854, -0.0404715, -0.039173, -0.0542594, -0.051004, -0.0317219, -0.040371, -0.0407181, -0.0206642, -0.0450981, -0.0385346, -0.021636, -0.0226697, -0.0142997, -0.019271, -0.0208464, -0.0276896, -0.0225258, -0.00838691, -0.014887, 1.42137e-05, 0.0100338, -0.00318235, -0.00383148, 0.00547906, 0.0138756, 0.0211583, 0.0149498, 0.0168732, 0.00747056, 0.0283747, 0.00807564, -0.00228128, 0.00342086, -0.0110724, 0.00274869, -0.00908247, -0.0231871, -0.0794996, 0.00706185, 0.0535608, 0.0915425, 0.230224, -0.596968, -0.471079, 0.0620554, -0.0171646, 0.312551, -0.0717306, -0.871038, 0.391734, -0.356247, 0.0118207, 0.393514, -0.0694183, 0, 0, 0.128511, 0.142362, 0.471126, -0.00783201, -0.205377, -0.136333, -0.0787576, -0.701425, -0.510825, 0.182665, -0.327118, -0.43997, 0.217956, -0.00233268, 0.0404661, -0.104182, -0.172576, -0.0504134, -0.0176285, -0.00773804, -0.0166695, -0.0298901, 0.0183, 0.0306251, -0.0165226, -0.0291249, 0.00234155, -0.0101982, -0.0186401, -0.0219169, -0.0304088, -0.0310161, -0.0124325, -0.0196556, -0.0233638, -0.0118846, -0.0328442, -0.0274516, -0.00290547, -0.0171366, -0.0177174, 0.00468367, -0.0413558, -0.0341101, -0.00367526, -0.00869261, -0.0208481, -0.0162911, -0.0119667, -0.00759023, 0.00797415, -0.00198518, 0.00330011, 0.00776892, 0.0234058, 0.012583, 0.021571, 0.0247238, -0.102833, -0.0763228, 0.0139222, -0.00349591, 0.0100509, -0.0133588, -0.0274278, -0.0802441, -0.0272539, 0.0241585, 0.101616, 0.263918, 0.0972339, -0.606801, -0.0889712, -0.667896, 0.0218763, -0.455695, -0.306789, 0.00350781, -0.155689, -0.46357, -0.0947435, -0.102665, 0, 0, 0.539175, 0.087963, 0.627938, -0.700615, -0.0902762, -0.0549459, -0.0860612, -0.00533496, 0.050373, -0.465637, 0.210664, 0.377927, 0.284223, 0.00642663, -0.0150038, -0.0876824, -0.328577, -0.046662, -0.0244414, -0.0212319, -0.0267685, 0.00303327, 0.0176938, 0.0323879, 0.00936426, 0.00458806, 0.00354627, -0.00486554, -0.00787692, -0.00477254, -0.0116544, -0.0125197, -0.00572854, -0.0168636, -0.00196477, -0.0133117, -0.0128892, -0.0203308, -0.00480916, -0.0169551, 0.00420686, -0.00993791, -0.00550043, -0.0273383, -0.0108454, 0.00888238, 0.00269615, -0.000268578, 0.00012457, -0.00793111, 0.0179279, 0.0236951, 0.0246179, 0.00206328, 0.0211937, 0.0268581, 0.0522184, 0.016452, -0.012135, -0.00858795, 0.00120803, -0.011006, -0.0587442, -0.0354608, -0.0414513, -0.0565358, -0.0169749, 0.00102747, 0.112798, -0.0320465, 0.534518, 0.0943536, -0.30273, -0.0437734, -0.934278, -0.109102, -0.155802, -0.0227225, -0.434801, -0.0531123, 0.146039, 0.439123, 0, 0, -0.0942964, 0.261633, 0.0482, 0.236219, -0.19691, 0.1151, -0.733259, -0.0243337, -0.834529, 0.00570497, 0.187009, -0.795929, 0.472999, 0.0252695, 0.0121947, -0.0518501, -0.136402, -0.0117358, -0.00976623, 0.00696897, 0.00902133, 0.000441746, 0.0217002, 0.0124912, 0.00485659, 0.0113356, -0.00829011, -0.0193383, -0.0204947, -0.00608662, -0.0190126, -0.0146416, -0.00941714, -0.00827809, -0.0170268, -0.0176692, -0.0270147, -0.0147372, -0.00734007, -0.0112941, -0.00641252, -0.00245783, -0.00537172, -0.00587149, -0.00869482, -0.00435832, 0.0107919, -0.0167084, -0.00127327, 0.0165896, 0.0192186, 0.0130747, 0.0285815, 0.0163599, 0.0230127, 0.0155468, 0.0245412, -0.0147325, 0.00855763, 0.0268653, 0.0186923, -0.0310038, -0.18751, -0.0592252, -0.0468672, -0.102733, 0.00147867, 0.0863372, 0.042922, -0.403526, 0.578581, -0.144902, -0.542699, -0.305707, -0.210116, -1.08386, -0.00555899, -0.347914, -0.00430652, 0.144005, 0.114113, -0.0127881, 0, 0, 0.188091, -0.674474, -0.148973, 0.0373107, 0.0899469, -0.0795542, -0.0179742, 0.316678, -0.172253, 0.00261324, 0.00124863, 0.393631, -0.565922, -0.0401476, 0.016973, -0.0108418, -0.0171418, -0.00225157, 0.0108706, 0.00523828, -0.0108039, 0.0134939, 0.014142, -0.00790292, -0.0166772, -0.0265912, -0.0241853, -0.0513839, -0.0477891, -0.00696847, -0.0354327, -0.0166183, -0.0232465, -0.0166491, -0.0171298, -0.0264992, -0.0244162, -0.0256434, -0.0183537, -0.0188171, -0.0193965, -0.029769, -0.0229885, -0.0259827, 0.00409849, -0.00690733, -0.0136515, -0.019543, -0.0229246, -0.0250224, -0.0210944, 0.000573057, -0.00510287, 0.00110229, 0.0190815, 0.00185917, -0.00410449, -0.0692421, -0.0286555, 0.0100566, 0.0244438, 0.00340259, -0.0822208, -0.0295891, -0.0396404, -0.0811885, 0.00380186, 0.0211173, 0.0284762, -0.10843, -0.160239, -0.00675668, -0.333202, -0.012584, -0.0372156, 0.159653, -0.16331, -0.0391159, -0.166875, -0.216615, 0.969372, -0.0461182, 0, 0, -0.0445588, 0.0722387, -0.466244, 0.062515, -0.0229374, -0.0802714, -0.0392472, -0.247969, -0.878215, -0.0260774, 0.301746, -0.017404, 0.051848, -0.331488, 0.0312752, 0.00441938, -0.00683359, -0.00092288, 0.00231327, -0.00763841, 0.0052679, 0.022876, 0.00340332, -0.000378609, 0.00333589, -0.0105034, -0.0251821, -0.101806, -0.0173183, -0.000750721, 0.00194591, -0.00763061, -0.00224815, -0.00446346, -0.000617503, -0.00666703, -0.0311539, -0.0306042, -0.0249979, -0.0266558, -0.0142911, -0.00787967, -0.0185795, -0.0251786, -0.00914451, -0.0231562, -0.0132076, -0.0347546, -0.0241394, -0.0215615, -0.0257452, -0.0896875, -0.00652265, -0.00375718, -0.0098543, -0.00796308, -0.0146076, -0.0213663, 0.0224742, 0.0179985, 0.0186189, 0.0175267, 0.0103381, 0.0150635, 0.00386416, -0.0642843, -0.00993035, 0.0283, 0.0345942, 0.456871, 0.512384, -0.024037, -0.0521942, -0.0159726, -0.0591839, -0.147348, -1.127, -1.39899, -0.0253397, -0.0437733, -0.171457, 0.345886, 0, 0, 0.514996, 0.561644, -0.205406, -0.0445193, -0.0564601, 0.134198, -0.806295, 0.0604271, 0.01829, 0.390995, -0.0965297, -0.494484, -0.0388082, 0.0366032, -0.00247165, 0.00181455, -0.0201387, -0.0138291, -0.0116446, -0.00831171, -0.0238302, 0.0112975, 0.00858317, 0.0101952, -0.0150952, 0.0170953, 0.0167433, -0.0139944, 0.0157044, 0.0179192, -0.00277275, 0.0125672, 0.00537864, 0.015491, -0.0175429, -0.0191748, -0.0124961, -0.00243863, 0.00482448, 0.0135138, -0.0114904, 0.0076344, -0.00892325, -0.0127097, -0.0236242, -0.0166666, 0.00205768, -0.0148039, -0.0097373, -0.00245467, -0.035396, -0.116932, -0.00303728, 0.00855751, 0.00829363, -0.0036094, 0.041409, -0.00655241, 0.0350246, 0.0317059, 0.00833973, 0.0386699, 0.0208045, -0.00801076, -0.000897668, -0.0513927, 0.00387146, 0.0852568, 0.0488904, 0.294637, 0.201512, -0.000120289, -0.393394, -0.0157558, -0.00244999, -0.402398, -0.814856, -0.131263, -0.145991, 0.304087, 0.198481, 0.475031, 0, 0, 0.487087, 0.0603728, -0.883903, 0.068626, -0.100645, -0.710072, -0.0384402, 0.0261081, 0.231246, -0.966398, 0.0382473, 0.277781, 0.318576, 0.0584894, 0.0323774, -0.00584073, -0.0237864, -0.00916544, 0.00955872, -0.0194183, 0.00175243, 0.0110246, 0.00400472, 0.0176626, -0.00276332, 0.0235265, 0.0244138, 0.025107, 0.0211633, 0.0249371, 0.0192607, 0.00114321, -0.0030998, 0.00603173, -0.00150025, -0.0214698, -0.00877634, 0.0050816, -0.00102392, -0.00649798, -0.00891615, -0.0024796, -0.0218676, -0.00332546, -0.0083349, -0.00132874, -0.00977819, 0.00373243, -0.0246856, -0.00864188, 0.00307244, -0.027256, -0.0050129, -0.00911769, 0.0172625, -0.0076875, 0.0167398, 0.0113858, -0.00919036, 0.00599864, -0.00648537, 0.00322963, 0.00456096, -0.0142341, 0.00246186, -0.111221, -0.0421933, 0.064048, 0.0805112, -0.132208, 0.167088, -0.159077, 0.403968, 0.0588063, -0.107258, 0.0513844, -0.465332, 0.1714, 0.649241, -0.921967, 0.105163, 0.366365, 0, 0, -0.469833, -0.370531, -0.642967, 0.135673, 0.288769, -0.0547151, 0.377978, 0.105099, -0.0329028, -0.0231325, 0.176367, 0.027945, -0.461634, 0.0230517, 0.016432, 0.00309921, -0.0185723, 0.0193853, -0.00342293, 0.00743381, 0.0365908, -0.00916713, 0.00687202, 0.00677377, -0.00434318, -0.00383762, 0.0296625, 0.0250713, 0.013498, 0.0163719, -0.00432996, -0.0151757, -0.00412119, -0.00233308, -0.0149814, -0.0195038, -0.018703, -0.00779119, 0.00134584, -0.0105926, -0.0120387, -0.00582776, -0.00559996, -0.0187136, -0.00618335, -0.00631535, -0.00266539, -0.0131128, -0.00240585, -0.00651355, -0.00809424, 0.00592583, -0.000794058, 0.00295955, 0.0206501, 0.0204275, -0.000720609, 0.0182077, 0.0262744, 0.0213577, 0.0194451, 0.0222895, 0.00829639, -0.00560416, -0.00772383, -0.0921394, 0.00353198, 0.0639885, 0.128508, -0.80228, 0.145663, 0.00329372, 0.269778, -0.0289202, -0.170571, 0.096877, -0.26603, 0.294007, -0.281785, 0.168197, 0.0768515, -0.0717519, 0, 0, -0.160468, -0.167454, 0.0722472, -0.523188, 4.54279e-05, -0.146478, 0.317777, -0.16829, 0.0389416, -0.0118359, -0.0276925, 0.799409, 0.231584, 0.0798621, 0.0132745, -0.0121159, -0.0507466, -0.00242222, 0.00435814, -0.003323, -0.0269588, -0.00261544, -0.00668584, -0.0076098, 0.00762192, 0.00354107, 0.00362186, 0.0187954, 0.0133044, -0.00639171, -0.00646299, -4.74649e-05, -0.0501537, -0.011273, -0.00566724, -0.0127273, -0.0225186, -0.0163958, -0.0095025, -0.0168432, -0.0135319, 0.0158495, -0.0128207, -0.00989353, 0.000543495, 0.00821998, -0.0130646, -0.0183809, -0.00446084, -0.00931883, 0.00353106, 0.0161154, 0.0077575, -0.00348148, 0.0251439, 0.0263258, 0.0109449, 0.0249446, 0.0293472, 0.00663873, 0.0141063, 0.0225947, 0.0101636, 0.023357, -0.00447109, -0.0597867, 0.0250747, 0.10092, 0.0829268, 0.268094, 0.584376, -0.00277389, -0.0505228, -0.0480325, -0.0877994, 0.207063, 0.371125, -0.338053, 0.0426316, -0.353855, -0.271876, 0.632716, 0, 0, 0.412453, -0.102212, 0.229594, -0.22968, -0.311309, -0.121611, 0.0202696, -0.033427, -0.0184133, -0.217968, 0.0138354, 0.125228, -0.583038, 0.0679278, 0.0342704, -0.0320054, -0.00759838, 0.00381854, 0.002775, 0.0207677, 0.0191766, 0.013438, 0.000236473, -0.00328369, -0.00738024, 0.0102559, -0.0215117, -0.0102612, 0.00551281, 0.00353404, -0.0219512, -0.0359688, -0.128106, -0.0201671, -0.0208436, -0.0125701, -0.0051259, -0.0138929, 0.00123723, -0.0046848, -0.00363812, -0.00292529, -0.0146595, -0.0200346, -0.00489605, -0.00161395, 0.00209404, -0.00183548, -0.00854808, 0.0156032, 0.015115, 0.0257162, 0.0221907, 0.0327419, 0.0382704, 0.0319965, 0.0477585, 0.0428388, 0.0286873, -0.000407853, 0.0166121, 0.0130654, 0.0325611, 0.0186069, 0.00599734, -0.0564511, 0.0301147, 0.0528055, 0.0713452, -0.102496, -0.216865, -0.0305222, 0.469911, 0.00430173, -0.0643096, -0.095333, -0.0427622, -0.106464, -0.0057675, -0.153411, 0.0228727, 0.465027, 0, 0, 0.490903, 0.010703, -0.752234, -0.504479, 0.180321, -0.0485829, -0.138849, -0.433403, 0.199836, 0.265529, 0.175806, 0.17827, 0.167614, 0.080165, 0.0255901, 0.0227498, 0.0087754, -0.0155624, 0.00986508, 0.036745, -0.00847819, 0.0259095, 0.0101354, 0.0333574, 0.0135552, 0.0186992, 0.00975935, -0.00324023, -0.00583006, -0.00470353, -0.00799836, -0.0234486, -0.042102, -0.0235798, -0.0213232, -0.0165996, -0.0153903, -0.0181036, -0.0127498, -0.00283848, -0.0166476, 0.009912, -0.0141243, 0.00333999, -0.00667299, -0.000324977, -0.00458827, 0.00894143, 0.0118349, 0.00970227, 0.0244991, 0.034225, 0.0306877, 0.0120308, 0.032865, 0.0168842, 0.0364536, 0.0513214, 0.0404303, -0.0100964, -0.0889457, -0.048743, 0.00334361, 0.00501566, -0.000210375, -0.0431099, 0.0284169, 0.0553588, 0.0531871, 0.164237, 0.252138, -0.00169722, -0.0196325, 0.0169965, 0.477422, -0.0940251, -0.350764, 0.12493, 0.120849, 0.255756, -0.197851, 0.326726, 0, 0, -0.0784585, -0.236337, -0.0839466, -0.75086, -0.209065, -0.163513, -0.112358, -0.0578356, -0.0266611, -0.0130958, -0.794319, 0.0898842, 0.29911, 0.0477039, 0.0198492, 0.0064856, 0.00949076, -0.00311521, 0.018277, -0.00321821, 0.0185347, 0.00702986, 0.00307391, 0.0367352, 0.0281243, -0.00730909, 0.00482063, -0.0144217, 0.00248845, -0.013463, -0.0196508, -0.0148965, -0.0297674, -0.021188, -0.03971, -0.0298531, -0.0251435, -0.0333612, -0.0198779, -0.0233053, -0.0136967, -0.00975937, -0.00581612, -0.0288246, -0.00726462, 0.00466704, -0.0147616, -0.0121361, 0.00742774, 0.0038624, -0.00647247, -0.00128619, 0.0100818, -0.00402962, 0.0104923, 0.0166883, 0.00556791, 0.0561834, 0.0346407, 0.000394754, -0.139566, -0.0932976, 0.0117126, 0.0158053, 0.0110353, -0.0943963, 0.0197573, 0.0897992, 0.0872092, 0.190855, -0.0376771, 0.167528, -0.0378886, -0.248822, -0.0460158, -0.0873755, 0.370288, 0.308448, 0.0418887, 0.228399, 0.162243, -0.0243572, 0, 0, -0.764274, -0.0353231, -0.133929, -0.216805, -0.0186686, 0.146013, -0.632863, -0.133204, 0.242889, -0.0349871, -0.298309, 0.143928, 0.472838, 0.0517977, 0.0453507, 0.0235543, -0.0235899, -0.0248391, -0.00324955, -0.00413335, -0.016879, -0.0112718, -0.019926, 0.0407289, 0.0225074, 0.0168289, 0.0220807, -0.00396754, 0.016518, -0.0249675, -0.0278254, -0.00740648, -0.0106253, -0.00674725, -0.0200977, -0.0198226, -0.0271197, -0.0178634, -0.0137568, -0.0110771, -0.0081122, 0.00800865, -0.021834, -0.0234804, -0.00739671, -0.0129857, -0.00294249, -0.0124249, -0.00199898, -0.00773251, -0.0034825, 0.000357097, 0.00375253, 0.0123657, 0.000886265, -0.142164, -0.0360362, 0.0291989, 0.0535296, -0.000693916, -0.0130461, 0.0120733, 0.0170278, 0.0043947, -0.0191736, -0.0882283, 0.0115019, 0.0787707, 0.0851559, 0.229709, 0.390121, -0.00913994, 0.0947911, -0.160598, -0.0669977, -0.00140256, -0.153701, 0.387987, 0.226922, 0.181253, 0.0975307, 0.314734, 0, 0, -0.0373773, -0.469394, 0.117187, 0.335853, 0.0650852, -0.818217, -0.0947168, -0.927015, -0.0261038, 0.0405434, -0.131505, -0.232911, 0.292449, 0.110209, 0.0500038, -0.0260945, -0.0181653, 0.00429683, -0.00827872, -0.0165463, -0.00559769, -0.034869, 0.00175833, 0.00119155, -0.00418814, 0.0289915, 0.0289137, 0.0017992, -0.0121535, -0.0376676, 0.00343028, 0.0068486, 0.00257965, 0.0069954, -0.00566578, -0.0146335, 0.00191801, -0.013504, -0.0176088, -0.014067, -0.00601168, 0.0158617, -0.00223197, -0.00256362, 0.00432481, 0.00824641, 0.00502719, 0.0188477, 0.00944534, -0.00290501, 0.00844489, 0.0216324, 0.0158925, 0.0267472, 0.0205894, -0.032283, 0.0284659, 0.0206014, 0.0125566, 0.0231035, 0.0242407, 0.0315298, 0.0413509, 0.00573452, -0.0125175, -0.0838344, -0.00150454, 0.0433353, 0.0676132, 0.260597, 0.484359, 0.0118737, -0.218791, -0.00237501, -0.0266996, -0.624673, -0.965155, -0.929037, 0.115582, -0.137254, 0.0205733, 0.362542, 0, 0, -0.137627, 0.130039, 0.190539, 0.00270871, 0.119092, -0.0352038, 0.0902964, -0.16766, -0.387323, -0.00494962, -0.0144727, 0.167822, -0.0334478, 0.031537, 0.060517, -0.0174261, 0.00603633, -0.0106068, -0.00588314, 0.00323905, 0.00969178, -0.00479701, 0.000866037, 0.00145422, 0.044815, 0.0253685, 0.023114, 0.0096506, 0.0126632, 0.00278587, -0.00442133, 0.00283768, -0.00172954, -0.0113011, -0.0218328, -0.0181368, -0.00220104, -0.0159413, -0.00859192, -0.00122952, 0.00455762, 0.0094058, 0.00967053, 0.000454905, -8.09632e-05, 0.00781213, 0.0155594, 0.0353065, 0.011578, 0.00196034, 0.0356725, 0.0222273, 0.030128, 0.0434793, 0.017285, 0.0268347, 0.0393556, 0.0107186, -0.0115325, -0.0107037, 0.0206162, 0.0183161, 0.0148265, -0.00808365, -0.00323057, -0.071564, -0.0115766, 0.0378367, 0.0900356, 0.363974, -0.0184806, 0.116281, -0.0510033, -0.927388, -0.0664264, -0.0667923, -0.458808, 0.30189, 0.0123816, -0.200028, 0.0947212, 0.4383, 0, 0, -0.352443, -0.0271867, -0.851348, -0.229274, -0.359254, -0.686877, -0.0781748, -0.0513121, 0.076785, -0.0427866, -0.350012, 0.63066, 0.305153, 0.0192058, 0.0422972, 0.0175127, -0.008154, 0.0111108, 0.0107197, 0.0207262, 0.0271889, 0.00402334, 0.0275648, 0.022447, 0.0469009, 0.00995687, 0.0341901, 0.00194895, -0.00151019, -0.0117528, -0.00804144, -0.00996144, -0.000955305, -0.00358441, -0.0209711, -0.00208539, -0.006522, -0.0121445, -0.00569893, -0.0103581, -0.0015039, -0.00945226, -0.000208602, -0.009544, -0.00461673, 0.0019743, 0.0108522, 0.000495774, 0.00222119, -0.00688415, 0.0172314, 0.0259007, 0.0156379, 0.00632626, 0.0233718, 0.0339504, -0.00692452, -0.0934351, 0.00444249, 0.0188699, 0.0130893, 0.034656, 0.0216284, 0.0112667, -0.00725486, -0.0558914, 0.0147992, 0.0946649, 0.122226, 0.251983, 0.425912, -0.41859, -0.285137, -0.0211333, -0.0581453, -0.46358, -0.28841, 0.168473, -0.00401958, -0.0375516, -0.0437869, -0.0362788, 0, 0, -0.0645841, 0.129516, -0.000334656, -0.126008, -0.016745, -0.248455, -0.71293, -0.751575, 0.0402889, -0.0486334, -0.204228, -0.4821, -0.0776748, 0.0158469, -0.0172952, -0.00318626, -0.0154577, -0.0196312, -0.0250324, 0.00310171, 0.00205304, -0.00686807, 0.0171125, 0.0430856, 0.0385478, -0.000123976, 0.0188028, 0.00641699, 0.00899384, 0.00161723, -0.00526521, -0.00218792, 0.00872821, -0.0104607, 0.00327255, -0.00662007, 0.00548795, -0.00233601, 0.0160168, 0.0117109, -0.0100239, 0.0013978, -0.0130098, -0.0153823, -0.00726189, -0.00261646, -0.00792061, 0.00120774, 0.0151148, 0.00104446, 0.0115268, 0.00172457, 0.0293624, 0.0198722, 0.0196464, -0.00540982, 0.000527376, 0.00499695, 0.0704101, 0.0316748, 0.0269637, 0.021452, 0.0230729, 0.00855679, -0.0090394, -0.0697754, -0.0152994, 0.0523278, 0.099182, 0.522185, 0.14118, 0.423743, -0.287687, -0.0452139, -0.117948, 0.29623, -0.107898, -1.25971, -0.361605, 0.330154, -0.0503144, -0.0867547, 0, 0, 0.567689, 0.364659, -0.827, -0.336805, -0.0440169, 0.194644, 0.156453, -0.399852, 0.0959205, -0.208557, 0.0154619, 0.277865, 0.246065, -0.0944117, 0.0132952, -0.00410205, -0.0309, -0.0382915, -0.0088056, 0.00648376, -0.00919404, 0.0142526, 0.0109588, 0.0465022, 0.0420259, 0.0556818, 0.0264411, 0.0152993, 0.015819, 0.0234409, 0.0198746, 0.0215849, 0.0059685, 0.0246895, 0.0065808, 0.00788566, 0.0152905, 0.0141323, 0.0144586, 0.0109033, -0.00700836, 0.00112261, -0.00434008, -0.0193208, -0.00684244, -0.0063874, -0.00301004, 0.0193706, 0.0236506, 0.0286668, 0.0228308, 0.0142812, 0.0475344, 0.0292229, 0.0309475, 0.0345772, 0.0462741, 0.0380922, 0.0742386, 0.0339316, 0.0239275, -0.00352686, 0.00189106, 0.00911325, -0.0306712, -0.0749828, 0.0085082, 0.03496, 0.109663, 0.227144, 0.479447, 0.0344832, -0.432594, 0.329296, -0.0484913, -0.129035, -0.0289525, -1.0997, -0.331487, 0.0396423, 0.0188143, 0.267985, 0, 0, 0.474068, 0.565743, 0.291803, 0.127826, -0.340298, 0.00850821, -0.038647, -0.0220303, -0.00341015, -0.281647, 0.11739, 0.179101, 0.271529, 0.0588565, 0.062436, -0.0137906, -0.0125144, -0.0386517, -0.0117843, -0.00301774, 0.0168448, 0.0124584, 0.0109708, 0.0317212, 0.0398139, 0.0198839, 0.0387523, 0.0243939, 0.0189793, 0.0131817, 0.0152631, 0.0106504, 0.0108413, 0.00846334, 0.00553531, -8.81631e-05, 0.00572586, -0.00381525, 0.00333184, -0.000931626, -0.0105723, -0.000518522, -0.00977453, 0.00400062, 0.0082163, 0.0169763, 0.00836601, 0.0197233, 0.00422216, 0.0174541, 0.0352514, 0.0492317, 0.0339615, 0.025172, 0.0362551, 0.0385417, 0.0489231, 0.0606944, 0.045074, 0.0413339, 0.028913, 0.0118882, 0.010215, 0.0248234, -0.0431525, -0.0527368, 0.0655685, 0.0319777, 0.0751326, 0.0576953, 0.258106, -0.547353, -0.57591, -0.00322955, -0.0300857, -0.0489761, -0.045763, 0.0934043, -0.0213741, -0.483446, 0.574264, 0.604101, 0, 0, 0.468868, -0.255142, 0.603, 0.302396, 0.0445697, -0.0926592, -0.0954019, -0.0477437, 0.149467, -0.0652414, -0.474575, 0.200247, 0.058097, 0.0435526, 0.0506328, 0.03344, -0.0139223, -0.0167207, -0.011589, -0.00726116, 0.00135071, -0.00240937, 0.0148891, 0.0560764, 0.0323766, 0.0134203, 0.0293246, 0.0068901, 0.0176646, 0.00468647, 0.00672726, 0.00450808, 0.00127631, 0.00376283, -0.00163642, 0.00440259, -0.0128817, -0.0626749, -0.0966205, -0.0186973, -0.020446, -0.0128209, 0.0102599, -0.0131158, 0.00186081, -0.00553744, 0.0193822, 0.00764306, 0.00235896, 0.00811266, 0.0100111, 0.0177341, 0.0247059, 0.0450321, 0.0393368, 0.0222615, 0.0261999, 0.0493461, 0.0531924, 0.0213252, -0.003802, 0.0106541, 0.0153315, -0.011353, -0.00671507, -0.0837262, 0.000950656, 0.0674277, 0.0169532, -0.130655, 0.189444, -0.0849854, -0.0611762, -0.63774, -0.0940997, -0.0832936, -0.806225, 0.417756, -0.20905, -0.0170449, -0.0271846, 0.038101, 0, 0, 0.256711, 0.0347792, 0.0689306, -0.316921, -0.114573, 0.029973, -0.389396, 0.269678, -0.231083, -0.00521652, 0.424923, 0.0724808, -0.0132379, 0.0733719, 0.0333046, 0.000477433, -0.0296686, 0.00641167, 0.00144191, -0.000704228, 0.00387595, -0.0164826, 0.00398204, -0.00417994, 0.00719027, 0.0244619, 0.0194756, 0.0136711, 0.0410231, 0.0183037, 0.00331396, 0.0129695, 0.00910799, -0.000376153, -0.0034614, -0.0155515, -0.0503086, -0.143599, -0.0252458, -0.0128705, -0.0042836, -0.00922437, -0.00507965, -0.0127228, 0.00665546, -0.0134515, 0.004969, 0.0122294, 0.00365113, -0.0191163, -0.00337505, 0.0376108, 0.028305, 0.0289388, 0.0248571, 0.0362369, 0.00522413, 0.0130496, 0.0223266, 0.0184807, 0.0247069, 0.00393174, 0.0209758, 0.0185673, -0.00533953, -0.0576861, -0.00136327, 0.088228, 0.104526, 0.162885, 0.440797, -0.212865, -0.0405528, -0.0373798, 0.53411, 0.0759461, -0.134568, 0.143933, 0.185544, 0.313432, 0.543655, 0.524747, 0, 0, 0.285554, 0.0938175, 0.0140067, 0.16869, -0.0379049, -0.884688, -0.843638, -0.0111433, 0.00680549, -0.765245, -0.773447, -0.929118, 0.241136, 0.243169, 0.0161056, -0.00362786, -0.0238943, -0.0331173, -0.0237155, 0.00412338, -0.000916272, 0.00211043, -0.00876187, -0.0425422, 0.0231618, 0.0484422, 0.0539192, 0.0207796, 0.0275439, 0.030618, 0.0394435, 0.0106166, 0.00785114, 0.0230352, 0.00930019, 0.00907458, -0.00385633, -0.0183487, 0.017975, 0.014572, -0.00925986, -0.0134209, -0.00584805, -0.000981601, 0.0095234, 0.000741293, 0.000746297, 0.012422, 0.0186821, 0.0164952, 0.00824205, 0.0103191, 0.0227753, 0.0123428, 0.02712, 0.0209312, 0.0351064, 0.0387539, 0.00845314, 0.00490201, 0.0247286, 0.0167781, 0.0199789, 0.0180647, -0.00741062, -0.0596374, 0.00773187, 0.0455633, 0.113919, 0.479269, 0.42, -0.00933152, -0.116784, 0.00675849, 0.230054, 0.436784, -0.915557, -0.430598, -0.0562791, 0.196232, -0.0363856, 0.154724, 0, 0, -0.0968438, 0.181096, -0.547241, -0.153808, -0.657329, -0.0163168, -0.073959, -0.00150703, -0.677391, 0.0179591, -0.0876846, 0.285846, 0.173226, -0.00890754, 0.0272661, 0.00990272, -0.00823185, -0.0271108, 0.00147029, 0.014102, 0.00624024, 0.00404098, 0.00680453, 0.00191996, 0.00529038, 0.0440742, 0.0327103, 0.027005, 0.0332788, 0.0255052, 0.0168161, 0.0174478, 0.0235832, 0.0294489, 0.00813896, 0.0242016, 0.00108998, 0.0163029, 0.0319206, 0.00297364, -0.00909728, 0.00773821, -0.00431425, 0.00620089, 0.00130829, 0.0164394, 0.00849434, 0.00819515, 0.0179435, 0.019757, 0.0140886, 0.017393, 0.0367399, 0.0181791, 0.0391924, 0.026121, 0.0613554, 0.0159572, -0.0047708, 0.00307815, 0.00275157, 0.00814405, 0.015414, 0.017374, -0.0231903, -0.0619798, -0.0017419, 0.0666946, 0.144328, 0.270529, -0.0127177, -0.163692, 0.333461, 0.21393, -0.820249, -0.0634003, 0.426521, 0.159808, -0.144702, -0.0551698, -0.526596, -0.768945, 0, 0, 0.213148, 0.48736, -0.569437, -0.927744, 0.266751, -0.254585, -0.240123, -0.169883, 0.0121259, 0.0129578, 0.0356094, -0.102713, 0.283517, 0.0958008, 0.0312318, -0.00852535, -0.0156626, 0.0128949, 0.0206644, 0.0128095, 0.00308875, 0.0205631, 0.0212793, 0.0145538, 0.00720527, 0.013139, 0.0383556, 0.0168298, 0.0202389, 0.0155021, 0.0166551, 0.0177693, 0.025016, -0.00443178, 0.0143339, 0.00464025, 0.00261175, 0.0222239, 0.0118307, -0.00849822, -0.00289024, -0.00334881, -0.00839379, -0.00677251, 0.00468507, 0.00838756, -0.0115318, 0.00225866, 0.00335871, 0.0101392, 0.0241945, 0.0212308, 0.0217558, 0.00984429, 0.0368488, 0.0437519, 0.0286635, 0.0431574, 0.0167878, 0.0201938, -0.00371585, 0.0220182, 0.00366472, 0.00967318, -0.00561887, -0.050555, 0.00109159, 0.048834, 0.128418, 0.477207, -0.0244346, -0.00461017, 0.159593, 0.0103105, -0.0981908, -0.434868, -0.0726663, -0.475869, 0.132743, -0.488725, -0.0275488, 0.170171, 0, 0, -0.665591, -0.548791, -0.155653, -0.362718, 0.274657, -0.358869, 0.27101, -0.0274195, -1.01241, 0.177572, 0.312746, 0.170328, 0.569206, 0.0811808, 0.0508739, -0.0145117, 0.0172258, -0.0285097, -0.0011462, -0.0294622, -0.012438, -0.00879043, -0.0700524, -0.0474107, 0.0160687, 0.028375, 0.0248763, 0.0183554, 0.0245794, 0.0410419, 0.0229063, 0.0152291, 0.0276783, 0.0161111, 0.000125525, 0.0165175, -0.00518105, 0.0131917, 0.022362, -0.0044441, -0.0210421, -0.0175121, -0.0221705, -0.00459055, -0.00907556, -0.00243828, 0.00351176, 0.00288084, 0.00409118, -0.00401313, 0.0119452, -0.011205, -0.0153165, 0.0312393, 0.0134687, 0.0299988, 0.015703, 0.0318357, 0.0442829, 0.0170794, 0.0172913, 0.0274013, 0.0236871, 0.00813454, 0.00787659, -0.0701278, 0.0311139, 0.0820691, 0.109051, 0.511492, -0.390849, -0.0243418, -0.0447951, -0.0630402, -0.0103672, -0.632792, 0.0412574, -0.276422, -0.360918, -0.281214, -0.415826, -0.0603878, 0, 0, 0.518902, -0.00917041, 0.619085, -0.330217, -0.375147, 0.0438818, -0.244806, -0.0734778, -0.427721, 0.00830995, 0.0557591, -0.0960946, 0.234338, 0.0679191, 0.0275283, 0.0004589, -0.0211333, -0.0151344, -0.0029067, -0.0189316, -0.010333, -0.015831, -0.120008, -0.0729659, 0.0100656, 0.0290164, 0.0404255, 0.027373, 0.0398937, 0.0276715, 0.0258076, 0.0256952, 0.0241144, 0.0270005, 0.0113276, 0.0123119, 0.019841, 0.0147927, 0.00230023, -0.0313803, -0.000366566, 0.00788887, -0.00411659, -0.0033315, -0.00275881, -0.00253768, -0.00774468, -0.00490827, -0.012624, 0.0152882, 0.0139935, 0.00633588, -0.12172, -0.0054805, 0.044812, 0.0401061, 0.0602699, 0.0364276, 0.0325051, 0.03415, 0.000879029, 0.00199077, 0.0024097, 0.0237892, 0.00265999, -0.0773409, -0.0157191, 0.0603549, 0.0095209, 0.116711, -0.218494, 0.0714444, -0.0292072, 0.00166384, -0.205378, -0.0674769, -0.272531, 0.251286, 0.3325, 0.110168, 0.489984, 0.137594, 0, 0, 0.282465, -0.538209, 0.355629, -0.165941, -0.10476, 0.218079, 0.0588997, -0.608306, 0.0148293, 0.083089, -0.0878697, 0.654587, -0.234143, 0.0796892, -0.0218722, -0.00147216, 0.00678093, 0.0126253, -0.00458454, 0.00224688, 0.0110782, 0.0462114, -0.00330523, 0.0103406, 0.0114749, 0.0432847, 0.0195602, 0.0329562, 0.0208327, 0.0355254, 0.0274344, 0.0244025, 0.0297783, 0.0055709, 0.0168503, 0.00905441, 0.00689504, 0.0263081, -0.0198348, -0.116094, -0.0244719, -0.00485548, 0.00637747, -0.0158007, 0.000219259, 0.0183549, 0.0122216, 0.0120023, 0.000128761, 0.000350897, 0.015845, 0.0299388, -0.00166408, 0.0248874, 0.0384204, 0.0205082, 0.039931, 0.0332823, 0.0226812, -0.00102848, -0.0640693, -0.117577, -0.00588403, -0.00527716, -0.00662141, -0.0618872, 0.00357748, 0.0985185, 0.0904423, 0.988234, 0.694153, -0.566178, -0.280197, -0.711804, -0.029273, 0.00788846, -0.0556212, -0.105823, 0.218209, -0.137276, 0.645705, 0.769655, 0, 0, -0.146334, -0.067852, 0.122968, 0.0920886, -1.22275, -0.604528, 0.0883226, -0.0617761, 0.00302242, -0.0171526, 0.00209907, 0.152934, 0.349062, 0.0730026, 0.0523228, 0.00946134, 0.008632, 0.0056303, 0.00829549, 0.000851382, 0.00792869, 0.0108204, 0.0182463, 0.024733, 0.00673304, 0.0287155, 0.0406378, 0.0312994, 0.0413864, 0.0418182, 0.0369506, 0.0370539, 0.0177369, 0.0095003, 0.0233036, -0.000700384, 0.00304831, -0.000604246, -0.000470082, -0.0200878, -0.030102, -0.00458233, -0.0219724, -0.0161205, -0.0115416, 0.00946259, 0.0109567, 0.00133263, 0.0144626, 0.0172334, 0.0141964, 0.0435944, 0.0174563, 0.00951074, 0.0260016, 0.0408961, 0.0269595, 0.0174142, 0.0330308, 0.0132602, -0.0626115, -0.0850521, -0.00432851, 0.00342182, -0.00718088, -0.0788448, 0.0248571, 0.0769142, 0.100981, 0.128756, 0.401656, 0.32818, 0.0394028, 0.0143098, -0.468516, 0.0170011, -0.222279, 0.011424, 0.314834, 0.00530675, -0.0793681, -0.0743263, 0, 0, -0.203357, 0.125817, 0.00957074, 0.0235344, -0.115017, -0.0473813, -0.0863566, -0.0580612, -0.234228, -0.012387, 0.0049144, -0.114936, 0.339114, 0.0942062, 0.0198327, 0.0198272, -0.0211762, 0.00396855, -0.00652219, 0.000736711, -0.0042167, 0.0175907, 0.00407946, 0.0252384, 0.0225823, 0.0222216, 0.0423247, 0.0197619, 0.0311646, 0.0360531, 0.0245384, 0.0139913, 0.0277021, 0.0220406, 0.00680013, -0.00630977, -0.0116784, 0.00634272, 0.013125, -0.0243763, -0.0899977, -0.0261413, -0.0216455, -0.0149805, -0.0088644, 0.0149323, 0.0229613, 0.00862829, -0.00807939, -0.00365056, -0.00528057, 0.0249938, 0.0277041, 0.0201506, 0.0372014, 0.0196213, 0.0313157, 0.0409834, 0.0545967, 0.0113809, 0.00484547, -0.00540078, 0.0105732, 0.030557, -0.0214714, -0.0741834, 0.0222679, 0.0543104, 0.132901, 0.259305, 0.317816, -0.657305, -0.0437738, -0.0445398, -0.086219, -0.611445, -0.111475, -0.104546, -0.160757, -0.335575, 0.45673, -0.2539, 0, 0, 0.246283, -0.0719589, -0.994892, -0.0653045, -0.202739, 0.586723, 0.496493, -0.0202195, 0.0177574, -0.697526, -0.0665636, 0.10782, 0.164539, -0.00468724, 0.0609013, 0.00453745, -0.0193451, 0.0232621, -0.012537, -0.0113523, 0.00887641, 0.0123339, -0.0225764, 0.0158771, 0.021631, 0.0289703, -0.0069234, 0.0138869, 0.0365123, 0.0567751, 0.0536065, 0.040065, 0.0569518, 0.0323697, 0.038831, -0.00906053, -0.00586934, 0.0188526, 0.0443188, -0.0297209, -0.0594459, 0.00131982, -0.00289056, -0.0311063, 0.0121785, 0.00486981, 0.0104585, 0.0139205, 0.0159479, 0.0154084, 0.0242853, 0.00713038, 0.0332106, 0.0263474, 0.0414948, 0.0362934, 0.0517032, 0.0503997, 0.0293981, 0.0109158, 0.00803413, 0.0311589, 0.011869, 0.0188506, 0.000807375, -0.0863201, 0.0358411, 0.0289555, 0.0814848, -0.132811, 0.646543, 0.0457598, -0.338453, -0.000841631, -0.0208658, -0.89382, -0.0434878, -0.00844001, -0.227162, -0.491601, -0.0264096, 0.426428, 0, 0, 0.00372199, 0.286412, 0.0183136, -0.748288, -0.900502, -0.547032, -0.335702, 0.198175, -0.455579, -0.199892, 0.0343837, 0.510331, 0.403131, 0.00173834, 0.0563006, 0.00926242, -0.0376458, 0.000531692, -0.020524, 0.00411826, -0.00584741, 0.0271808, 0.0205406, 0.0516302, 0.0166334, 0.0403845, 0.0242353, 0.0138148, 0.0357146, 0.0538779, 0.0396519, 0.0217732, 0.0192475, 0.0372445, 0.033549, -0.0111482, -9.72724e-05, -0.00545523, 0.0282897, 0.00682413, -0.00387146, -0.00322451, 0.0109015, 0.0044109, 0.00170595, -0.0103847, 0.0214663, 0.0138048, 0.0270979, 0.0297033, 0.0267648, 0.00504753, 0.0334063, 0.0334338, 0.0358756, 0.0339381, 0.073705, 0.0445921, 0.0449998, 0.0265696, 0.0268274, 0.0032604, 0.0129649, -0.00482095, -0.0015237, -0.0925757, 0.00952142, 0.0382732, 0.0819697, -0.0235272, 0.508963, 0.0564022, -0.0471959, -0.084742, -0.311661, -0.0683124, -0.0250345, 0.262046, -0.227139, 0.132215, 0.121797, 0.222294, 0, 0, -0.0685417, 0.102452, -0.384681, -0.5105, 0.0625694, -0.767467, -0.490105, -0.0196521, -0.0220972, -0.0544016, -0.0925155, 0.164404, 0.292962, 0.0364602, 0.043254, 0.00677732, 0.000481197, 0.00550375, -0.0162834, -0.0218137, 0.00267434, 0.0198052, 0.0115759, 0.0345459, 0.00594962, 0.0100011, 0.0402728, 0.00499571, 0.00570409, 0.0221696, 0.0140657, 0.0199732, 0.00587501, 0.00390848, 0.0183446, -0.000766243, -0.00315815, 0.0146144, 0.00998808, 0.00706028, -0.00910795, -0.0143504, -0.0158773, 0.0208904, 0.00847494, 0.0113172, 0.00142912, 0.0135999, 0.0101645, 0.0298408, 0.0154765, 0.0283115, 0.0290653, 0.0213282, 0.0366113, 0.0433472, 0.0049613, 0.0197069, 0.00720653, 0.0269874, 0.0189334, 0.00773565, 0.0177883, -0.00125039, -0.0336874, -0.0685633, 0.0089281, 0.0702273, 0.108329, -0.00939061, 0.471226, 0.201797, -0.0651638, -0.139713, 0.0257269, -0.256268, -0.332082, -0.222797, -1.15328, 0.0753911, 0.240489, 0.046298, 0, 0, -1.01876, 0.593438, -0.0613108, 0.0399549, -0.00248154, 0.00151671, -0.0227075, -0.020484, 0.0782725, -0.0353005, -0.317624, 0.143725, 0.173164, 0.0386616, 0.000459277, -0.0069127, -0.00398356, -0.0203686, -0.0099291, -0.00940823, -0.00884635, 0.00959506, 0.00245873, 0.0171696, -0.0172876, -0.00896671, 0.0268062, 0.0247659, 0.0284506, 0.0336323, 0.00295828, 0.00683507, 0.0150019, 0.0132954, -0.000158027, 0.00462248, -0.00675035, 0.00426709, 0.0183899, 0.00604282, -0.0165633, -0.00178701, -0.00516577, -0.00854568, 0.000532439, 0.00316793, 0.000416409, 0.00931765, 0.00667184, 0.0162837, 0.0192017, 0.0269953, 0.0312393, 0.0223733, 0.0284207, 0.0270282, 0.0303332, 0.00618649, 0.00773094, 0.0234434, 0.0285912, 0.0387504, 0.0300578, 0.0187718, -0.0253459, -0.0666918, 0.0257714, 0.0632744, 0.0913376, 0.220237, 0.168844, 0.00724778, -1.04144, -0.806747, -0.0835278, -0.241282, -0.68873, -0.107098, 0.316019, 0.228518, -0.0225996, 0.00139005, 0, 0, 0.265442, -0.431818, -0.19796, 0.0899678, -0.155252, 0.00215466, -0.0127861, -0.014815, -0.00744308, -0.149457, 0.0308132, 0.257472, 0.409868, 0.0524471, 0.00284593, -0.038289, -0.0234874, -0.00842504, -0.037795, -0.00572459, -0.0044657, 0.00522617, 0.0211807, 0.0163307, 0.0322805, -0.00507713, 0.0191668, 0.0155016, 0.0219234, 0.0261445, 0.0327369, 0.0294832, 0.0279805, 0.018769, 0.00713055, 9.14024e-05, -0.00525414, 0.0174684, 0.0295901, 0.00835123, 0.0126229, 0.0109286, 0.0169899, -0.00105976, -0.00733655, 0.00550842, 0.00875289, 0.0216723, 0.025018, 0.0210606, 0.0188668, 0.0209328, 0.0473916, 0.0443143, 0.0473115, 0.0373988, 0.0505317, 0.0232437, 0.032629, 0.0138938, 0.0279738, 0.0304912, 0.0260407, 0.0078039, -0.0433301, -0.0693908, 0.0244546, 0.0438885, 0.0919842, 0.327894, 0.273599, 0.124876, -0.248923, 0.198562, -0.0649173, -0.0372624, -0.580333, -0.0546666, -0.434312, -0.142273, -0.0792224, 0.0379805, 0, 0, -0.0392965, 0.225817, 0.20675, 0.872424, -0.528384, -0.11823, -0.183989, -0.899002, -0.0839549, 0.00255322, -0.663782, -0.595813, 0.038134, 0.00281667, 0.032252, 0.00980012, -0.02638, -0.00764993, -0.0224687, -0.00277274, 0.00262855, 0.0081975, 0.020377, 0.0491066, -0.00103837, 0.00702655, 0.0216671, 0.0181806, 0.0405864, 0.0277735, 0.0233342, 0.02551, 0.0049571, 0.00612553, 0.0128414, -0.00189984, -0.00878928, 0.0102623, 0.0262409, 0.015908, 0.00127095, 0.0212317, 0.00241519, 0.0141779, -0.00150357, 0.0222476, 0.0165319, 0.0251957, 0.0263582, 0.02577, 0.0328827, 0.028204, 0.0428357, 0.0502387, 0.0477567, 0.0375045, 0.0607335, 0.0242411, 0.0379453, 0.0336771, -0.00554722, 0.00947403, -0.00164044, -0.0119889, -0.0436832, -0.071733, -0.0250161, 0.0115439, 0.126686, 0.103895, 0.483381, -0.525958, -0.26388, -0.0265104, -0.0898388, 0.0465105, -0.0623347, -0.00808122, 0.0777356, -0.655722, 0.408428, 0.11097, 0, 0, -0.0992612, -0.0577615, 0.0998013, 0.184513, 0.0839962, -0.100963, -0.334062, -0.48245, -0.682128, -0.0489821, 0.161905, -0.515182, -0.395989, 0.0858467, -0.0092678, -0.0252763, -0.0210047, -0.0165576, -0.00378644, -0.00472319, -0.0145516, 0.00882283, 0.0248209, 0.0360069, 0.0177375, 0.0127712, 0.0290282, 0.0247178, 0.014408, 0.00610567, 0.0202675, -0.00123078, 0.0169757, 0.00303772, -0.000792381, -0.0109072, -0.0073194, -0.00222963, 0.0205946, 0.0095631, -0.0108242, -0.00409214, 0.00556835, -0.00715496, 0.0135197, 0.00900683, 0.0239569, 0.0159028, 0.0215783, 0.0103107, 0.0224353, 0.0329219, 0.036147, 0.0310862, 0.0406853, 0.0419585, 0.0305494, 0.0274095, 0.0248481, 0.00995389, 0.00492556, 0.0145553, 0.0108515, -0.00138313, -0.0264144, -0.087739, -0.0254928, 0.0484336, 0.0387579, 0.166979, -0.0803578, -0.33586, -0.0763215, -0.119233, -0.138447, -0.387759, -0.0820572, 0.0993884, -0.263588, -0.495411, 0.46916, 0.116763, 0, 0, -0.245613, -0.113703, 0.0118875, 0.0296809, -0.306364, -0.652218, -1.04544, -0.104275, -0.0697636, -0.513971, 0.247091, -0.232931, -0.163896, 0.067742, -0.0052943, -0.0340804, -0.0415384, -0.0359707, -0.0392279, -0.034303, -0.0277956, 0.000509346, 0.00542796, 0.0375028, 0.0467412, 0.029736, 0.0528485, 0.0153851, 0.0219206, 0.0209196, -8.74976e-05, 0.0106978, 0.0187066, 0.012489, 0.00293893, 0.0138717, 0.0121022, 0.00119003, 0.0147978, 0.00681378, -0.0199709, 0.000741018, 0.00289361, -0.000441383, -0.00623635, -0.00568748, 0.000441086, -0.00259727, 0.0027266, 0.00276687, 0.0115583, 0.0114489, 0.0110668, 0.0306436, 0.034739, 0.0151023, -0.00806518, 0.0314449, 0.0379138, 0.0272257, 0.0173208, 0.0131249, 0.0134361, 0.00881735, -0.0235424, -0.0549539, -0.0121637, 0.0795593, 0.0311638, -0.916532, 0.219572, -0.0106419, -0.0752538, -0.0810391, -0.699387, 0.047105, -0.234077, -0.100498, -0.15074, 0.23341, 0.19322, -0.0637937, 0, 0, 0.48725, 0.392468, -0.438235, -0.317283, -0.0171185, 0.34949, -0.140288, -0.0522922, -0.0512312, -0.300687, -0.073301, 0.153624, -0.430166, 0.0641011, 0.0121458, -0.0107392, -0.0394629, -0.0434752, -0.0306192, -0.0388062, -0.0202683, -0.0131464, 0.0287843, 0.0403245, 0.0111263, 0.0417975, 0.0333501, 0.0123118, 0.0231154, 0.028805, 0.0243247, 0.0209184, 0.0209486, 0.0159613, 0.00985221, 0.0148775, 0.00584829, 0.00849495, 0.0333869, 0.0118524, -0.0042414, 0.00161843, -0.0200042, 0.0102891, 0.00341239, 0.0107654, -0.00325793, 0.0079865, 0.02274, 0.0109389, 0.0243055, -0.0528661, 0.0239813, 0.0521115, 0.0282039, 0.00616737, 0.0349059, 0.0355754, 0.0348149, 0.0152696, -0.00463567, -0.00699199, 0.021354, 0.00864328, -0.000506012, -0.0638973, 0.00274869, 0.0519887, 0.101596, 0.227448, -0.0328892, 0.0100104, -0.0289931, -0.0425224, -0.07143, 0.156121, 0.461259, 0.230923, 0.0574211, -0.422733, -0.210901, 0.305272, 0, 0, 0.698438, -0.0179089, -0.0329178, -0.147486, -0.39451, -0.905195, -0.0526279, -0.951764, -0.7772, -0.0305975, 0.00290463, -0.369503, 0.257101, 0.112471, -0.0163157, -0.0461737, -0.0333175, -0.0131869, -0.00221316, -0.0207256, -0.0115925, -0.0140795, 0.0138267, 0.0235427, 0.00958422, 0.034613, 0.0221475, 0.0207126, 0.0310965, 0.0237246, 0.0160446, 0.024808, 0.00609561, 0.00803639, -0.00236023, -0.000470391, 0.000138559, 0.0113679, 0.0273229, 0.0157627, -0.0101368, -0.0192116, -0.114096, -0.0133968, 0.010782, -0.0100384, 0.0170714, 0.0172206, 0.026701, 0.0221853, 0.0145396, 0.00173154, 0.0321881, 0.0410344, 0.0595316, 0.0453291, 0.0470307, 0.0277713, -0.0200028, -0.00750242, -0.021993, 0.00976053, 0.00774041, -0.0116371, -0.00828527, -0.078889, -0.016722, 0.0242837, 0.114925, 0.198791, -0.211402, -0.409901, 0.0567135, -0.0438111, -0.0998713, -0.916552, -0.0278018, 0.0544017, -0.000724608, 0.0258834, 0.263195, -0.00487059, 0, 0, 0.650388, 0.0502527, -0.3966, -0.227395, 0.164546, -0.0151218, -0.150059, -0.0625071, -0.910209, 0.323126, -0.0198686, -0.160591, 0.292992, 0.112411, 0.01494, 0.00156862, -0.0104466, 0.000284653, 0.00365058, -0.0148429, 0.016378, 0.0156945, 0.01143, -0.0116057, 0.0171195, 0.014455, 0.0265213, 0.028409, 0.0236109, 0.0223981, 0.00744843, 0.0124471, 0.00517267, -0.0139489, -0.00906393, -0.0085006, 0.010049, 0.0115086, 0.0143984, 0.00513182, -0.0134541, 0.00621427, -0.0141339, -0.0144157, 0.00370398, 0.0075745, -0.00247875, 0.00536246, 0.00966803, 0.0210591, 0.0133996, 0.0186089, 0.00483288, 0.0265817, 0.0385742, 0.0391408, 0.0249114, 0.0245516, -0.00630635, 0.00856558, 0.011787, 0.0196222, 0.0134523, 0.0053506, -0.0403884, -0.07246, -0.00837899, 0.0669184, 0.0598687, 0.0180422, -0.267703, -0.599739, -0.0436834, -0.252406, -0.113267, -0.115737, -0.0925753, 0.307375, -0.530153, -0.123443, 0.0390831, -0.306688, 0, 0, -0.0992588, -0.0432217, 0.0116146, -0.892421, 0.0376836, -0.0587217, 0.0104119, -0.435426, -0.59713, -0.0428942, -0.114366, 0.224791, 0.255811, 0.0789262, 0.025715, -0.0278015, -0.0212532, -0.00625368, -0.00852712, -0.0158437, -0.000784756, -0.00906278, 0.00682517, 0.00341883, 0.0164873, 0.0323086, 0.037833, 0.0397897, 0.0288629, -0.00464241, -0.0180699, 0.0151558, 0.00523114, -0.00116467, 0.00617694, -0.00636661, 0.0194085, 0.00614572, 0.00415125, -0.00734223, -0.00758562, 0.00094133, -0.00900085, -0.00259053, -0.0165463, -0.00990726, -0.00585919, -7.47582e-05, 0.00132341, 0.00389527, 0.0196897, 0.0055524, 0.0293384, 0.0104818, 0.0222517, -0.00258699, 0.000517312, 0.00729022, -0.00606104, 0.00306823, 0.000680861, 0.0187049, 0.0326464, 0.0107634, -0.0242755, -0.0724112, -0.00807652, 0.0548901, 0.0615586, 0.237624, 0.214874, -0.525542, -0.0726131, -0.0791607, -0.953991, 0.0304966, -0.161459, 0.0554141, 0.0659353, -0.626632, -0.0274396, -0.0220616, 0, 0, 0.252981, 0.119262, -0.487641, -0.00440314, 0.846196, -0.0445451, -0.603208, -0.0663854, -0.475704, -0.366212, -0.0081511, 0.636296, 0.0338984, 0.0892305, 0.0352755, -0.031756, -0.0146112, -0.000945929, -0.0119249, 0.0181721, -0.00164474, 0.0136734, 0.00427929, -0.00391912, 0.026385, 0.0417624, 0.0396396, 0.0149601, 0.0225314, 0.0153201, -0.0728673, 0.0131518, 0.0174683, 0.0192791, 0.00693971, 0.021672, -0.00481689, 0.0036607, 0.0318913, 0.01254, 5.97464e-05, 0.00500501, -0.0143941, -0.0136788, 0.00285695, 0.0066917, 0.00101396, 0.0168939, 0.0146894, 0.0162697, 0.0177935, 0.0240098, 0.0292357, 0.0420592, 0.0387172, 0.0139875, 0.0295275, 0.0253007, 0.0152955, 0.00151956, 0.03223, 0.0309112, 0.00679317, 0.0162584, -0.0250995, -0.0835707, -0.0104919, -0.0453841, 0.0496603, 0.141457, 0.244237, -0.0270638, -0.344467, -0.0741511, -0.0764805, 0.269878, -0.70601, -0.346114, -0.19298, -0.203243, -0.0909161, 0.942457, 0, 0, 0.147294, -0.265537, -0.1897, -0.511291, 0.0689881, -0.0716478, -1.11696, -0.0868843, 0.199653, -0.0305479, -0.562658, 0.148485, 0.276318, 0.0597316, 0.0503935, 0.00280234, -0.0186057, -0.00417385, -0.00663404, -0.000977621, -0.0187945, -0.00865054, -0.00960446, 0.0197299, 0.0283923, 0.0319663, 0.0363241, 0.0263971, 0.00492914, 0.0327512, 0.00897477, 0.0287664, 0.0296527, 0.0115264, -0.00308887, -0.00526423, -0.0109135, 0.00353958, 0.00028131, 0.00904816, 0.000398505, 0.0133579, -0.00160332, -0.00587989, -0.00224084, -0.00461353, 0.0108367, 0.0172536, 0.0166681, 0.00791887, 0.0264345, 0.029616, 0.0145689, 0.0296446, 0.0295071, 0.0351783, 0.0255033, 0.0163904, 0.00171528, -0.00336701, -0.00550662, -0.0104476, 0.0192379, 0.0102066, -0.0104267, -0.0528547, 0.0147369, 0.0505313, 0.114657, 0.541913, 0.161698, 0.011645, -0.0825849, -0.0139871, -0.0243065, 0.124075, -0.65216, -0.444474, -0.000166958, 0.00171508, 0.102054, 0.824594, 0, 0, 0.37049, -0.0747746, -0.536755, 0.0171254, -0.00872758, 0.0420486, 0.141932, -0.0844821, -0.0369042, -0.406891, -0.0328018, 0.114066, 0.667828, 0.0900975, 0.0427256, 0.00855759, -0.0106994, 0.0282212, -0.0107679, 0.000370212, -0.00309426, -0.00220198, 0.000678782, 0.013579, 0.0141079, 0.0198859, 0.0378348, 0.00714675, -0.00228261, 0.00699546, 0.0143833, 0.00943658, 0.00646608, -0.0074553, 0.00111332, -0.0111301, -0.0119289, -0.0119794, 0.00521208, -0.00559213, -0.00905372, 0.010236, -0.0211618, -0.00371007, -0.0153717, 0.00271828, -0.00970323, -0.00339922, 0.00596341, -0.00578064, -0.00662926, 0.0223744, 0.00716139, 0.0193817, 0.0261908, 0.0260653, 0.0236919, 0.0226802, 0.0228638, 0.00677364, -0.0278927, 0.0140775, 0.00705368, -0.00545304, -0.0309438, -0.0550137, 0.00918721, 0.0719731, 0.13495, 0.336243, -0.655015, -0.0138721, -0.0608827, -0.0556678, -0.597054, -0.0369204, -0.0280419, 0.0167134, -0.0249108, 0.0243958, 0.196598, -0.0767991, 0, 0, 0.508009, 0.204864, -0.00102597, 0.100752, -0.162777, -0.138688, 0.124386, -0.555079, -0.43677, -0.45401, 0.00678712, -0.722982, 0.332607, 0.0455621, 0.0766814, 0.030014, -0.0053889, 0.00244459, -0.00131996, 0.0163149, 0.0153865, 0.0172206, 0.00930722, 0.0184641, 0.0187241, 0.0126064, 0.0131554, 0.0112952, 0.00834291, -0.0114998, 0.00610614, -0.0149451, 0.000544372, 0.000341439, 0.00413242, -0.011546, -0.00623196, -0.00724011, -0.00512765, -0.0174467, 0.000387605, 0.0170659, -0.0194602, -0.0107724, -0.00664718, -0.00386877, -0.00962072, -0.00733194, -0.00165608, -0.00120531, -0.00554244, 0.00491247, 0.016125, 0.0156923, 0.0378768, 0.0246114, 0.00854811, 0.0209285, 0.0322715, 0.0187322, 0.0115713, -0.00758071, 0.00288845, 0.0031621, -0.021632, -0.0536008, 0.00958466, 0.0553998, 0.0516421, -0.362804, -0.359848, -0.0264263, -0.20272, -0.845671, -0.0426512, -0.52084, -0.122174, -0.306362, 0.138053, 0.367971, -0.162902, 0.373646, 0, 0, 0.302924, 0.391026, -0.00252572, -0.486831, 0.336509, -0.0184343, -0.737242, -0.076899, -0.00236726, -0.00481738, -0.459048, 0.178736, 0.0429152, 0.070081, 0.0561737, 0.00445502, 0.00192745, -0.0363317, 0.0237657, 0.000389743, -0.00226704, 0.0162825, 0.0167321, 0.0102697, 0.0055537, 0.0392202, 0.0194265, 0.0106656, 0.0168429, 0.0272514, 0.0116789, 0.0142826, 0.0265024, -0.00316317, -0.00282664, 0.001752, 0.0120221, -0.00347178, 0.0202143, 0.00380035, -0.00845979, 0.00626398, 0.000301709, -0.0114703, -0.00960216, 0.00683605, -0.0102642, 0.00193381, -0.00151242, -0.00227414, 0.00939337, 0.0216271, 0.0253337, 0.030132, 0.0255289, 0.028849, 0.0324897, 0.0148348, 0.00805963, 0.0215771, 0.0217038, 0.00826188, -0.0112783, -0.00884127, 0.00963741, -0.0698863, -0.000911436, 0.0584377, 0.0495732, 0.152062, 0.155532, -0.281291, -0.0319968, -0.039209, -0.724134, 0.202795, -0.118387, 0.519834, -0.287817, -0.210424, 0.0220311, -0.191867, 0, 0, 0.147941, -0.202518, 0.216698, 0.553142, -0.000751373, 0.187399, -0.200728, -0.115244, -0.00410032, 0.115761, -0.020317, -0.0569364, 0.325236, 0.156513, 0.0553213, -0.00183581, 0.000264857, 0.0263455, 0.0176887, 0.00296128, 0.0158467, 0.0198843, 0.01865, 0.0423801, 0.01508, 0.0365796, 0.0244824, 0.023025, 0.00885249, 0.030193, 0.0071067, 0.00250964, 0.00936573, -0.00439593, -0.0118937, -0.00987785, -0.00214435, -0.00119421, -0.00902052, 0.0126012, -0.00332406, 0.00244584, 0.000398122, -0.017264, 0.00147485, 0.0045104, 0.01209, 0.0125444, 0.00761086, 0.0235679, 0.0288435, 0.021589, 0.0245993, 0.0336267, 0.0343748, 0.0255069, 0.0514035, 0.0422775, 0.0179762, 0.00764799, 0.00210643, 0.0227221, 0.00109554, 0.00307974, -0.0145877, -0.100545, -0.0135027, 0.0763785, 0.0977904, 0.232717, -0.453753, -0.0146168, -0.300841, -0.0149778, -0.0609157, -0.417796, -0.15992, 0.31539, -0.555992, -0.0660149, 0.174614, 0.375595, 0, 0, 0.0048925, -0.0727787, 0.220221, -0.2742, -0.197246, 0.099072, 0.130148, -0.299925, 0.0496679, -0.0657984, 0.430617, 0.0772272, -0.0770005, 0.12025, 0.0285729, 0.0219299, 0.0156657, 0.00166793, 0.0201434, 0.0189229, 0.0167526, 0.014749, 0.0203046, 0.0384957, 0.00602998, 0.00531123, 0.0309888, 0.00829964, -0.0020951, 0.0170037, 0.00135018, 0.00985851, -0.00364041, -0.00533802, -0.00421934, -0.0116825, -0.00861549, -0.00165728, 0.00746814, 0.0100959, -0.0120587, 0.000958517, -0.0128713, -0.0252199, 0.000903971, 0.0106806, 0.0167861, 0.00263956, 0.00442836, 0.00588989, -0.00797339, 0.026586, 0.030576, 1.83148e-05, 0.0264816, 0.0178928, 0.0171183, 0.0382019, 0.047579, -0.00355839, 0.0212228, 0.0116234, 0.00844567, 0.00993479, -0.0121323, -0.059014, -0.0269984, 0.0615179, 0.117055, 0.585423, 0.684157, 0.104034, -0.0666031, -0.181917, -0.501975, 0.0882885, 0.000872051, -0.279825, -0.0781391, 0.0266521, -0.0373671, -0.647814, 0, 0, 0.205179, 0.51881, -0.0930054, -0.913266, 0.0177702, 0.0817123, -0.0750177, -0.0750085, -0.0591276, -0.336073, 0.0356377, -0.0207855, 0.494404, 0.260488, 0.0311779, 0.00716285, -0.0428417, -0.0134084, -0.00456318, 0.0119121, 0.00260633, 0.00800102, 0.00545888, 0.0110277, 0.0151164, 0.00408046, 0.0303743, 0.0038119, 0.0104763, 0.0141529, 0.00492205, 0.00189908, -0.00971155, -0.00391769, -0.00318654, -0.0295093, -0.00602918, 0.00200748, 0.0206489, 0.018856, -0.00319855, -0.0059589, -0.00907259, -0.00706556, -0.00538776, 0.00395761, -0.000255632, -0.0021112, -0.00680653, 0.00253119, 0.0159636, 0.00263989, 0.0166313, 0.015642, 0.046997, 0.0160729, -0.000377109, -0.00923344, 0.022694, 0.00256616, 0.0140611, 0.035498, 0.00737457, -0.00882354, -0.000860461, -0.0865804, 0.00250541, 0.0874217, 0.095041, 0.316977, -0.0837397, 0.0160659, -0.920193, -0.213996, -0.899954, 0.0726788, -0.213273, -0.1066, 0.140314, -0.0283222, 0.372705, -0.796626, 0, 0, -0.0616946, -0.0528123, -0.343999, -0.105818, 0.341667, -0.0718149, -0.774621, -0.0482305, 0.00824347, -0.013803, -0.0190731, -0.134246, -0.525725, 0.0524497, 0.0274307, -0.0634279, -0.00673913, -0.00912215, -0.0157237, 0.00423373, -0.00561372, 0.00593607, -0.0179501, -0.0205505, 0.00743718, 0.0290636, 0.0141892, -0.0268146, 0.00949452, 0.0305913, 0.0254939, 0.0177529, 0.0261567, 0.0131034, 0.00332722, -0.0642759, 0.0055598, 0.00888913, 0.0219003, 0.022078, 0.000949581, 0.00899298, -0.0122688, -0.013102, -0.00549916, 0.00672696, 0.00742359, 0.00766048, -0.0070827, -0.0232299, -0.0139869, -0.0101353, 0.00117657, 0.0190632, 0.0477351, 0.0180093, 0.0200713, -0.146685, -0.021362, 0.0260758, 0.0233105, 0.0236924, 0.00061149, -0.0141818, -0.0167007, -0.060627, 0.0229977, 0.0180896, 0.128399, 0.579606, 0.00977599, -0.0095084, -0.0308555, -0.448108, -0.310554, 0.0325797, 0.179778, 0.300556, 0.0368795, -0.58094, -0.209208, 0.154696, 0, 0, 0.255315, -0.0406512, -0.230167, -0.332973, -0.0112345, -0.870278, 0.0019795, -0.725944, -0.0334865, -0.0282818, -0.16167, 0.152365, 0.281829, 0.11045, 0.0282845, 0.00339088, 0.0182179, -0.0121355, 0.00279382, 0.000947757, -0.0192798, 0.00386849, -0.00598062, -0.0201355, 0.00432094, 0.0156285, -0.0141063, -0.145358, -0.00524363, 0.00756815, 0.00382831, 0.00971231, 0.0208189, 0.00684301, -0.00317452, -0.0148968, 0.0114113, 0.00123251, 0.0155732, 0.0135194, -0.00606739, -0.00128068, -0.00182919, -0.00637761, 0.00528734, 0.00100327, 0.0131439, 0.0158229, -0.000574546, 0.0168126, 0.00210108, -0.145558, -0.0223899, 0.00388094, 0.0314101, 0.0214451, 0.0311446, -0.0147641, -0.0275539, 0.00340645, -0.0118624, 0.018025, -0.0112756, -0.0109004, -0.017222, -0.0434719, 0.0186645, 0.0149186, 0.156765, 0.229706, -0.101683, -0.149892, -0.488184, -0.0438473, -0.825203, 0.0389099, -0.708294, 0.0571473, 0.0691868, -0.0144293, 0.540723, 0.23486, 0, 0, -0.0928547, 0.50587, 0.0707275, -0.169811, -0.350696, 0.0351496, -0.211081, -0.29953, -0.129368, 0.0260164, -0.037699, -0.768449, 0.272574, 0.0805908, 0.0594952, -0.0224014, 0.00346892, -0.00913754, 0.00687869, -0.00577127, -0.000543885, 0.00709044, -0.0189893, 0.0129917, -0.00449072, -0.00469534, 0.00918129, -0.0296887, -0.014911, -0.0134627, -0.0289903, -0.00652499, -0.00547025, 0.00302434, -0.0127285, -0.0054028, -0.0179597, -0.0113184, 0.00138318, -0.00791181, -0.0323312, -0.00282967, -0.0197283, -0.00689776, -0.0100884, -0.00223654, -0.00919393, 0.00480874, -0.0163191, 0.000827612, -0.0179603, -0.0171037, 0.00867333, 0.0142669, 0.0256921, 0.0281514, -0.00368459, 0.0366529, -0.014798, -0.00953416, 0.00784489, 0.0138663, -0.015232, 0.0125552, -0.0052914, -0.0924471, 0.023098, 0.0582814, 0.16249, 0.165202, 0.421409, -0.0220976, -0.0870507, -0.0733237, -0.795946, -1.08217, -0.164853, -0.611269, -0.26316, -0.304925, -0.180882, -0.0677237, 0, 0, -0.910419, -0.117218, 0.0474061, -0.29463, -0.0191715, -0.0996418, 0.0569543, -0.175528, 0.0715727, -0.377865, -0.147969, 0.482716, 0.66595, 0.0518676, 0.00902027, 0.0107757, -0.0245539, -0.000296237, -0.0156798, -0.00611591, -0.0085915, -0.0140093, -0.00868501, 0.0224025, 0.0049486, -0.00764806, 0.0183608, 0.0149972, 0.0021157, -0.0323493, -0.142581, -0.0102873, -0.0111501, -0.0113885, -0.0149581, -0.0117386, -0.0291949, -0.0238726, -0.00428805, -0.0128977, -0.0322522, -0.0143668, -0.0192972, -0.00446865, -0.0226752, -0.00447593, -0.012509, -0.00683194, -0.0200869, -0.0164418, -0.0171585, 0.00434479, 0.011083, -0.0213636, -0.00536442, 0.0120941, 0.00448111, 0.00510821, 0.0449189, -0.00878883, 0.0170015, 0.0109607, -0.00126038, 0.0164789, 0.00415484, -0.0918901, 0.00902795, 0.0922005, 0.138101, 0.239967, -0.221346, 0.388596, -0.7431, -0.0349885, -0.0938241, 0.0380496, 0.0367916, -0.229838, -0.0758836, -0.0336475, 0.0777599, 0.0485188, 0, 0, 0.213743, 0.00595661, 0.144031, -0.390639, -0.174438, -0.06933, -0.264102, -0.0901481, -0.994106, -0.232529, 0.0159191, 0.453388, 0.192463, 0.211309, 0.00880768, -0.00661944, -0.0230726, -0.00609252, -0.0221554, -0.0236559, -0.0158812, -0.00865571, -0.00562447, -0.0151174, -0.0404961, 0.0488262, 0.0126605, -0.00213763, -0.029521, -0.0232662, -0.0181078, -0.0124374, -0.0110903, -0.0359884, -0.00651208, -0.0104924, -0.0279929, -0.0185733, 0.00261495, -0.00480519, -0.0211443, 2.51862e-05, -0.0150982, -0.0132034, -0.00137705, -0.00410795, 0.00617117, 0.0131772, -0.00955002, -0.00693258, -0.00935779, 0.00662442, 0.00309521, -0.00201117, 0.0100092, 0.00828049, 0.0430313, 0.0276785, 0.0280697, -0.0104084, 0.0134777, 0.00640698, -0.00900979, 0.023296, -0.000282168, -0.0705523, 0.00779331, 0.0591927, 0.199269, 0.131718, -0.92801, 0.030978, -0.0165735, -0.0206176, -0.056327, -0.0955485, -0.509421, 0.251624, -0.0792704, 0.504801, 0.562525, 0.659468, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  //jetasymmetry_norm
  // std::vector<float> TH2_contents = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.507465, 0.0796012, 0.183024, 0.212938, 0.730266, 0.131144, 0.232972, 0.231465, -0.0574707, 0.133946, -0.844814, -0.0787219, -0.0152037, -0.0163551, 0.027989, 0.0308469, 0.0019693, -0.0149127, 0.0013724, -0.0195291, -0.0107795, -0.00870356, -0.0143934, 0.000882332, -0.136103, 0.00433799, -0.00352195, -0.0217966, -0.00320519, -0.0108043, -0.00943697, -0.00156682, -0.0122157, -0.0118559, 0.0103557, -0.00153315, -0.000658004, 0.00636435, -0.00746871, -0.00851731, -0.00256419, -0.0101753, 0.00703516, -0.00249335, 0.0021446, 0.000466199, -0.00481683, 0.00902269, 0.00208355, -0.00196148, 0.00127837, -0.0105149, -0.0234657, -0.100749, -0.028911, -0.0148835, -0.0168364, -0.000106435, -0.0111054, -0.0323909, -0.0118812, 0.0119817, 0.00458793, -0.0136839, -0.0176551, -0.0317254, -0.00311071, 0.0082967, -0.0537145, -0.940951, 0.173691, 0.0271618, 0.118512, 0.110611, -1.12147, -0.700733, -0.231008, -0.848315, 0.220186, 0.458185, -0.265471, -0.360494, 0, 0, 0.404848, -0.0104684, 0.684818, 0.23301, 0.150628, -0.0574562, 0.231239, 0.217279, 0.13974, -0.717883, 0.0779257, 0.109014, 0.230435, 0.00597394, 0.0151825, -0.00798514, -0.00523494, -0.00859097, -0.000101764, -0.00119164, -0.00539382, 0.0156037, 0.0209037, 0.0115576, -0.0263665, -0.0181522, -0.00317165, -0.0162885, -0.0149336, -0.024433, -0.0134773, -0.00662291, 0.0129602, -0.00818216, -0.000674648, 0.00541831, -0.00750815, -0.00468904, -0.00638552, -0.00747796, 0.00800984, 0.00706084, 0.013477, 0.00655266, -0.0162287, -0.00729501, 0.0015935, -0.00481111, -0.00979848, 0.00530281, -0.00681717, -0.020569, -0.00352835, -0.0429322, -0.0133893, -0.0124006, -0.0175662, -0.00868263, 0.000598477, -0.0146299, 0.00167058, 0.0192479, 0.005642, -0.00253209, -0.00339242, -0.0207422, -0.00348661, -0.0230884, -0.17085, -0.203951, -0.367511, -0.0109471, 0.0611342, 0.0597771, -0.638761, 0.405634, 0.127306, 0.477882, 0.0722736, 0.169819, -0.164424, -0.144408, 0, 0, -0.193409, 0.0704684, -0.905197, 0.0316707, -0.268832, 0.0887115, 0.523357, -0.0403159, 0.167245, 0.103278, 0.37897, 0.0876719, 0.0115831, 0.058714, 0.00258173, -0.0243898, 0.00191969, -0.011517, -0.0374223, -0.0169825, -0.00612449, -0.000925643, -0.00575197, -0.0337002, -0.00419126, -0.010567, 0.00601481, -0.0171573, -0.00961165, -0.0015392, -0.00126578, -0.00669872, 0.00470078, -0.00334335, -0.00989288, -0.00277547, -0.013744, -0.00626298, 0.00149204, 0.00993829, -0.00414957, -0.0238606, -0.0109437, -0.00469799, -0.0146127, -0.011803, -0.0103703, -0.00296007, 0.0115326, -0.00746594, -0.00389898, -0.0122945, -0.0138343, -0.0211843, -0.0127381, -0.0143903, -0.0179773, -0.00125161, 0.00654168, -0.00795082, -0.0143928, -0.00428543, 0.0213084, 0.00229586, 0.0137565, 0.0376298, 0.00897149, 0.0461807, -0.00340644, 0.0946458, -0.396783, 0.253587, 0.0853267, 0.30113, -0.515035, 0.487289, 0.227264, -0.0751977, -0.315966, 0.152529, -0.0480716, -0.138713, 0, 0, -0.2688, -0.654623, -0.23011, -0.0726812, 0.0934335, -0.368647, 0.232625, -0.320506, 0.206524, 0.120227, 0.0893648, -0.10617, -0.177776, 0.00797075, -0.0442123, -0.00599339, 0.00603135, -0.0176438, -0.0213556, -0.0132834, -0.0108779, -0.00391507, -0.0185094, -0.0414886, -0.0160773, 0.00678968, -0.00173394, 0.00165083, -0.0122033, -0.0108461, -0.00135812, 0.00193834, 0.0161502, 0.00706841, -0.00541672, 0.00226648, 0.000662838, -0.00387916, 0.00623334, 0.00895252, 0.0103183, 0.0128436, -0.00211187, -0.0047729, -0.00919391, -0.0182098, 0.00910894, -0.00247426, -0.000151817, -0.00980482, -0.017452, -0.0148396, -0.00732512, -0.00419328, -0.0230181, -0.00817394, -0.000553289, -0.0207177, -0.0216196, -0.039174, -0.000592323, -0.0167229, -0.0035481, 0.00287897, -0.0220521, 0.0157952, -0.0372135, -0.0218068, -0.00484392, -0.722917, -0.438572, 0.12319, 0.119938, -0.100626, 0.120016, -0.747964, 0.666395, -0.255783, 0.211562, -0.291922, -0.274441, -0.0240356, 0, 0, 0.105532, 0.388309, 0.0527036, -0.00572825, 0.372368, 0.0310079, -0.142079, -0.683721, -0.0129152, 0.131062, 0.0194483, -0.638714, 0.108353, 0.0217421, -0.00133042, 0.000229962, -0.00493446, 0.00862471, -0.00108092, -0.0162092, 0.000923624, -0.00421173, 0.0121794, -0.0231656, -0.0242077, -0.00506742, -0.0110009, -0.0013527, -5.16256e-05, -0.022936, -0.00201295, -0.00434465, 0.00145365, -0.0063126, -0.0142526, 0.000960517, -0.00537796, -0.0142161, -0.00150329, 0.00709749, 0.00366862, -0.000641899, -0.00349017, -0.00438784, -0.0114969, -0.0192256, 0.00674492, -0.00603516, -0.00824241, -0.0111463, -0.0241455, 0.0093198, -0.0130724, -0.000724664, -0.00479745, -0.0178786, 0.00249408, -0.0110747, -0.0596441, -0.000994173, -0.0053209, -0.00144051, -0.00937248, -0.0127653, 0.00279908, 0.0186407, -0.0102118, -0.00757754, -0.0642962, -0.152892, -0.0458228, 0.605385, 0.402215, 0.591594, 0.00759949, 0.0294184, 0.311682, 0.0249182, 0.255878, -0.402368, -0.240697, 0.465084, 0, 0, 0.456733, 0.371257, 0.181129, -0.348751, -0.492062, 0.279185, -0.310194, 0.245196, -0.00971937, -0.38973, -0.664234, 0.313976, 0.210575, 0.0463579, -0.00230211, -0.00815445, -0.00401876, -1.66109e-05, -0.00542834, -0.0108393, 0.00363981, 0.00618873, 0.00420081, 0.0191333, 8.70024e-05, -0.0489405, -0.0673538, -0.0136006, -0.0217874, -0.0341552, -0.00656403, -0.0118893, -0.0174082, -0.0145179, -0.0204429, -0.012476, -0.013048, -0.00779579, -0.0151059, -0.0217071, 0.00190326, -0.005365, -0.00145606, -0.00809699, -0.00720109, -0.0250863, -0.0372309, -0.0104194, -0.0157014, -0.0117378, -0.00677177, -0.0212321, -0.00845211, -0.0225443, -0.0168747, -0.00845916, -0.0170339, -0.0192602, 0.00722784, -0.0133231, -0.0105022, 0.00269707, 0.00522469, -0.025989, -0.00211084, 0.010628, -0.0179059, -0.0120644, 0.00980931, 0.325011, -0.0939149, 0.400783, -0.14431, -0.00894576, 0.103275, -0.190297, 0.203432, 0.281942, 0.0402529, 0.735564, -0.151162, -0.0750299, 0, 0, 0.0901673, 0.0726722, 0.19501, 0.172631, 0.50351, -0.205109, -0.00816477, 0.165634, -1.10114, 0.103787, 0.0917308, 0.182363, 0.169248, 0.0340761, -0.0117205, 0.0232691, 0.0132761, -0.0757986, -0.058835, -0.0180756, 0.00251094, -0.00613381, 0.00941034, -0.000863279, -0.00181801, -0.00504352, -0.00887719, -0.0167313, -0.00507549, -0.0119106, -0.0042828, -0.0131565, -0.00536848, -0.0180556, -0.00629445, -0.0023061, 0.00230139, 0.0015366, -0.0175102, -0.118522, -0.0244602, -0.00165338, -0.00783065, -0.0200815, -0.0114842, -0.00795278, -0.0241526, -0.0266009, -0.0187046, -0.0238406, -0.0292679, -0.0169372, -0.0293133, -0.0259649, 0.00302882, 0.0101774, -0.0413945, -0.072486, -0.0106455, -0.00211014, 0.00743852, 0.0229691, 0.016491, 0.0252052, -0.00356182, -0.0120923, 0.00221182, -0.0288034, 0.00386962, 0.14774, -0.22319, 0.0231582, 0.0921232, 0.0284628, -0.37662, 0.0747781, -0.777846, -0.0824901, 0.0339186, -0.000999132, 0.328868, -1.07322, 0, 0, 0.248837, -0.898874, -0.113178, -0.288112, -0.352743, 0.118058, -0.250011, 0.116239, -0.446066, 0.119006, 0.339942, 0.104172, 0.125118, 0.0134225, 0.0514595, 0.0216726, 0.0181226, -0.217771, -0.16766, -0.00426379, -0.00678534, -0.000138611, 0.00758755, 0.00188911, 0.00166715, -0.00561095, 0.00064784, 0.00839696, 0.0110279, 0.00276051, 0.00598958, -0.00159488, -0.0116647, -0.0247525, -0.00219241, 0.00204039, 0.0183151, 0.00926324, 0.00805677, 0.00737878, 0.00542505, 0.00415139, -0.00246142, -0.00248533, -0.0178295, -0.0147452, 6.72561e-05, -0.00385706, -0.00471262, -0.00613036, -0.00556882, 0.00159633, -0.013359, -0.00175565, -0.0114332, -0.0128465, 0.000846742, -0.0089984, 0.028267, -0.00227565, -0.00723246, -0.00549154, 0.0224517, 0.00741744, -0.00131762, -0.00379517, -0.00574153, 0.0310587, -0.0450525, 0.00102435, -0.0800447, 0.115267, 0.12083, 0.0712979, -0.023491, 0.525223, 0.191792, 0.259961, -0.79798, -0.136291, -0.174271, 0.269844, 0, 0, 0.584099, -0.0893128, 0.0170604, 0.286254, 0.569247, 0.146777, 0.119617, -0.491715, 0.189969, -0.448681, -0.145837, 0.165747, -0.112999, 0.0182769, 0.0153422, 0.0417382, 0.00120913, -0.0832189, -0.0888126, -0.00786404, 0.0258338, -0.00384909, -0.0124028, -0.0138164, 0.00411737, 0.0102018, 0.000700018, 0.000326264, 0.00203166, 0.014328, 0.00686572, -0.00300505, -0.0149385, -0.00606302, 0.00122414, 0.00533491, -0.00201323, 0.0145133, -0.00873432, 0.00464461, -0.00448346, 0.00701487, 0.000911633, -0.0247622, -0.0126293, -0.000462322, -0.00890812, -0.0122403, -0.0104576, -0.0101473, -0.000680515, -0.00863839, -0.00812902, -0.0246704, -0.00144535, -0.00887119, -0.00393555, 0.010759, 0.00728176, 0.00744885, 0.016402, 0.0314985, -0.0011376, -0.00446241, 0.00713654, 0.031985, -0.00153117, 0.0375063, -0.0541968, 0.0387885, 0.590376, 0.0611124, 0.056837, -0.0854827, 0.146747, 0.0615332, -0.0904183, 0.652132, 0.0905602, 0.120625, -0.253781, -0.242811, 0, 0, -0.423029, -0.398683, 0.182991, 0.32952, -0.01693, 0.778755, 0.27361, 0.148696, 0.16082, 0.0934281, 0.0631339, 0.0747339, 0.0678742, -0.0187473, 0.0154728, 0.0106189, 0.00830067, -0.0010942, -0.00857716, 0.00345042, -0.0281243, -0.0251142, 0.0212891, 0.00939189, 0.00577745, -0.00450644, 0.000662273, -0.0140219, -0.0181085, -0.00774129, -0.00315378, -0.0136244, -0.0380917, -0.0174782, -0.00330299, -0.00526791, -0.00610472, 0.00181315, -0.0118601, 0.00123035, -0.00609716, -0.0151161, -0.00962135, -0.0216655, -0.0135237, -0.0175273, -0.0282853, -0.0201313, -0.0251359, -0.031572, -0.0114983, -0.00450421, -0.0182685, -0.020007, -0.0134262, -0.00362053, -0.0120905, -0.00618415, -0.00199895, 0.00228962, -0.0072358, 0.0104823, 0.00394639, 0.00398799, 0.0198361, 0.0120439, -0.00746682, 0.0106609, -0.0107519, 0.18979, -0.385552, 0.00486197, 0.0955249, -0.240778, 0.37038, 0.530576, -0.0902613, -0.505238, -0.813943, 0.0672219, -0.153569, -0.155037, 0, 0, 0.452523, -0.00458618, 0.0152003, -0.0153162, -0.0796731, 0.349414, 0.307178, 0.176006, 0.171279, 0.127342, -0.625605, -0.264915, 0.0625971, 0.0168788, 0.00946198, 0.00349269, 0.0161956, -0.0391673, -0.0135792, 0.000959731, -0.0227746, -0.00635603, 0.000199195, -0.00708142, 0.0190234, -0.0290701, -0.00471766, -0.0152848, -0.00400813, -0.00512425, -0.00942244, -0.0266074, -0.0995338, -0.0474388, -0.0331773, -0.0108453, -0.00732286, -0.00302823, -0.00859253, -0.00307311, -0.00939811, -0.0256612, -0.00310166, -0.0166486, -0.00876744, -0.0281526, -0.0297663, -0.0325851, -0.0250883, -0.0368855, -0.024156, -0.00886939, -0.0138615, -0.0251246, -0.0153763, -0.0145473, -0.025285, -0.0156286, 0.0191443, 0.0103784, -0.0139047, 0.00588848, 0.0136644, -0.00990056, -0.00339406, 0.0244269, 0.0024923, -0.0184435, -0.0285005, -0.120123, 0.296093, -0.199939, -0.217366, 0.0859074, 0.150356, 0.0809544, 0.119162, 0.144831, 0.00439509, -0.106714, 0.0101154, -0.0372904, 0, 0, -0.0830443, 0.69578, 0.406326, -0.0338487, -1.11473, 0.0312861, 0.0605072, -0.671109, 0.202698, 0.148029, 0.12789, -0.0802132, 0.0556198, 0.0123526, 0.0336623, -0.0327231, -0.00485251, 0.00625317, 0.00152891, 8.52052e-05, -0.00758534, -0.0167253, -0.00863777, -0.0198229, 0.000216189, 0.0102661, -0.00257356, -0.0181715, -0.00408095, -0.00244318, 0.0117464, -0.0127218, -0.0370041, 0.001999, -0.0120475, -0.00391904, -0.00200466, -0.00136742, 0.00455679, 0.0120843, 0.0056123, -0.00811707, -0.0072579, -0.0127966, 0.00268827, -0.00778032, -0.00548835, -0.0172081, -0.0261136, -0.0128478, -0.0216188, -0.00489259, -0.00468361, -0.013713, -0.00803653, -0.0177185, 0.00221131, -0.00804633, -0.0355752, -0.0476015, -0.0127049, -0.0794839, -0.0215982, -0.00917764, -0.0037352, 0.00189165, -0.0130382, -0.0234516, -0.00674752, -0.0252974, 0.0598054, 0.140011, -0.579897, -0.412036, 0.997212, -0.493679, 0.239045, 0.00983235, -0.13999, 0.000266088, -0.109387, 0.310899, 0, 0, 0.112511, 0.972695, 0.12371, 0.01389, -0.380999, -0.729759, 0.170332, 0.240079, 0.200682, 0.437696, 0.115552, 0.346171, 0.422237, 0.0540096, 0.00373391, 0.00192172, 0.00209588, -0.00194842, -0.0139833, 0.00874938, 0.0192085, 0.00719754, -0.00667794, 0.0130758, -0.0123364, -0.00754074, -0.0259373, -0.0319638, -0.00329259, 0.00486902, 0.00885602, -0.00228109, -0.0120428, -0.0162009, -0.00183495, -0.0138153, 0.00479197, -0.00988451, -0.0042201, -0.00586009, -0.00373122, -0.00522287, -0.00551452, -0.00127634, -0.00757426, 0.00296333, -0.0102404, -0.0134089, -0.000949749, -0.0141249, -0.00613948, 0.00510292, -0.0118895, -0.0218675, -0.0166773, 0.00219765, -0.00225429, -0.000738212, -0.0296897, -0.0405418, -0.0324582, -0.199749, -0.0794256, 0.0122776, 0.00863389, -0.0103936, -0.0232436, -0.0282257, -0.057769, 0.0410363, -0.137771, 0.529227, -0.36285, 0.137257, 0.184641, -0.996701, -0.311512, 0.390202, 0.0228104, 0.019066, -0.0050971, 0.0546847, 0, 0, -0.456187, -0.156135, -0.840614, -0.256079, -0.315965, -0.470229, 0.0836215, -0.873065, 0.208779, 0.103943, 0.0649832, 0.309658, 0.065571, -0.0767791, 0.0161487, 0.0113887, -0.026364, -0.000928436, 0.0165488, -0.00704153, 0.00702222, 0.0105503, -0.00863292, -0.018571, -0.0228166, -0.0570282, -0.0389174, -0.0238494, -0.050862, -0.0316999, -0.0120556, -0.0217039, -0.0203752, -0.0210811, -0.0214721, -0.0434535, -0.0254311, -0.00474111, -0.0129031, -0.00145137, -0.0121821, -0.0111736, -0.02201, -0.017269, -0.0156415, -0.0189303, -0.00782446, -0.0176109, -0.012501, -0.0367858, -0.0302841, -0.0330414, -0.0184113, -0.0155653, -0.0191292, -0.0298249, -0.0277198, -0.0306485, -0.0291278, -0.0118204, -0.00280048, -0.0315281, -0.021078, -0.0120957, 0.00271056, -0.0101209, -0.0144147, -0.0126156, -0.0540657, 0.27779, -0.0924864, -0.402997, 0.0953008, -0.109196, 0.139463, -0.556564, -0.0462593, -0.460418, -0.150012, 0.262676, -0.0755589, -0.0286506, 0, 0, 0.141723, -0.464516, 0.217748, 0.462801, 0.804173, 0.135939, 0.388447, -0.600569, 0.234436, -1.12881, -0.257329, 0.143756, -0.633982, -0.0930873, -0.0370146, -0.00551938, -0.0281607, -0.0106104, -0.0281944, -0.0161915, -0.0328599, -0.0152007, -0.0105469, -0.0237199, -0.0244804, -0.0394805, -0.0287565, -0.0160339, -0.0339419, -0.0316736, -0.0137835, -0.0180637, -0.032522, -0.0289584, -0.0365962, -0.0133692, -0.00542766, -0.00836985, -0.0311413, -0.0141715, -0.0120714, -0.020151, -0.00818119, 0.00177177, -0.0144357, -0.00534632, -0.0174252, -0.0221479, -0.0210638, -0.0213739, -0.0102684, -0.0191013, -0.0146387, -0.0128986, -0.0213729, -0.0348688, -0.0545937, 0.0026122, -0.0212884, 0.0161595, 0.00358823, 0.0161107, 0.0183385, -0.00211863, 0.0258166, 0.00776422, -0.0224813, -0.022392, -0.0228126, -0.48118, -0.560569, 0.0334464, 0.0961975, 0.0874295, -0.996084, -0.651801, -1.00416, -0.336764, 0.180023, 0.0936661, -0.0248707, -0.163696, 0, 0, 0.158773, -0.961879, -0.0390384, 0.122011, 0.268303, 0.0882348, 0.0255704, -0.846111, -0.34745, 0.110535, -0.895851, 0.0572245, -0.332221, -0.0296099, -0.0517357, 0.000909608, -0.00989064, -0.0105084, -0.0159132, -0.0246403, -0.0113524, -0.0139094, -0.0181662, -0.026228, -0.0135233, -0.00785019, -0.00748394, -0.0098986, -0.035096, -0.0311628, -0.0434805, -0.0339281, -0.0140187, -0.0250465, -0.0287851, -0.0159931, -0.0119046, -0.00744726, -0.0216514, -0.00854652, 0.00681203, 9.61479e-05, 0.010665, -0.0035811, 0.00315475, -0.00034733, -0.0143868, -0.00728838, -0.0134464, -0.0115685, -0.0103948, 0.0145376, 0.00587499, 0.00445154, -0.00820417, -0.0026358, -0.00659971, 0.00286424, -0.0115629, 0.00855968, 0.00521098, 0.0057157, 0.017182, 0.00529393, 0.0115685, 0.00603839, 0.000267604, -0.0413702, -0.0556007, 0.0724915, 0.0268458, 0.0745947, 0.0952326, 0.0857576, 0.186836, 0.162012, -0.420824, 0.0684674, 0.0161761, -0.741894, 0.622015, -0.36239, 0, 0, -0.146219, -0.166172, -0.354392, 0.303608, 0.189742, 0.0466855, 0.180402, 0.202292, 0.176355, 0.0970095, 0.0555272, -0.111476, -0.311692, 0.00161617, -0.00512263, -0.0180595, 0.00860495, -0.0300931, -0.0318901, -0.0218078, 0.000550667, 0.00891109, 0.0105518, -0.00284521, -0.00954916, -0.0171936, -0.0174284, 0.00228638, -0.0282348, -0.0445389, -0.0543819, -0.121257, -0.0225991, -0.0199847, -0.0309234, -0.0280404, -0.0181595, -0.0194736, -0.0193032, -0.0133765, -0.0136791, 0.00969771, -0.00995556, 0.0185805, -0.0146668, -0.0050731, -0.00600675, -0.00370777, -0.00499578, -0.00643066, -0.0015086, 0.00909296, 0.0079026, -0.00842592, -0.00219962, -0.0108772, -0.00205748, 0.0256072, 0.0188153, -9.25175e-05, 0.00702072, 0.000367276, 0.0059808, 0.00174413, 0.0151386, -0.0178747, 0.000287906, -0.0112974, 0.0129834, -0.100854, 0.126291, 0.0131759, -0.17108, 0.0932396, 0.273501, -0.107578, 0.247053, -0.00510846, 0.289184, 0.327354, -0.0366151, 0.507169, 0, 0, 0.0803317, -0.00696958, -0.174965, -0.476357, -0.214226, 0.307341, 0.200282, 0.787724, 0.200016, 0.12031, 0.0795716, -0.703017, 0.356421, -0.00240331, -0.0375966, -0.0506823, -0.0396719, -0.0204859, 0.0142038, -0.0042861, -0.0162062, 0.000436767, 0.0307259, -0.00121523, -0.0365486, -0.0403762, -0.0327446, -0.0161005, -0.047031, -0.0455954, -0.0528378, -0.0507513, -0.0294126, -0.0374937, -0.0348449, -0.0109731, -0.0360245, -0.0307376, -0.0227681, -0.0148103, -0.00129274, -0.0166163, -0.00883225, -0.0157151, -0.0164311, -0.00621051, -0.013413, -0.00052852, 0.00956173, -0.00362234, -0.00882557, 0.00165093, 0.00171899, 0.00886903, -0.00916197, 0.00130197, -0.0129174, 0.013793, -0.00759866, -0.00911192, 0.001167, -0.0126347, -0.00122097, -0.0115093, -0.0100729, -0.00932753, 0.00518416, 0.000670483, 0.0107682, 0.084265, -0.634217, -0.438818, 0.245696, 0.110569, 0.633186, 0.135934, -0.829022, 0.513644, -0.289247, 0.0954459, 0.281259, -0.168951, 0, 0, 0.00514551, 0.132503, 0.659284, 0.171604, -0.136055, -0.00108891, 0.0818014, -0.620278, -0.406899, 0.347107, -0.255872, -0.456446, 0.052281, -0.0580267, 0.0124119, -0.0967266, -0.15061, -0.0324344, -0.00378576, -0.00197461, -0.0133632, -0.0332094, 0.0148319, 0.0193165, -0.0229573, -0.0416496, -0.0140691, -0.0109171, -0.0253488, -0.0284546, -0.0289514, -0.030758, -0.0100773, -0.0167162, -0.0173844, -0.00210656, -0.0236542, -0.0195648, -0.00405922, -0.00923267, -0.0047556, 0.00740323, -0.0295933, -0.0222147, 0.00253694, -0.00651687, -0.019383, -0.016825, -0.0124285, -0.00802827, 0.00292087, -0.00578489, -0.00872968, -0.00435918, -0.000906928, -0.00292248, 0.000897713, 0.0101938, -0.116783, -0.0826465, 0.0116447, -0.00507013, 0.00605234, -0.0157752, -0.0143706, -0.0101288, -0.0290676, -0.0272558, 0.020096, 0.113961, -0.00417635, -0.582818, 0.0685556, -0.624735, 0.271504, -0.333928, -0.0809417, 0.0914115, -0.067816, -0.419235, -0.167667, -0.198642, 0, 0, 0.370917, 0.0785737, 0.836152, -0.64647, -0.0109127, 0.0930428, 0.0732249, 0.264996, 0.273525, -0.391337, 0.338851, 0.337389, 0.109534, -0.0497563, -0.0415624, -0.0800899, -0.310753, -0.028612, -0.0106946, -0.0155469, -0.0234962, -0.000398707, 0.0142278, 0.02106, 0.00276023, -0.00837157, -0.0128841, -0.00558828, -0.0146592, -0.0114248, -0.0101688, -0.0122567, -0.00335733, -0.0139159, 0.00414563, -0.00354778, -0.00350952, -0.0123862, -0.00596071, -0.00904975, 0.017458, -0.00725792, 0.00670205, -0.0153595, -0.00467792, 0.0110967, 0.00419651, -0.000811158, -0.000342912, -0.00836901, 0.0128247, 0.0197976, 0.0123326, -0.00999616, -0.00306644, 0.011134, 0.0309249, 0.00203927, -0.027495, -0.0153754, -0.00104086, -0.0125684, -0.0624704, -0.037823, -0.0285824, 0.0153869, -0.0188078, -0.0492256, 0.0304512, -0.146889, 0.392693, 0.161104, -0.182165, 0.0805022, -0.918224, 0.0902025, 0.119236, 0.0628835, -0.375976, 0.0251463, 0.0537184, 0.285198, 0, 0, -0.193305, 0.250744, 0.182265, 0.459792, -0.126849, 0.289716, -0.68677, 0.240833, -0.799374, 0.145542, 0.312691, -0.801933, 0.272632, -0.0319654, -0.0150972, -0.0439594, -0.113475, 0.00697547, 0.00418733, 0.0128178, 0.0124139, -0.00298136, 0.0182205, 0.00138163, -0.00171794, -0.00171108, -0.0245267, -0.0200505, -0.0271908, -0.0127301, -0.017538, -0.0143791, -0.00705473, -0.0053046, -0.0110087, -0.00794837, -0.0177693, -0.00674731, -0.00848869, -0.00334325, 0.00669846, 0.000242402, 0.00683234, 0.00637166, -0.00251393, -0.00217307, 0.0123043, -0.0172421, -0.0017401, 0.0161409, 0.014109, 0.00921768, 0.0162485, 0.00412842, -0.00129061, -4.09592e-06, 0.00380778, -0.0287031, -0.00712416, 0.0198351, 0.0164041, -0.0325345, -0.190726, -0.0615292, -0.0340709, -0.034332, -0.000388609, 0.0318014, -0.0342544, -0.474295, 0.432683, -0.0927452, -0.463627, -0.215474, -0.0171606, -1.10262, 0.318428, -0.290793, 0.0993222, 0.238555, 0.0243645, -0.118378, 0, 0, 0.0582126, -0.677283, -0.0401264, 0.224911, 0.185033, 0.0645809, 0.153178, 0.674526, 0.00360232, 0.14202, 0.107262, 0.352631, -0.624968, -0.0937306, -0.0104478, -0.00260989, 0.00895039, 0.0166393, 0.025115, 0.0110771, -0.00747793, 0.0100261, 0.0106881, -0.0187887, -0.0231109, -0.0391486, -0.0401616, -0.0520729, -0.0542986, -0.0136061, -0.0339828, -0.0163564, -0.0209171, -0.0137007, -0.0111122, -0.0168658, -0.0151461, -0.0177419, -0.0194895, -0.0109267, -0.00645687, -0.0271427, -0.0110006, -0.0139873, 0.0103592, -0.00472768, -0.0121756, -0.0200751, -0.0233813, -0.0254527, -0.026002, -0.0032364, -0.0170319, -0.0109456, -0.00512843, -0.0134821, -0.0242582, -0.0824398, -0.0437587, 0.00314152, 0.0221427, 0.00181748, -0.0858541, -0.0319657, -0.0267471, -0.0111452, 0.00193025, -0.0301443, -0.0476312, -0.21421, -0.237852, 0.0538262, -0.217905, 0.115745, 0.197977, 0.41908, 0.109282, 0.0450541, -0.0801659, -0.151869, 0.810727, -0.148143, 0, 0, -0.149005, 0.0629851, -0.397976, 0.254674, 0.0623006, 0.0637515, 0.128198, -0.0435821, -0.852342, 0.109341, 0.439576, -0.0463116, -0.0912317, -0.368807, 0.0034688, 0.0127783, 0.0195323, 0.0179932, 0.016437, -0.0018744, 0.00864789, 0.0193761, -1.40438e-05, -0.011347, -0.00322869, -0.0232683, -0.0411421, -0.102459, -0.0240361, -0.00742991, 0.00345193, -0.00736629, 0.00013136, -0.00147853, 0.00550115, 0.00316263, -0.0219478, -0.0227429, -0.0261261, -0.0188284, -0.00128404, -0.00519411, -0.00653749, -0.0131733, -0.00296641, -0.0210122, -0.0117311, -0.0352784, -0.0245956, -0.0219934, -0.0306295, -0.0931533, -0.0184347, -0.0157466, -0.0333768, -0.023154, -0.0345488, -0.0352428, 0.00657606, 0.011029, 0.0163309, 0.0159192, 0.00633834, 0.0125775, 0.0173415, 0.0070476, -0.0117764, -0.0233222, -0.0419659, 0.284021, 0.372604, 0.0354918, 0.111692, 0.111916, 0.170642, 0.0434002, -1.16837, -1.43394, 0.0760999, 0.0352572, -0.238201, 0.201932, 0, 0, 0.349381, 0.548167, -0.103777, 0.128282, 0.0258535, 0.311805, -0.772536, 0.34863, 0.234626, 0.584404, -0.00086918, -0.509356, -0.169556, -0.0212644, -0.0293681, 0.0101518, 0.00587397, 0.00484255, 0.00228245, -0.00255161, -0.020548, 0.00783727, 0.00514817, -0.000889175, -0.0215392, 0.00397434, 9.68743e-05, -0.0147105, 0.00876086, 0.0111152, -0.00127382, 0.0128368, 0.00777634, 0.0185357, -0.0115278, -0.00946887, -0.0031127, 0.00565104, 0.00366178, 0.0216641, 0.00155357, 0.010362, 0.00323724, -0.000550742, -0.0175364, -0.0145084, 0.00355708, -0.0153386, -0.0102002, -0.00289498, -0.0402318, -0.120294, -0.0149911, -0.00358008, -0.01566, -0.018867, 0.0203342, -0.020639, 0.0189313, 0.0246426, 0.00607482, 0.0370291, 0.0167633, -0.0104402, 0.0125157, 0.020922, 0.00199972, 0.0307752, -0.0287277, 0.141036, 0.0904644, 0.0608674, -0.288506, 0.112161, 0.241235, -0.268708, -0.754536, -0.0551649, -0.0571079, 0.411868, 0.101936, 0.317265, 0, 0, 0.324524, 0.0512216, -0.869054, 0.26189, -0.0221858, -0.664671, 0.129146, 0.304984, 0.492825, -0.961726, 0.148178, 0.240189, 0.139214, -0.000599905, 0.0045413, 0.00243284, 0.00212942, 0.00959453, 0.0237846, -0.0137227, 0.0051206, 0.00756527, 0.000585311, 0.00649624, -0.009288, 0.0103225, 0.00764184, 0.0243625, 0.0141825, 0.0180862, 0.0207928, 0.00140987, -0.000722324, 0.00904813, 0.00461299, -0.0117866, 0.000642397, 0.0132322, -0.00217985, 0.00149148, 0.0041618, 0.00022058, -0.00986597, 0.00894904, -0.00215176, 0.000863158, -0.0082965, 0.00318768, -0.0251415, -0.00907946, -0.00195626, -0.0309595, -0.016943, -0.0210426, -0.00690423, -0.0228826, -0.00383571, -0.00295506, -0.0245962, -0.000888695, -0.00871698, 0.00164479, 0.000584116, -0.0166484, 0.0159204, -0.0434675, -0.0439792, 0.0106312, 0.000553152, -0.235167, 0.0592218, -0.107785, 0.646729, 0.196413, 0.110825, 0.286591, -0.291138, 0.274011, 0.820889, -0.915517, 0.0161361, 0.220222, 0, 0, -0.527789, -0.375963, -0.597302, 0.341062, 0.4012, 0.0933097, 0.618139, 0.405443, 0.172557, 0.112695, 0.300922, -0.00229678, -0.534867, -0.0340594, -0.0109741, 0.0114472, 0.00748194, 0.0386859, 0.01062, 0.0132854, 0.0400761, -0.0125574, 0.00344284, -0.00427309, -0.0108575, -0.0166886, 0.0128046, 0.0243269, 0.00656951, 0.00957826, -0.00283338, -0.0149134, -0.00174615, 0.000658239, -0.00895074, -0.00980118, -0.00937855, 0.00025507, 0.000187167, -0.00263607, 0.000998072, -0.00313664, 0.0066013, -0.00662865, 1.32108e-05, -0.00413439, -0.00117306, -0.0136484, -0.00287215, -0.00695207, -0.013067, 0.002096, -0.0127748, -0.00911067, -0.00359717, 0.0048019, -0.0209428, 0.00377006, 0.0103172, 0.0143653, 0.0171553, 0.0206745, 0.00430476, -0.00803953, 0.00559794, -0.0229309, 0.00166087, 0.0105747, 0.0449984, -0.825739, 0.0397767, 0.0644896, 0.489336, 0.0972857, 0.0320448, 0.342261, -0.0269036, 0.407357, -0.207035, 0.264746, -0.00989524, -0.171035, 0, 0, -0.252243, -0.174639, 0.209388, -0.436956, 0.0872885, -0.0128231, 0.547446, 0.0577526, 0.259665, 0.125562, 0.0752567, 0.746471, 0.0640558, 0.0195796, -0.0140465, -0.00389458, -0.0255465, 0.0164654, 0.0185107, 0.00246608, -0.0236871, -0.00602808, -0.0100688, -0.0184988, 0.00102929, -0.00940506, -0.0128097, 0.0180555, 0.00637732, -0.0130332, -0.00496961, 0.000218873, -0.0478884, -0.00830851, 0.000420491, -0.00295763, -0.0132305, -0.00841929, -0.0106486, -0.00893697, -0.000514895, 0.0185993, -0.000708002, 0.00230009, 0.006782, 0.0104328, -0.0115878, -0.0189136, -0.00492618, -0.00975611, -0.00149994, 0.0122467, -0.00432573, -0.0154742, 0.000789921, 0.0106098, -0.00951343, 0.0104115, 0.0133421, -0.000252981, 0.0118284, 0.0209793, 0.00616454, 0.0208508, 0.00889434, 0.0118881, 0.0231635, 0.0456518, 0.00278998, 0.117641, 0.437943, 0.0580519, 0.113652, 0.0756894, 0.135036, 0.477097, 0.817835, -0.280069, 0.151146, -0.300452, -0.33053, 0.458084, 0, 0, 0.258049, -0.10996, 0.38686, -0.0903657, -0.251228, 0.0159384, 0.198088, 0.229269, 0.190125, -0.109232, 0.121182, 0.0921239, -0.639756, 0.00831158, 0.00638323, -0.0239496, 0.0187472, 0.0228243, 0.0169053, 0.0266967, 0.0226034, 0.00997046, -0.0031701, -0.0142202, -0.0138747, -0.00277689, -0.0375318, -0.01098, -0.00136105, -0.00317379, -0.0204811, -0.035712, -0.126027, -0.0172293, -0.0148488, -0.00279888, 0.00432752, -0.00589616, 7.86803e-05, 0.00331925, 0.00950947, -0.000226324, -0.00256936, -0.00796588, 0.00130853, 0.000577324, 0.00359349, -0.00237721, -0.00901151, 0.0151549, 0.0100259, 0.021811, 0.0099344, 0.0203133, 0.0136045, 0.0161937, 0.0265552, 0.0280519, 0.0126925, -0.00725133, 0.0143286, 0.011465, 0.0284734, 0.0161123, 0.0195033, 0.0154779, 0.0281941, -4.69488e-05, -0.00793454, -0.20898, -0.289245, 0.0286111, 0.724074, 0.134825, 0.164264, 0.107051, 0.269104, -0.0281937, 0.0977092, -0.083442, -0.0595257, 0.308331, 0, 0, 0.327923, 0.0019805, -0.720545, -0.414863, 0.283291, 0.100402, 0.0112374, -0.279413, 0.454741, 0.441493, 0.300302, 0.143606, 0.00878702, 0.0198656, -0.00206299, 0.0312613, 0.0355557, 0.00307648, 0.0240953, 0.0427668, -0.00514442, 0.0223993, 0.00669511, 0.0220188, 0.00692376, 0.00555753, -0.00677271, -0.00396415, -0.0126264, -0.0113563, -0.00650729, -0.0231885, -0.0398175, -0.0206522, -0.0153313, -0.00686825, -0.00603441, -0.010141, -0.0138921, 0.00518041, -0.00367169, 0.0126457, -0.00202759, 0.0156966, -0.000479483, 0.00186913, -0.00309882, 0.00839385, 0.011362, 0.00925659, 0.019363, 0.0302874, 0.0183295, -0.000148583, 0.00832763, 0.00131283, 0.0154791, 0.0364143, 0.0242529, -0.0168736, -0.090992, -0.0502458, -0.000628408, 0.00255429, 0.0132123, 0.0298362, 0.0264994, 0.00237816, -0.0247489, 0.0261069, 0.136411, 0.0591942, 0.149884, 0.14917, 0.838332, 0.108652, -0.139244, 0.22347, 0.237504, 0.359542, -0.262469, 0.184822, 0, 0, -0.179198, -0.242927, 0.0332168, -0.705803, -0.140064, -0.0325254, 0.0423449, 0.198226, 0.180125, 0.124127, -0.772541, 0.0578202, 0.122396, -0.0107833, -0.00764908, 0.0148617, 0.03629, 0.0157593, 0.0326257, 0.00257147, 0.0219593, 0.00358422, -0.000342333, 0.0253596, 0.0213976, -0.0201152, -0.0116306, -0.0151375, -0.00436474, -0.0200572, -0.0181772, -0.0146341, -0.0274535, -0.0182532, -0.0338307, -0.0202529, -0.0158803, -0.0255223, -0.021012, -0.015451, -0.000681831, -0.0070789, 0.00638249, -0.0168642, -0.0010748, 0.00687209, -0.0132873, -0.0126722, 0.00695684, 0.0034193, -0.0114533, -0.00508856, -0.00202927, -0.0160157, -0.0135136, 0.00111996, -0.0147816, 0.0412073, 0.0185534, -0.00645421, -0.141498, -0.09473, 0.00770746, 0.0133175, 0.0246089, -0.0253598, 0.0178559, 0.0350896, 0.00675551, 0.0495666, -0.126618, 0.238742, 0.128471, -0.151196, 0.187027, 0.116789, 0.816726, 0.423063, 0.150325, 0.329924, 0.0686176, -0.12871, 0, 0, -0.790043, -0.0436484, -0.0231587, -0.0751625, 0.0669419, 0.32547, -0.568876, 0.102374, 0.506941, 0.099192, -0.224013, 0.110274, 0.272492, -0.00691804, 0.0171648, 0.0320725, 0.0023312, -0.00637585, 0.0107958, 0.00165103, -0.0135735, -0.0146549, -0.0232639, 0.0293095, 0.0158174, 0.00371131, 0.00534691, -0.00469093, 0.00956895, -0.0314848, -0.0263641, -0.00714211, -0.0082658, -0.00376917, -0.0140983, -0.0101231, -0.0178753, -0.0098988, -0.014898, -0.00312451, 0.00497636, 0.0107372, -0.00983192, -0.0114541, -0.00120772, -0.0108194, -0.00145058, -0.0129609, -0.00246547, -0.00817049, -0.00847834, -0.00345153, -0.00828269, 0.000182285, -0.0228914, -0.1553, -0.0555437, 0.0146054, 0.0371486, -0.00753543, -0.015263, 0.0104745, 0.0130016, 0.00193484, -0.00600557, -0.0187217, 0.00961597, 0.0246148, 0.00485413, 0.0838111, 0.261641, 0.0512976, 0.284092, -0.0515052, 0.16092, 0.221995, 0.122022, 0.50957, 0.354616, 0.278882, 0.00911812, 0.174113, 0, 0, -0.142608, -0.473973, 0.260076, 0.577445, 0.158002, -0.789751, 0.0630608, -0.907179, 0.180801, 0.185225, -0.0395476, -0.255478, 0.116641, 0.0482321, 0.0216924, -0.0179895, 0.00789975, 0.0233117, 0.00569581, -0.010834, -0.00225424, -0.0381713, -0.00165343, -0.00979405, -0.0107035, 0.015717, 0.012068, 0.00107162, -0.0189066, -0.0441, 0.00493853, 0.00711677, 0.00497068, 0.0100147, 0.000421968, -0.00488272, 0.0114384, -0.00550405, -0.0187456, -0.00613836, 0.0071046, 0.0186116, 0.0100106, 0.00972027, 0.0105869, 0.0104593, 0.00653103, 0.0182948, 0.0089735, -0.00334512, 0.00338925, 0.0177427, 0.00371169, 0.0143907, -0.00365643, -0.0471015, 0.00765304, 0.00612988, -0.00318734, 0.016099, 0.0219401, 0.0299003, 0.0372284, 0.00327138, 0.000739901, -0.0139928, -0.00336626, -0.00904173, -0.0113904, 0.111034, 0.34717, 0.073593, -0.0837113, 0.127281, 0.211062, -0.540709, -0.953802, -0.922821, 0.231688, -0.0659497, -0.0616399, 0.216807, 0, 0, -0.231899, 0.120287, 0.342809, 0.184051, 0.216721, 0.115876, 0.280319, 0.0585534, -0.25716, 0.133406, 0.0898762, 0.133465, -0.164925, -0.0260477, 0.0319221, -0.00924896, 0.0327439, 0.0081259, 0.00812514, 0.00906624, 0.0130866, -0.00820219, -0.00254268, -0.00953427, 0.0379791, 0.0121408, 0.00636334, 0.00891732, 0.00574043, -0.00391695, -0.00292488, 0.00310478, 0.000651205, -0.00833666, -0.015844, -0.00842067, 0.00728017, -0.0079611, -0.00973909, 0.00680231, 0.0178134, 0.0121381, 0.0220592, 0.012776, 0.00615364, 0.0100241, 0.017079, 0.0347446, 0.0111052, 0.00151809, 0.0304804, 0.0183354, 0.0177766, 0.0309215, -0.00688228, 0.011111, 0.0183224, -0.00361279, -0.0269019, -0.0174767, 0.0183237, 0.0167074, 0.010809, -0.0105129, 0.0101515, -0.00078698, -0.0134195, -0.0142643, 0.00937275, 0.202146, -0.109196, 0.184368, 0.113088, -0.917951, 0.16163, 0.141977, -0.282488, 0.415931, 0.117747, -0.133912, 0.00653495, 0.284463, 0, 0, -0.423232, -0.0355822, -0.832336, -0.0898858, -0.303356, -0.637845, 0.0824857, 0.206523, 0.305548, 0.090308, -0.281191, 0.582686, 0.127617, -0.0376906, 0.0141936, 0.0259807, 0.0181768, 0.0302547, 0.0249619, 0.026655, 0.0306426, 0.00058798, 0.0240652, 0.0112282, 0.0400513, -0.00307203, 0.0172581, 0.00122126, -0.00833604, -0.0183584, -0.00655043, -0.00969775, 0.00142729, -0.000596852, -0.0149771, 0.00778961, 0.00291816, -0.00413354, -0.00684945, -0.00239971, 0.0116719, -0.00677096, 0.0120588, 0.00265392, 0.0015896, 0.00417345, 0.0123648, -4.7221e-05, 0.00175273, -0.00732251, 0.0121317, 0.0219949, 0.00346021, -0.00578448, -0.000940132, 0.0181176, -0.0270212, -0.10629, -0.0111753, 0.0118945, 0.0108137, 0.0330215, 0.0175839, 0.00878999, 0.0060732, 0.0160803, 0.0129071, 0.0397111, 0.0391812, 0.103442, 0.294124, -0.383127, -0.16153, 0.106085, 0.171934, -0.343576, -0.0565759, 0.270827, 0.099639, 0.0419931, -0.120815, -0.139356, 0, 0, -0.166841, 0.119768, 0.127523, 0.0320555, 0.0690333, -0.130769, -0.662898, -0.684058, 0.261299, 0.0836482, -0.119971, -0.497337, -0.203136, -0.0408619, -0.043792, 0.00510939, 0.0106793, -0.00106937, -0.011294, 0.0089281, 0.00542222, -0.0102662, 0.0136485, 0.0316403, 0.0317528, -0.0130228, 0.00212269, 0.00568606, 0.00209619, -0.00507779, -0.00377003, -0.00192216, 0.0111339, -0.00749373, 0.00941502, 0.00321006, 0.0150422, 0.00575449, 0.0148411, 0.0198468, 0.00303947, 0.00410847, -0.000899428, -0.00325632, -0.00107206, -0.00042739, -0.00643614, 0.000664354, 0.0146403, 0.000602607, 0.00645572, -0.00208926, 0.0170201, 0.00759848, -0.00457699, -0.0206398, -0.0197201, -0.00925335, 0.0537666, 0.0246116, 0.0246569, 0.0198384, 0.0190228, 0.00608674, 0.0042647, 0.00113791, -0.0171354, -0.000500687, 0.0178424, 0.341586, 0.035708, 0.510584, -0.164521, 0.0788744, 0.0975233, 0.586211, 0.182747, -1.28246, -0.295163, 0.440089, -0.126817, -0.184434, 0, 0, 0.396314, 0.352882, -0.804874, -0.216864, 0.0393822, 0.381717, 0.358005, -0.236743, 0.328749, -0.0985118, 0.12298, 0.240271, 0.0765664, -0.144965, -0.0140264, 0.00418598, -0.00517303, -0.0200829, 0.0051615, 0.0123298, -0.00586269, 0.0107822, 0.00751573, 0.0350194, 0.0352082, 0.042063, 0.00963596, 0.0145619, 0.00887474, 0.0166, 0.0214076, 0.021857, 0.0083676, 0.0277618, 0.0127435, 0.0178593, 0.0249379, 0.0223563, 0.0132848, 0.0190327, 0.00609476, 0.00383254, 0.00787664, -0.00724325, -0.000649987, -0.00420661, -0.00151822, 0.0188173, 0.0231722, 0.0282127, 0.0177031, 0.0104196, 0.0349742, 0.0168366, 0.00645567, 0.0187349, 0.0251009, 0.0233726, 0.0575355, 0.0268531, 0.0216276, -0.00510103, -0.00207521, 0.00664184, -0.0176575, -0.00446639, 0.00662781, -0.0169965, 0.0275475, 0.0815498, 0.342711, 0.0975815, -0.334484, 0.502057, 0.183947, 0.06581, 0.287413, -1.10843, -0.26191, 0.125567, -0.0632572, 0.132364, 0, 0, 0.312927, 0.552231, 0.457025, 0.331797, -0.282747, 0.166433, 0.128903, 0.243763, 0.208316, -0.181764, 0.235701, 0.144413, 0.0985668, -0.000253293, 0.0337894, -0.00558321, 0.0137007, -0.02045, 0.00214084, 0.00277311, 0.0202637, 0.0089942, 0.00752768, 0.0204006, 0.0330107, 0.0067269, 0.0217456, 0.0236499, 0.0120134, 0.00640936, 0.0167892, 0.0109196, 0.0132521, 0.011487, 0.0116916, 0.0098066, 0.0152824, 0.00426325, 0.00217087, 0.0071026, 0.00248377, 0.00218696, 0.00237551, 0.0163653, 0.0145026, 0.0192084, 0.00987485, 0.0191698, 0.00375276, 0.017005, 0.0300613, 0.045237, 0.0215641, 0.0128345, 0.0116372, 0.0226387, 0.0276963, 0.0456544, 0.0288245, 0.0342047, 0.0266019, 0.0102897, 0.00621582, 0.0223135, -0.0303064, 0.0194755, 0.0635817, -0.0198292, -0.00442743, -0.0677945, 0.141827, -0.519744, -0.50258, 0.126315, 0.206849, 0.163779, 0.265126, 0.189183, 0.0804783, -0.440753, 0.447448, 0.43253, 0, 0, 0.308296, -0.261571, 0.808024, 0.537938, 0.135697, 0.0494239, 0.0622563, 0.211061, 0.393671, 0.064731, -0.418943, 0.164936, -0.0858327, -0.0147029, 0.0223045, 0.0420405, 0.0122554, 0.0018962, 0.00233892, -0.00149495, 0.00471752, -0.00582272, 0.0114326, 0.0444885, 0.025622, 0.000346738, 0.0124722, 0.00615882, 0.0107076, -0.00202907, 0.00824047, 0.00477563, 0.00366422, 0.00677242, 0.004476, 0.0143418, -0.00350194, -0.0550737, -0.0976658, -0.0108059, -0.00752023, -0.0101487, 0.0226558, -0.000961889, 0.00810752, -0.00335478, 0.0209075, 0.00709618, 0.00189043, 0.00766768, 0.00494761, 0.0138593, 0.0124195, 0.0324555, 0.0146457, 0.00660776, 0.00543295, 0.0344669, 0.0368166, 0.0143329, -0.00603963, 0.00905754, 0.011312, -0.0137743, 0.00662024, -0.0138763, -0.000915638, 0.0138412, -0.0583015, -0.233798, 0.0795112, -0.0291741, 0.101157, -0.590659, 0.127197, 0.121784, -0.743093, 0.541947, -0.12673, 0.0641946, -0.105551, -0.0729321, 0, 0, 0.119332, 0.0258488, 0.205647, -0.193384, -0.0373294, 0.191259, -0.282976, 0.614752, -0.0677257, 0.133102, 0.575796, 0.0409288, -0.147464, 0.0134517, 0.00544352, 0.00880358, -0.0039089, 0.0254666, 0.0155534, 0.00510006, 0.00725125, -0.0198478, 0.000562712, -0.0151066, 0.000600463, 0.0112459, 0.00278443, 0.0129349, 0.0339065, 0.0114972, 0.00482203, 0.0132393, 0.0115146, 0.00262103, 0.00263984, -0.00580972, -0.0412845, -0.136654, -0.0263737, -0.00493225, 0.00885548, -0.00654245, 0.007128, -0.000564064, 0.0129321, -0.0112862, 0.00647276, 0.0116801, 0.003182, -0.0195492, -0.00837143, 0.0336604, 0.0159754, 0.0165559, 0.000509921, 0.0203692, -0.0151184, -0.00131489, 0.00643072, 0.011508, 0.0224052, 0.00234579, 0.016934, 0.0160728, 0.00801424, 0.0141489, -0.00322525, 0.0335973, 0.0227905, 0.0249154, 0.307633, -0.164854, 0.125346, 0.0877266, 0.908868, 0.316647, 0.147388, 0.244137, 0.308932, 0.421985, 0.419304, 0.361663, 0, 0, 0.145022, 0.0843777, 0.143698, 0.380051, 0.0460275, -0.866631, -0.816386, 0.257609, 0.220702, -0.732604, -0.749459, -0.931204, 0.0723077, 0.17377, -0.0112917, 0.00466412, 0.00201869, -0.0148108, -0.00995847, 0.00995571, 0.00244292, -0.00131839, -0.0121378, -0.0530479, 0.0164675, 0.0349169, 0.0366642, 0.0200382, 0.0205194, 0.0237291, 0.0410059, 0.0108858, 0.0102547, 0.0261026, 0.0154796, 0.01906, 0.00560915, -0.010388, 0.0167971, 0.022731, 0.00381355, -0.0107504, 0.00635017, 0.0113218, 0.0158179, 0.00293774, 0.00224373, 0.0118725, 0.0182059, 0.0160465, 0.00318743, 0.00647251, 0.010512, 0.00015965, 0.00271909, 0.00529789, 0.0141591, 0.0240249, -0.00722703, -0.00197782, 0.0224268, 0.0151718, 0.0159411, 0.0155714, 0.00591535, 0.0120488, 0.00585294, -0.00692558, 0.0314885, 0.303762, 0.288759, 0.0510943, 0.0359337, 0.137601, 0.530536, 0.758209, -0.888046, -0.380721, 0.0419405, 0.295098, -0.11401, 0.0312169, 0, 0, -0.195574, 0.170903, -0.489333, -0.0007724, -0.627435, 0.137721, 0.0874363, 0.269864, -0.608853, 0.1595, 0.00891245, 0.248017, 0.0136354, -0.0642345, -0.000432224, 0.0183073, 0.0180969, -0.00869057, 0.0155822, 0.0199923, 0.00962349, 0.000605561, 0.00337559, -0.00907364, -0.00128699, 0.0306051, 0.0158025, 0.0262591, 0.0262151, 0.0186505, 0.0183444, 0.0177188, 0.0260244, 0.0325355, 0.0143112, 0.0343368, 0.0106025, 0.0245445, 0.0307265, 0.0110393, 0.00397828, 0.010466, 0.00790279, 0.0185927, 0.00755156, 0.0186703, 0.0100034, 0.00764798, 0.0174677, 0.0193069, 0.00900468, 0.0135195, 0.0243092, 0.00592569, 0.0145046, 0.0104082, 0.039877, 0.00155149, -0.0202454, -0.00378919, 0.000499213, 0.00655144, 0.0113942, 0.0148824, -0.0100762, 0.00952784, -0.00360318, 0.013145, 0.0596471, 0.119787, -0.103966, -0.112682, 0.564031, 0.371698, -0.776339, 0.146128, 0.891279, 0.261403, -0.0556854, 0.0229187, -0.564731, -0.793658, 0, 0, 0.0805303, 0.474524, -0.514368, -0.914677, 0.377261, -0.137858, -0.107688, 0.0557272, 0.227152, 0.153803, 0.145261, -0.129111, 0.108924, 0.0346285, 0.00342663, -0.000274127, 0.0104689, 0.0320726, 0.0350468, 0.0186923, 0.00646141, 0.0170711, 0.017801, 0.00342161, 0.000615368, 6.90743e-05, 0.0213554, 0.0160914, 0.0132644, 0.00871429, 0.0181833, 0.0180404, 0.0274605, -0.00144676, 0.0205441, 0.0145818, 0.0121387, 0.0305136, 0.0106599, -0.000524844, 0.0102672, -0.000650982, 0.0037732, 0.00545955, 0.0109494, 0.0106008, -0.0100528, 0.00171471, 0.00288972, 0.00969334, 0.0190599, 0.0173427, 0.00950468, -0.00230879, 0.0122167, 0.0277691, 0.00784667, 0.028366, 0.000978042, 0.0132093, -0.00595368, 0.0204037, -0.000308575, 0.0072004, 0.00773115, 0.0218236, -0.000774972, -0.00381911, 0.0449146, 0.301945, -0.1146, 0.0561036, 0.360099, 0.141615, 0.122107, -0.308441, 0.229457, -0.429957, 0.250635, -0.446469, -0.105886, 0.045012, 0, 0, -0.702147, -0.552685, -0.0476603, -0.247463, 0.385857, -0.258473, 0.492528, 0.236909, -1.01504, 0.341306, 0.451742, 0.135897, 0.355751, 0.0208248, 0.0225391, -0.00631026, 0.0442304, -0.010116, 0.0129288, -0.023825, -0.00911755, -0.0121819, -0.0732196, -0.057863, 0.00942084, 0.0151085, 0.00809672, 0.0176158, 0.0175752, 0.0340834, 0.0244438, 0.0154995, 0.0301292, 0.0191577, 0.00624873, 0.0265766, 0.00427185, 0.0214081, 0.0211791, 0.00356188, -0.00812415, -0.0148526, -0.0101725, 0.00766837, -0.00289704, -0.000248814, 0.00501333, 0.00233655, 0.00362184, -0.00445276, 0.00687201, -0.0149696, -0.0271231, 0.0188288, -0.0106079, 0.0142266, -0.0048516, 0.0172048, 0.0280456, 0.0101162, 0.0150063, 0.0257783, 0.0196345, 0.00566553, 0.0214078, 0.000758709, 0.0291914, 0.0277476, 0.0269809, 0.332162, -0.447149, 0.0351684, 0.12037, 0.0587313, 0.231384, -0.550643, 0.380497, -0.213039, -0.294405, -0.221807, -0.462884, -0.160887, 0, 0, 0.352861, -0.0177214, 0.826166, -0.209085, -0.320636, 0.207346, -0.113187, 0.178333, -0.30614, 0.148509, 0.167544, -0.122687, 0.0664346, 0.00830336, -0.000177063, 0.00878489, 0.00485298, 0.00351259, 0.0111435, -0.0132331, -0.00700551, -0.0191985, -0.123005, -0.0831379, 0.00345694, 0.0157416, 0.0233914, 0.0266269, 0.0327848, 0.0208024, 0.0273495, 0.0259684, 0.0265568, 0.0300798, 0.0175194, 0.0223294, 0.0295316, 0.0230221, 0.00114046, -0.023591, 0.0128242, 0.0106171, 0.00810288, 0.00894293, 0.0034591, -0.000348431, -0.00625996, -0.00544833, -0.0130855, 0.0148401, 0.00891006, 0.00250449, -0.132251, -0.0174491, 0.0199907, 0.0241792, 0.0388134, 0.0217316, 0.016451, 0.02707, -0.00136912, 0.000407884, -0.00155862, 0.0212819, 0.0161212, -0.00700427, -0.0175543, 0.00712348, -0.0651838, -0.0157805, -0.290724, 0.136797, 0.138653, 0.131845, -0.0112648, 0.141139, -0.0355231, 0.360894, 0.471183, 0.201921, 0.369957, 0.0159194, 0, 0, 0.14227, -0.542195, 0.529014, -0.0150997, -0.0266598, 0.408821, 0.24345, -0.501852, 0.23043, 0.233686, 0.0087078, 0.60591, -0.33832, 0.0194164, -0.0482456, 0.00683775, 0.0335082, 0.0317978, 0.00944204, 0.00806831, 0.0144777, 0.0426317, -0.00669975, -0.000745386, 0.00485706, 0.0298259, 0.00286774, 0.032206, 0.0138541, 0.0286037, 0.0289788, 0.0246753, 0.0322342, 0.00858591, 0.0230759, 0.0190396, 0.0164627, 0.0346309, -0.020969, -0.108986, -0.0115992, -0.00216174, 0.0187257, -0.00367979, 0.00645574, 0.02059, 0.0137362, 0.0114531, -0.000338723, -9.06518e-05, 0.0107523, 0.0260175, -0.0136344, 0.0125533, 0.013751, 0.00488135, 0.0188861, 0.018631, 0.00677977, -0.00786771, -0.0661716, -0.118971, -0.00981952, -0.00771333, 0.00671516, 0.00962752, 0.00170629, 0.0433712, 0.0097494, 0.752341, 0.537573, -0.539717, -0.155735, -0.674348, 0.20786, 0.233364, 0.252056, -0.0274965, 0.344996, -0.0659735, 0.513134, 0.580377, 0, 0, -0.239654, -0.0758966, 0.266596, 0.289596, -1.24218, -0.5426, 0.278001, 0.193215, 0.216115, 0.119506, 0.108203, 0.119015, 0.165553, 0.0131031, 0.0239489, 0.0178622, 0.0354085, 0.0246704, 0.0225036, 0.00666471, 0.0113176, 0.0073618, 0.0147784, 0.0134891, 0.000146234, 0.0154446, 0.0236002, 0.0305504, 0.0342673, 0.0348545, 0.0385092, 0.0373301, 0.0201641, 0.0125271, 0.0295687, 0.00918832, 0.0125794, 0.0075003, -0.00162665, -0.0122077, -0.0173036, -0.00188784, -0.00997207, -0.00400357, -0.00537845, 0.0116782, 0.0124694, 0.000789181, 0.0139885, 0.0167844, 0.00911192, 0.0396211, 0.0052568, -0.00263833, 0.00162728, 0.0249571, 0.00617718, 0.00298785, 0.0169685, 0.00632317, -0.064717, -0.0864975, -0.00827016, 0.000964355, 0.00614817, -0.00862281, 0.0229462, 0.0228515, 0.0195085, -0.00516469, 0.27211, 0.409192, 0.219127, 0.146134, -0.338684, 0.244516, 0.0311012, 0.100021, 0.451678, 0.0883936, -0.15353, -0.173334, 0, 0, -0.290443, 0.116101, 0.138695, 0.208643, -0.0378117, 0.101792, 0.072878, 0.197939, -0.0715397, 0.124935, 0.111316, -0.140974, 0.156958, 0.033123, -0.00766518, 0.0283144, 0.00480889, 0.0229772, 0.00747709, 0.00654937, -0.000868612, 0.0141089, 0.000659792, 0.013989, 0.0158918, 0.00903451, 0.0252594, 0.0190212, 0.0241154, 0.0291279, 0.0260783, 0.0142613, 0.030153, 0.025105, 0.0129642, 0.00352343, -0.00228727, 0.0145036, 0.0119527, -0.0165307, -0.0779897, -0.0235052, -0.00964113, -0.00284949, -0.00268456, 0.0171599, 0.024492, 0.00808088, -0.00854303, -0.00409035, -0.0102674, 0.0210914, 0.0153817, 0.00787349, 0.0125609, 0.00400798, 0.0104452, 0.0262228, 0.0381991, 0.00445668, 0.00258841, -0.00697198, 0.00657253, 0.028033, -0.0083342, -0.00360611, 0.0203619, 0.00138239, 0.0490661, 0.109895, 0.196019, -0.636403, 0.121568, 0.0796361, 0.137003, -0.524521, 0.178005, -0.0261073, -0.0734116, -0.280662, 0.339382, -0.333702, 0, 0, 0.110043, -0.0799681, -0.994239, 0.103738, -0.133187, 0.835193, 0.75731, 0.246066, 0.23398, -0.655469, 0.0322699, 0.0752283, 0.00613013, -0.0602498, 0.0322961, 0.0128974, 0.00668863, 0.042636, 0.00137754, -0.00560987, 0.0122685, 0.00887011, -0.0259053, 0.00473034, 0.0149468, 0.0156961, -0.0231823, 0.0131505, 0.0294266, 0.0497114, 0.0551901, 0.040342, 0.0594725, 0.0354651, 0.0451912, 0.00074545, 0.00357702, 0.0271149, 0.0431104, -0.0219182, -0.0470347, 0.00403028, 0.00934395, -0.019174, 0.0184895, 0.00707531, 0.0119705, 0.0133702, 0.015473, 0.0149602, 0.0191502, 0.00329596, 0.0208222, 0.0139958, 0.0167524, 0.0204249, 0.0304201, 0.0355056, 0.0133922, 0.00399478, 0.0057699, 0.0295299, 0.00786322, 0.0163553, 0.0142437, -0.016668, 0.0339098, -0.0226997, 0.00145475, -0.235698, 0.494363, 0.109546, -0.224065, 0.129013, 0.218321, -0.870067, 0.268142, 0.0784171, -0.146728, -0.449583, -0.104838, 0.27386, 0, 0, -0.106002, 0.27531, 0.148556, -0.702765, -0.891822, -0.4761, -0.219924, 0.523815, -0.339917, -0.088642, 0.143906, 0.465898, 0.212268, -0.0541829, 0.0278194, 0.0176617, -0.0120979, 0.0194753, -0.00672207, 0.00995056, -0.0025048, 0.0236662, 0.0170649, 0.0400912, 0.00998179, 0.0269631, 0.00746627, 0.0130785, 0.0286343, 0.0468336, 0.0412146, 0.0220454, 0.0216783, 0.0403545, 0.0398768, -0.00136286, 0.00940393, 0.00260998, 0.0270998, 0.0149207, 0.00927306, -0.000526352, 0.0233053, 0.0167807, 0.0079517, -0.0082127, 0.0229947, 0.0132546, 0.0266178, 0.0292488, 0.0216173, 0.00122105, 0.0210155, 0.0209969, 0.0112666, 0.0181056, 0.0519767, 0.0297803, 0.0287513, 0.0195414, 0.0245209, 0.00167551, 0.0089548, -0.00725823, 0.0118813, -0.0234004, 0.00763914, -0.0138497, 0.00190379, -0.13938, 0.369499, 0.120837, 0.117554, 0.034209, -0.143511, 0.140117, 0.292607, 0.372597, -0.146702, 0.225791, 0.0314292, 0.0915597, 0, 0, -0.170366, 0.0929373, -0.305982, -0.421973, 0.155267, -0.731054, -0.401238, 0.246787, 0.185659, 0.0770779, 0.00357011, 0.130148, 0.117085, -0.0213994, 0.0151246, 0.0151559, 0.0270413, 0.0245415, -0.00242166, -0.016132, 0.0060456, 0.0163159, 0.00813071, 0.0231943, -0.000632065, -0.00302836, 0.0232412, 0.00426581, -0.00117107, 0.0153372, 0.0155899, 0.0202449, 0.00827389, 0.00691851, 0.0245794, 0.00912181, 0.00631397, 0.0228423, 0.00881941, 0.0151588, 0.00396746, -0.0116824, -0.00380213, 0.0334631, 0.0147629, 0.0135369, 0.00292757, 0.0130498, 0.0096923, 0.0293863, 0.0103856, 0.0243964, 0.0167266, 0.00903693, 0.0119849, 0.0273706, -0.0153759, 0.00524803, -0.00845425, 0.0199563, 0.0166447, 0.00614369, 0.013759, -0.00369641, -0.0207142, 0.00244241, 0.00704693, 0.0165003, 0.0263125, -0.126921, 0.33525, 0.275101, 0.0964794, -0.0279062, 0.276295, -0.0898866, -0.114476, -0.154717, -1.16924, 0.16427, 0.140561, -0.0656118, 0, 0, -1.01671, 0.579686, 0.0587478, 0.228033, 0.0845411, 0.158347, 0.14762, 0.245729, 0.307352, 0.098835, -0.245373, 0.110077, 0.0135818, -0.0193209, -0.0265162, 0.00135194, 0.022458, -0.00182074, 0.00402217, -0.0036545, -0.00551382, 0.00614064, -0.000955415, 0.0060087, -0.0237173, -0.0217515, 0.00999508, 0.0240216, 0.0214199, 0.0267233, 0.00446582, 0.00710324, 0.0174225, 0.0163336, 0.00596344, 0.0145639, 0.00268764, 0.0124111, 0.0172115, 0.0141331, -0.00358624, 0.000915042, 0.00704083, 0.00366454, 0.00677087, 0.0053697, 0.00191335, 0.00876987, 0.0062013, 0.0158351, 0.0140922, 0.0230853, 0.0188745, 0.0100694, 0.00398891, 0.0113014, 0.00948254, -0.00808068, -0.007938, 0.0164367, 0.0262808, 0.0371095, 0.02598, 0.0162768, -0.0122608, 0.00445659, 0.0238588, 0.00989636, 0.0105785, 0.0754624, 0.0608152, 0.0686849, -1.0486, -0.781631, 0.140351, -0.0715477, -0.587318, -0.0288828, 0.452986, 0.330053, -0.101335, -0.105717, 0, 0, 0.127108, -0.436722, -0.0953791, 0.287091, -0.081557, 0.159085, 0.159271, 0.252939, 0.203426, -0.0311945, 0.139957, 0.220478, 0.218088, -0.00630489, -0.0241939, -0.0302855, 0.00243636, 0.010349, -0.0242364, 5.05415e-05, -0.00111845, 0.0017867, 0.0177028, 0.00517896, 0.0255265, -0.0179121, 0.00248077, 0.014764, 0.0149373, 0.0192855, 0.0342892, 0.0297574, 0.0304321, 0.0218236, 0.0132966, 0.00998794, 0.00419806, 0.0257195, 0.0283988, 0.0164601, 0.025985, 0.0136651, 0.0294683, 0.0112426, -0.00114719, 0.00771533, 0.0102623, 0.0211178, 0.0245389, 0.0206099, 0.013759, 0.0170459, 0.0348332, 0.0317464, 0.0224309, 0.0215132, 0.0292723, 0.0087347, 0.0165729, 0.0069524, 0.0256648, 0.0288632, 0.0219788, 0.0053357, -0.0304864, 0.00155183, 0.0225444, -0.0085163, 0.0111772, 0.170346, 0.155889, 0.193488, -0.119053, 0.354333, 0.163508, 0.178113, -0.443607, 0.0281412, -0.375436, -0.0713834, -0.153396, -0.0730397, 0, 0, -0.144317, 0.215238, 0.361094, 1.21106, -0.48724, 0.0198489, -0.0417697, -0.871553, 0.110659, 0.141952, -0.628183, -0.607704, -0.10308, -0.0531648, 0.00441932, 0.0182038, -0.000532968, 0.0111387, -0.00869416, 0.00301953, 0.00599966, 0.00474786, 0.0169018, 0.0375952, -0.00757433, -0.00596454, 0.00494013, 0.0174411, 0.0334728, 0.0209036, 0.0248724, 0.0257831, 0.00735379, 0.00914221, 0.0190425, 0.007977, 0.000629332, 0.018455, 0.0250534, 0.0240776, 0.0144833, 0.0239961, 0.0147148, 0.026668, 0.00472217, 0.0244913, 0.0180529, 0.0246393, 0.0258784, 0.0253173, 0.0277045, 0.0242894, 0.0303319, 0.0375995, 0.0228655, 0.0216174, 0.0392676, 0.00971793, 0.0218066, 0.0266002, -0.00778094, 0.00787932, -0.00559273, -0.0144086, -0.0308442, -0.000968875, -0.026834, -0.0392371, 0.043311, -0.0270759, 0.346282, -0.497044, -0.136597, 0.100009, 0.132499, 0.280627, 0.243155, 0.0788073, 0.189903, -0.627269, 0.294971, -0.00785697, 0, 0, -0.197727, -0.0658932, 0.240466, 0.398735, 0.178563, 0.0398199, -0.217998, -0.34179, -0.614596, 0.0832511, 0.284929, -0.529445, -0.478151, 0.0252301, -0.035981, -0.0171645, 0.00498494, 0.00206238, 0.0102514, 0.00105775, -0.0112382, 0.00537105, 0.0213306, 0.0246392, 0.0110787, -0.000293954, 0.0121807, 0.0239736, 0.00747329, -0.000619344, 0.0218011, -0.00096476, 0.019401, 0.00604514, 0.0053252, -0.00111951, 0.00211318, 0.00586173, 0.0194137, 0.0176817, 0.00222861, -0.00139633, 0.0179067, 0.00507238, 0.0198391, 0.0112214, 0.025489, 0.0153514, 0.0211008, 0.00986471, 0.0173096, 0.0289893, 0.0237234, 0.0186775, 0.0159621, 0.0260032, 0.00969434, 0.0128414, 0.008913, 0.00303948, 0.00266832, 0.0129525, 0.0068498, -0.00382883, -0.0133436, -0.018195, -0.0273098, -0.00419936, -0.0381104, 0.0285236, -0.165354, -0.295351, 0.0833925, -0.00476403, 0.0720168, -0.250793, 0.217007, 0.195691, -0.186944, -0.453708, 0.35081, -0.0026836, 0, 0, -0.32808, -0.121352, 0.141308, 0.215901, -0.245852, -0.597758, -1.05336, 0.139165, 0.127865, -0.446391, 0.379135, -0.255498, -0.277629, 0.00813613, -0.0321147, -0.0260419, -0.0160939, -0.0177183, -0.0256895, -0.0286939, -0.0245268, -0.00291399, 0.0020037, 0.0261188, 0.0398927, 0.0164519, 0.035611, 0.0146476, 0.0149345, 0.0140955, 0.00141547, 0.010967, 0.0211361, 0.0155247, 0.00907935, 0.0239046, 0.0217193, 0.00930912, 0.0136236, 0.0149103, -0.00703879, 0.00344991, 0.0151991, 0.0118686, -4.01236e-05, -0.00350515, 0.00193806, -0.00313859, 0.00225791, 0.00232425, 0.0064871, 0.00759805, -0.0010561, 0.0182402, 0.0101571, -0.000441796, -0.0281388, 0.0168195, 0.0217755, 0.020193, 0.0150357, 0.0115244, 0.00942411, 0.00634666, -0.010433, 0.0170893, -0.0140055, 0.0253638, -0.0451425, -0.926435, 0.106855, 0.0497039, 0.0846448, 0.0383932, -0.625952, 0.281354, 0.0154588, -0.0217047, -0.0623516, 0.335349, 0.0970995, -0.163928, 0, 0, 0.324669, 0.380451, -0.366386, -0.193811, 0.0686272, 0.56081, 0.00954701, 0.205276, 0.150335, -0.203452, 0.0248191, 0.119685, -0.507679, 0.0046985, -0.0151448, -0.00250645, -0.0139632, -0.0253649, -0.0169595, -0.0332232, -0.0169742, -0.016523, 0.0252805, 0.0289095, 0.00451079, 0.0283578, 0.0164318, 0.0115766, 0.0161212, 0.0219283, 0.0258644, 0.0211903, 0.0233834, 0.0190075, 0.016035, 0.0249204, 0.015406, 0.0166733, 0.0321911, 0.0199894, 0.00889824, 0.0043297, -0.00797965, 0.0227313, 0.00966878, 0.0129838, -0.00176649, 0.00743944, 0.0222619, 0.0104927, 0.0191704, -0.0564721, 0.0117035, 0.0394497, 0.00377718, -0.00923991, 0.0139627, 0.0208916, 0.0187249, 0.00831879, -0.00687143, -0.00856069, 0.0173107, 0.00617303, 0.0129127, 0.00746416, 0.000879039, -0.000822699, 0.0200773, 0.081818, -0.122273, 0.0716159, 0.138905, 0.0819157, 0.155405, 0.414759, 0.937334, 0.338747, 0.167474, -0.375023, -0.274467, 0.165663, 0, 0, 0.51277, -0.0263845, 0.0907723, 0.00669375, -0.341688, -0.890349, 0.112485, -0.938655, -0.729866, 0.104192, 0.109093, -0.388052, 0.0861015, 0.0503686, -0.0428389, -0.0382358, -0.0076547, 0.00549689, 0.0118468, -0.0150376, -0.00826924, -0.0174529, 0.0103739, 0.0123119, 0.00297875, 0.021266, 0.00541264, 0.0199713, 0.0240477, 0.0168818, 0.0175718, 0.0250809, 0.00849503, 0.0110588, 0.00374775, 0.00942059, 0.009642, 0.0195695, 0.0261342, 0.0239312, 0.00292506, -0.0165568, -0.103226, -0.00124631, 0.0170844, -0.00786558, 0.0185933, 0.0166685, 0.0262211, 0.0217342, 0.0094534, -0.00208232, 0.0198119, 0.028506, 0.0343607, 0.0293222, 0.0258421, 0.0131981, -0.0352405, -0.0142973, -0.0241897, 0.00816537, 0.00375098, -0.0140577, 0.00502896, -0.00867037, -0.0185554, -0.0271369, 0.0324208, 0.0565613, -0.284287, -0.373908, 0.239431, 0.0804595, 0.120016, -0.897884, 0.288939, 0.146763, 0.103277, 0.110671, 0.161437, -0.111308, 0, 0, 0.469973, 0.0411889, -0.319425, -0.0876676, 0.26614, 0.139103, -0.00192692, 0.192285, -0.891133, 0.507099, 0.083909, -0.185286, 0.11711, 0.0503114, -0.012426, 0.00990384, 0.0158234, 0.0192235, 0.0177932, -0.00912072, 0.0197953, 0.0122192, 0.00798533, -0.0224509, 0.0104648, 0.00136809, 0.00971481, 0.0276621, 0.0166133, 0.0155641, 0.00896272, 0.0127167, 0.00756988, -0.0109924, -0.00299699, 0.00131092, 0.0196466, 0.0197113, 0.0132246, 0.0132148, -0.000436025, 0.00893798, -0.0020373, -0.00227778, 0.00996218, 0.00978594, -0.000986138, 0.00481682, 0.00919609, 0.0206084, 0.00831915, 0.0147308, -0.00721529, 0.0142272, 0.0139011, 0.0232286, 0.00417044, 0.010024, -0.021757, 0.00166067, 0.00951437, 0.0180115, 0.00944031, 0.0028884, -0.0275052, -0.00175132, -0.0102279, 0.0133575, -0.0185617, -0.102743, -0.335385, -0.575325, 0.121674, -0.155245, 0.103348, 0.0820822, 0.203062, 0.421896, -0.481253, -0.0509974, -0.0446212, -0.380843, 0, 0, -0.197725, -0.0514788, 0.141, -0.872965, 0.12821, 0.0886758, 0.186512, -0.281986, -0.51154, 0.0901854, -0.020594, 0.188758, 0.0849871, 0.018696, -0.00194144, -0.0197107, 0.00472987, 0.0125614, 0.00544391, -0.0101274, 0.00257488, -0.0124534, 0.00339615, -0.00759121, 0.00983669, 0.0189914, 0.0208413, 0.0390345, 0.0218295, -0.0112956, -0.016594, 0.0154262, 0.00762849, 0.00183015, 0.0123372, 0.00346603, 0.029095, 0.014305, 0.00298933, 0.000640444, 0.00550988, 0.00365077, 0.00315869, 0.00969303, -0.0104143, -0.00773419, -0.00437164, -0.000617444, 0.000855367, 0.00345216, 0.0145777, 0.00172399, 0.0169964, -0.00167891, -0.00203361, -0.0178602, -0.0197299, -0.0069926, -0.0215155, -0.00379904, -0.00156684, 0.0170956, 0.0285584, 0.0082879, -0.011176, -0.00169873, -0.00992598, 0.00193298, -0.0169969, 0.0907866, 0.102591, -0.496603, 0.0877421, 0.0405157, -0.942751, 0.26103, 0.111736, 0.147865, 0.176875, -0.595774, -0.105785, -0.12666, 0, 0, 0.11601, 0.109602, -0.42211, 0.175653, 1.00726, 0.105072, -0.534053, 0.187353, -0.364317, -0.278088, 0.0968671, 0.588157, -0.10674, 0.028425, 0.00736128, -0.0236981, 0.0115482, 0.0179697, 0.00199823, 0.0240861, 0.00171201, 0.0102051, 0.00085894, -0.0148487, 0.0196696, 0.0283232, 0.0226183, 0.014223, 0.0155412, 0.00853353, -0.0714738, 0.0134216, 0.0198948, 0.0223352, 0.0131046, 0.031782, 0.00463946, 0.0117998, 0.0306972, 0.0206826, 0.0132561, 0.00772545, -0.0023007, -0.0015318, 0.00910988, 0.0089012, 0.0025118, 0.016342, 0.0142151, 0.0158211, 0.012691, 0.0201111, 0.0168949, 0.0295184, 0.0140407, -0.00153951, 0.00869321, 0.0107625, -0.000491018, -0.00533711, 0.0299114, 0.0292826, 0.0028075, 0.0137695, -0.012011, -0.013709, -0.0123369, -0.0933073, -0.0280147, 0.00602922, 0.12924, 0.0322805, -0.231119, 0.0461764, 0.14912, 0.553964, -0.610229, -0.288836, -0.108988, -0.137392, -0.164148, 0.734696, 0, 0, 0.0218754, -0.271875, -0.0860625, -0.422906, 0.162246, 0.0737254, -1.13734, 0.161283, 0.45452, 0.104248, -0.516352, 0.114697, 0.102705, 0.000572892, 0.0220716, 0.0111478, 0.00744765, 0.0146806, 0.00736366, 0.00482508, -0.0154954, -0.0120425, -0.0129775, 0.00854091, 0.0216637, 0.0186535, 0.0193571, 0.0256516, -0.00194073, 0.025848, 0.0104914, 0.0290404, 0.0321083, 0.0145593, 0.00301465, 0.00457932, -0.00151504, 0.0116777, -0.00087613, 0.0171626, 0.0135994, 0.016101, 0.010647, 0.00636316, 0.0039803, -0.00242884, 0.0123492, 0.0167015, 0.0161929, 0.00747398, 0.0212887, 0.0256959, 0.00240396, 0.0172532, 0.00504944, 0.0193268, 0.00475044, 0.00197857, -0.0138601, -0.0101902, -0.00774043, -0.0120108, 0.015203, 0.0077325, 0.00285876, 0.0193486, 0.0128449, -0.00220699, 0.0321726, 0.358973, 0.0543295, 0.0733502, 0.0760461, 0.11416, 0.21404, 0.375543, -0.538835, -0.395812, 0.103893, 0.0845051, 0.0132767, 0.62944, 0, 0, 0.220672, -0.0827594, -0.477505, 0.201075, 0.0777501, 0.205226, 0.340954, 0.164338, 0.167706, -0.324423, 0.0696063, 0.0812903, 0.440959, 0.0292436, 0.0146105, 0.016951, 0.0155639, 0.0476891, 0.00317159, 0.00618074, 0.00025761, -0.00561604, -0.0027293, 0.0024575, 0.00747287, 0.00672896, 0.0208431, 0.00641529, -0.00910318, 0.000264495, 0.015908, 0.00970544, 0.00886638, -0.00447934, 0.00724256, -0.00134461, -0.00254015, -0.00396706, 0.00404893, 0.00240462, 0.00402241, 0.0129706, -0.00915147, 0.0085597, -0.00923241, 0.00491906, -0.00822143, -0.0039401, 0.0054932, -0.00621948, -0.0116093, 0.018482, -0.0049147, 0.00711381, 0.00181194, 0.0103533, 0.0029757, 0.00817914, 0.00695958, -0.000119002, -0.0300762, 0.0124755, 0.00306697, -0.00788878, -0.0179338, 0.017025, 0.00730556, 0.0181585, 0.0509637, 0.177705, -0.686899, 0.0462767, 0.101501, 0.0670618, -0.498621, 0.178531, 0.28862, 0.105774, 0.0765734, 0.10906, 0.100205, -0.175543, 0, 0, 0.343159, 0.194466, 0.126743, 0.299826, -0.0897381, -0.00381269, 0.32035, -0.434158, -0.317112, -0.378094, 0.113387, -0.731132, 0.151337, -0.0128055, 0.0476507, 0.0385859, 0.0210153, 0.0214244, 0.0127526, 0.022218, 0.0188005, 0.01374, 0.00586975, 0.00728893, 0.0120589, -0.00045664, -0.00343222, 0.0105607, 0.0014497, -0.0181072, 0.00761841, -0.0146827, 0.00293054, 0.00334077, 0.0102801, -0.00176466, 0.00321095, 0.000810621, -0.00627883, -0.00954531, 0.0135883, 0.019819, -0.00742899, 0.00141036, -0.000453509, -0.00168244, -0.0081388, -0.00787068, -0.00212273, -0.00164617, -0.010528, 0.0010865, 0.00394147, 0.00346885, 0.0132204, 0.00892171, -0.0118617, 0.00645227, 0.016221, 0.0117577, 0.00929914, -0.00914847, -0.00108177, 0.00070527, -0.00849699, 0.0185456, 0.00770226, 0.00241712, -0.0261796, -0.438404, -0.419013, 0.0329568, -0.0648615, -0.825614, 0.191214, -0.413647, 0.16382, -0.245602, 0.256498, 0.481031, -0.230335, 0.226724, 0, 0, 0.160492, 0.379021, 0.125051, -0.394023, 0.453105, 0.135272, -0.691447, 0.173982, 0.20958, 0.133557, -0.401771, 0.144058, -0.0989494, 0.0103445, 0.027696, 0.0128143, 0.0285259, -0.018086, 0.0381917, 0.00620039, 0.00108761, 0.0128052, 0.0132693, -0.000815501, -0.0010254, 0.0258138, 0.0027362, 0.00993157, 0.00989159, 0.0203851, 0.0131995, 0.0145528, 0.0289504, -0.000174347, 0.00327849, 0.011665, 0.0216385, 0.00460951, 0.0190338, 0.0118726, 0.00462418, 0.00898782, 0.0125754, 0.000703935, -0.00342692, 0.00904587, -0.00878324, 0.00139004, -0.00197914, -0.00271453, 0.00433298, 0.0177375, 0.0130398, 0.0177348, 0.00116581, 0.0130944, 0.0115954, 0.000445021, -0.00761442, 0.0145831, 0.0194089, 0.00666908, -0.0151925, -0.0112687, 0.0231923, 0.00101864, -0.00277426, 0.00530248, -0.0280954, 0.0153765, 0.048734, -0.237454, 0.135381, 0.0856597, -0.656745, 0.471874, 0.168841, 0.652966, -0.213696, -0.145167, -0.0602996, -0.278303, 0, 0, 0.0224515, -0.2094, 0.372314, 0.834031, 0.0864222, 0.373337, -0.0614263, 0.125215, 0.207479, 0.270901, 0.083413, -0.084681, 0.144969, 0.0919517, 0.0268666, 0.00647108, 0.0268192, 0.0457779, 0.0320291, 0.00878686, 0.0192623, 0.0163947, 0.0151807, 0.0309426, 0.0084386, 0.0232073, 0.00770925, 0.022282, 0.00195581, 0.023307, 0.00862047, 0.00277666, 0.0117729, -0.0014108, -0.00584406, -7.99608e-05, 0.0073374, 0.00690555, -0.0101672, 0.0207443, 0.00982768, 0.00515935, 0.012673, -0.0051611, 0.00771916, 0.00671512, 0.0136044, 0.0119948, 0.00713987, 0.0231161, 0.0236856, 0.0176995, 0.0123142, 0.0211874, 0.00980152, 0.00980353, 0.0301264, 0.0274986, 0.00214791, 0.000749363, -0.000144481, 0.0211065, -0.00286758, 0.000623105, -0.00135806, -0.0319775, -0.0153421, 0.0223426, 0.0165537, 0.0864622, -0.504239, 0.0454867, -0.179949, 0.11304, 0.168487, -0.287551, 0.113776, 0.430614, -0.509781, 0.0111774, 0.0799922, 0.228464, 0, 0, -0.104959, -0.0807807, 0.376288, -0.142937, -0.127214, 0.271179, 0.327116, -0.109658, 0.27267, 0.0640966, 0.582092, 0.0455356, -0.202553, 0.0577131, 0.000839399, 0.0304346, 0.0426289, 0.020633, 0.0345184, 0.0248412, 0.0201712, 0.0112769, 0.0168296, 0.0271008, -0.000552233, -0.00765773, 0.0141092, 0.00756734, -0.00891695, 0.0102058, 0.00285531, 0.0101275, -0.00126422, -0.00235571, 0.00187726, -0.00190243, 0.000804765, 0.00643872, 0.00630239, 0.0182188, 0.000977798, 0.003668, -0.000759258, -0.013215, 0.00714472, 0.0128989, 0.0183075, 0.00209541, 0.00395886, 0.00544589, -0.0129467, 0.0226776, 0.0182191, -0.0120165, 0.00209588, 0.00230601, -0.00346487, 0.0234807, 0.0312905, -0.0103803, 0.0189289, 0.0100253, 0.00445346, 0.00746137, 0.00113026, 0.0127197, -0.0288126, 0.00822805, 0.0343931, 0.397321, 0.528501, 0.171375, 0.0947913, -0.0755955, -0.380315, 0.331751, 0.326954, -0.21674, 0.0178053, 0.111503, -0.114913, -0.685483, 0, 0, 0.073433, 0.505702, 0.0229995, -0.89758, 0.10656, 0.251101, 0.0861931, 0.176386, 0.140761, -0.243758, 0.145292, -0.0495936, 0.291125, 0.190122, 0.00337411, 0.0155446, -0.0174318, 0.00527125, 0.0094637, 0.0177896, 0.00597736, 0.00455204, 0.00203451, -6.58653e-05, 0.00847471, -0.00887262, 0.0135048, 0.00308286, 0.00356851, 0.00737411, 0.00643255, 0.00216593, -0.00734983, -0.000931123, 0.00291638, -0.0199057, 0.00341565, 0.0101332, 0.0194678, 0.0270493, 0.00995485, -0.00326814, 0.00308607, 0.00516288, 0.000813756, 0.00616111, 0.0012403, -0.00265278, -0.00727078, 0.00208868, 0.0108703, -0.00117744, 0.00444171, 0.00341918, 0.0221239, 0.000513972, -0.0206063, -0.023282, 0.00679237, -0.00429768, 0.0117833, 0.0338622, 0.00338659, -0.011251, 0.0125534, -0.0169482, 0.00063622, 0.0328315, 0.0140078, 0.160725, -0.168424, 0.0780409, -0.906393, -0.111844, -0.875514, 0.312649, 0.0430418, -0.0283411, 0.258995, 0.0519853, 0.262125, -0.818378, 0, 0, -0.164267, -0.0609866, -0.260096, 0.0558968, 0.458712, 0.0735322, -0.735341, 0.210442, 0.222445, 0.123322, 0.0847887, -0.159716, -0.590239, -0.00630243, -0.000272007, -0.0556336, 0.0196293, 0.00963864, -0.00185404, 0.0100667, -0.00227032, 0.00249417, -0.0212947, -0.0312975, 0.000845762, 0.0157883, -0.00241539, -0.0275214, 0.00259345, 0.0237026, 0.0270353, 0.018024, 0.028604, 0.016141, 0.00947002, -0.0550163, 0.0151148, 0.0170707, 0.0207178, 0.0302973, 0.0141577, 0.0117242, -0.000149357, -0.000947857, 0.000701661, 0.00893654, 0.00893102, 0.00711359, -0.00754681, -0.023661, -0.0189301, -0.013904, -0.0108278, 0.00679914, 0.0228444, 0.0024207, -0.000571654, -0.158785, -0.0365786, 0.019051, 0.0210119, 0.0220753, -0.00334971, -0.0165962, -0.00349948, 0.0109838, 0.0210903, -0.0330201, 0.0448968, 0.392194, -0.0835511, 0.0509066, 0.13672, -0.376382, -0.142133, 0.263579, 0.564148, 0.41448, 0.144795, -0.546305, -0.272911, 0.031192, 0, 0, 0.118088, -0.0489305, -0.131705, -0.212339, 0.0750245, -0.849964, 0.17661, -0.651461, 0.17185, 0.10683, -0.0729062, 0.118462, 0.107465, 0.0484597, 0.000558717, 0.0117413, 0.0452489, 0.00656822, 0.0169244, 0.00676164, -0.0159824, 0.000433659, -0.00936603, -0.0308871, -0.00225009, 0.0025264, -0.0302476, -0.145979, -0.012044, 0.000833357, 0.00533716, 0.00998125, 0.0232534, 0.00986184, 0.00292847, -0.00514859, 0.0210219, 0.00935195, 0.0143981, 0.0216698, 0.00704815, 0.00142274, 0.0104183, 0.00585931, 0.0115554, 0.00320029, 0.0146599, 0.0152716, -0.0010417, 0.0163638, -0.00292275, -0.148811, -0.0341117, -0.00820037, 0.00690723, 0.00580386, 0.0102776, -0.0287342, -0.0426742, -0.00346314, -0.014082, 0.0164168, -0.0151897, -0.0133228, -0.00402774, 0.0294466, 0.0167652, -0.0360319, 0.0711642, 0.083808, -0.184709, -0.0980398, -0.399685, 0.0804186, -0.782502, 0.271326, -0.613257, 0.14975, 0.180464, 0.0670264, 0.416608, 0.102782, 0, 0, -0.192021, 0.492874, 0.207674, -0.0196691, -0.294052, 0.197247, -0.0735843, -0.109156, 0.0555986, 0.168678, 0.0641907, -0.775262, 0.0994696, 0.0202677, 0.0309279, -0.0142657, 0.0301083, 0.00962297, 0.0210668, 3.58343e-06, 0.00281656, 0.00364458, -0.0223304, 0.0018766, -0.0110041, -0.0175352, -0.0073413, -0.0303934, -0.0216453, -0.0200569, -0.0275308, -0.00626038, -0.00309842, 0.00603172, -0.00668402, 0.00443937, -0.00862822, -0.00330077, 0.000224463, 6.62862e-05, -0.0195622, -0.000130437, -0.00770041, 0.00533275, -0.00391618, -4.66355e-05, -0.00771137, 0.0042634, -0.0167789, 0.000385853, -0.0228836, -0.0208459, -0.00342089, 0.00206059, 0.00132507, 0.0124075, -0.0238468, 0.0219538, -0.0301166, -0.0163151, 0.0055811, 0.0122646, -0.0191304, 0.0100754, 0.00806302, -0.0232621, 0.0211904, 0.00515403, 0.0764653, 0.0269575, 0.290038, 0.0375496, 0.070808, 0.0471113, -0.746099, -1.10055, 0.107236, -0.577217, -0.186472, -0.247479, -0.246867, -0.167438, 0, 0, -0.920212, -0.124836, 0.18137, -0.167062, 0.0663951, 0.0413478, 0.241166, 0.0485479, 0.299229, -0.291361, -0.0577546, 0.439095, 0.439336, -0.00685212, -0.018186, 0.0191875, 0.00134158, 0.0186317, -0.00180959, -0.000343055, -0.00525811, -0.0173829, -0.0120612, 0.0111842, -0.00162653, -0.0204498, 0.00168794, 0.0142601, -0.00473493, -0.0388172, -0.141292, -0.0100236, -0.00879179, -0.00842432, -0.00892723, -0.00195914, -0.0199702, -0.0159568, -0.00544021, -0.00495973, -0.0194822, -0.0116988, -0.007264, 0.00779178, -0.0165814, -0.00229094, -0.0110314, -0.00737096, -0.020545, -0.0168759, -0.0220858, 0.000520981, -0.00104014, -0.0331411, -0.0289936, -0.00340396, -0.0158464, -0.00914367, 0.0286717, -0.0155749, 0.0147171, 0.00936366, -0.00521418, 0.0139895, 0.0176361, -0.0226626, 0.00714659, 0.0373704, 0.0538818, 0.0928518, -0.293312, 0.473294, -0.698679, 0.0904287, 0.12754, 0.270273, 0.374576, -0.162374, 0.0202955, 0.0462198, -0.00906006, -0.0636285, 0, 0, 0.0810608, -0.00272496, 0.290353, -0.280434, -0.102417, 0.0764063, -0.135846, 0.157132, -0.992854, -0.125817, 0.123486, 0.41063, 0.030256, 0.143689, -0.0183929, 0.00164764, 0.00286223, 0.0127256, -0.00837642, -0.0179849, -0.0125724, -0.0120477, -0.00901109, -0.0259241, -0.0467739, 0.0352959, -0.00391901, -0.00286235, -0.0361554, -0.0297949, -0.0166319, -0.0121744, -0.00873186, -0.033098, -0.000429523, -0.000700627, -0.0187568, -0.0106145, 0.00145481, 0.00319789, -0.00822775, 0.00273214, -0.00301346, -0.00105057, 0.00484948, -0.00192215, 0.00767673, 0.0126274, -0.010013, -0.00737091, -0.0143242, 0.00279192, -0.00893212, -0.0140216, -0.0139853, -0.00715914, 0.0219237, 0.0131066, 0.0120845, -0.0171834, 0.0112012, 0.00481712, -0.0129329, 0.0207899, 0.0131395, 0.000301797, 0.00591426, 0.00601964, 0.110523, -0.00255375, -0.934664, 0.0938624, 0.153472, 0.106667, 0.174197, 0.106787, -0.349591, 0.361262, 0.0165563, 0.62917, 0.436654, 0.481975, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  Asymm_map_ = new TH2F("Asymm_map","Asymm_map", 82, x, 72, y);
  for(size_t i{0}; i < (1+sizeof(x)/sizeof(x[0]))*(1+sizeof(y)/sizeof(y[0])); ++i){ 
    Asymm_map_->SetBinContent(i, TH2_contents.at(i));  
  }

}

void HiInclusiveJetSubstructure::analyze(const Event& iEvent, const EventSetup& iSetup) {
  int event = iEvent.id().event();
  int run = iEvent.id().run();
  int lumi = iEvent.id().luminosityBlock();

  jets_.run = run;
  jets_.evt = event;
  jets_.lumi = lumi;
  // std::cout << ReadJetAsymmMap(1.5, -0.14, *Asymm_map_) << " expect 0.0704101 from map" << std::endl;
  LogDebug("HiInclusiveJetSubstructure")<<"START event: "<<event<<" in run "<<run<<endl;

    // loop the events
  reco::Vertex::Point vtx(0,0,0);
  if (useVtx_) {
    edm::Handle<vector<reco::Vertex> >vertex;
    iEvent.getByToken(vtxTag_, vertex);
    if(vertex->size()>0) {
      jets_.vx = vertex->begin()->x();
      jets_.vy = vertex->begin()->y();
      jets_.vz = vertex->begin()->z();
      vtx = vertex->begin()->position();
    }
  }

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetTag_, jets);

  edm::Handle<pat::JetCollection> matchedjets;
  iEvent.getByToken(matchTag_, matchedjets);

  auto const& pfCandidates = iEvent.get(pfCandidateToken_);
    // edm::Handle<reco::PFCandidateCollection> pfCandidates;
    // iEvent.getByToken(pfCandidateLabel_,pfCandidates);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackTag_,tracks);
    // 
  // edm::Handle<reco::GenParticleCollection> genparts;
  // if(isMC_){
  //   iEvent.getByToken(genParticleSrc_,genparts);
  // }

  edm::Handle<CaloTowerCollection> towers;
  if(doTower){
    iEvent.getByToken(TowerSrc_,towers);
  }

  double pthat = 0;
  if(isMC_){
    edm::Handle<GenEventInfoProduct> hEventInfo;
    iEvent.getByToken(eventGenInfoTag_,hEventInfo);
    pthat = hEventInfo->qScale();
  }

  jets_.nref = 0;
  // std::cout << "Checking jets inside the range " << jetAbsEtaMax_-rParam << " for radius parameter " << rParam << std::endl;
  // if(doChargedConstOnly_) std::cout << "Doing only charged jet constituents!" << std::endl;
  // else std::cout << "Doing ALL jet constituents!" << std::endl;
  jets_.triggerJetInAcceptance = false;
    // std::cout << "Number of jets: " << jets->size() << std::endl;
  for(unsigned int j = 0; j < jets->size(); ++j){
    const pat::Jet& jet = (*jets)[j];
    auto pt = useRawPt_ ? jet.correctedJet("Uncorrected").pt() : jet.pt();
    if(pt < jetPtMin_) continue;
    // if(std::abs(jet.eta()) > jetAbsEtaMax_-rParam) continue;
    // std::cout << "Raw pt: " << jet.correctedJet("Uncorrected").pt() << " corrected: " << jet.pt() << " eta: " << jet.eta() << std::endl;
    //assume highest jet in event is also the trigger object, check if it's within the eta acceptance above
    if(j==0){ 
      jets_.triggerJetInAcceptance = true;
    }

    jets_.rawpt[jets_.nref] = jet.correctedJet("Uncorrected").pt();
    jets_.jtrawE[jets_.nref] = jet.correctedJet("Uncorrected").energy();
    // jets_.jtMapPt[jets_.nref] = jet.pt()*(1+ReadJetAsymmMap(jet.eta(), jet.phi(), *Asymm_map_));
    jets_.jtpt[jets_.nref] = jet.pt();
    jets_.jteta[jets_.nref] = jet.eta();
    jets_.jtphi[jets_.nref] = jet.phi();
        // jets_.jtsym[jets_.nref] = 0;
        // jets_.jtrg[jets_.nref] = 0;
        // jets_.jtdynkt[jets_.nref] = 0;
    // jets_.jtangu[jets_.nref] = 0;

        // jets_.jtdyn_pt1[jets_.nref] = 0;
    // jets_.jtdyn_var[jets_.nref] = 0;
    jets_.jtdyn_split[jets_.nref] = 0;
    // jets_.jtdyn_theta[jets_.nref] = 0;
    jets_.jtdyn_kt[jets_.nref] = 0;
    jets_.jtdyn_z[jets_.nref] = 0;

    fastjet::PseudoJet *sub1Gen = new fastjet::PseudoJet();
    fastjet::PseudoJet *sub2Gen = new fastjet::PseudoJet();
    fastjet::PseudoJet *sub1Hyb = new fastjet::PseudoJet();
    fastjet::PseudoJet *sub2Hyb = new fastjet::PseudoJet();

    // std::cout << jets_.neutralN[jets_.nref] << " Neutral number before IDR " << jets_.neutralSum[jets_.nref] << std::endl;

    IterativeDeclusteringRec(groom_type, groom_combine, jet, sub1Hyb, sub2Hyb);
    // std::cout << jets_.neutralN[jets_.nref] << " Neutral number after IDR " << jets_.neutralSum[jets_.nref] << std::endl;

    jets_.refpt[jets_.nref] = 0;
    jets_.refeta[jets_.nref] = 0;
    jets_.refphi[jets_.nref] = 0;
    jets_.refsym[jets_.nref] = 0.;

        // jets_.refrg[jets_.nref] = 0;
        // jets_.refdynkt[jets_.nref] = 0;
    // jets_.refangu[jets_.nref] = 0; 
        // jets_.refdyn_pt1[jets_.nref] = 0;
    // jets_.refdyn_var[jets_.nref] = 0;
    jets_.refdyn_split[jets_.nref] = 0;
    // jets_.refdyn_theta[jets_.nref] = 0;
    jets_.refdyn_kt[jets_.nref] = 0;
    jets_.refdyn_z[jets_.nref] = 0;
    jets_.jtdyn_isClosestToTruth[jets_.nref] = 0;
    jets_.refdyn_isClosestToReco[jets_.nref] = 0;
    jets_.jtdyn_refdyn_dR[jets_.nref] = 0;


    jets_.refsub11[jets_.nref] = 0;
    jets_.refsub12[jets_.nref] = 0;
    jets_.refsub21[jets_.nref] = 0;
    jets_.refsub22[jets_.nref] = 0;
    // std::cout << jets_.jtJetConstituent.size() << " " << jets_.refJetConstituent.size() << " sizes of consts" << std::endl;
    if(isMC_){
      const reco::GenJet * genjet = jet.genJet();
      if(!genjet) continue;

      if(jet.genParton()){
        const reco::GenParticle & parton = *jet.genParton();
        jets_.refparton_pt[jets_.nref] = parton.pt();
        jets_.refparton_flavor[jets_.nref] = parton.pdgId();
      }
      else {
        jets_.refparton_pt[jets_.nref] = -999;
        jets_.refparton_flavor[jets_.nref] = -999;
      }

      jets_.refpt[jets_.nref] = genjet->pt();
      jets_.refeta[jets_.nref] = genjet->eta();
      jets_.refphi[jets_.nref] = genjet->phi();
            //cout<<"jet daughters gen"<<genjet->numberOfDaughters()<<endl;
      if(dopthatcut) if(pthat<0.35*genjet->pt()) continue;

      IterativeDeclusteringGen(groom_type, groom_combine, *genjet, sub1Gen, sub2Gen);
      if(doHardestSplitMatching_){
        TruthRecoRecoTruthMatching();
        // if(!(jets_.jtdyn_isClosestToTruth[jets_.nref] && jets_.refdyn_isClosestToReco[jets_.nref])){
        //   std::cout << "Failed matching " << jets_.jtdyn_isClosestToTruth[jets_.nref] << " " <<  jets_.refdyn_isClosestToReco[jets_.nref] << std::endl;
        //   std::cout << "dR=" << jets_.jtdyn_deltaR[jets_.nref] << " kt=" << jets_.jtdyn_kt[jets_.nref] << " truth jtpt=" << jets_.refpt[jets_.nref] << std::endl;

        //   if(jets_.jtdyn_deltaR[jets_.nref] > 0.175 && jets_.jtdyn_deltaR[jets_.nref] < 0.2 && jets_.jtdyn_kt[jets_.nref]>2 && jets_.jtdyn_kt[jets_.nref]<30 && jets_.jtpt[jets_.nref]>150 ){
        //     std::cout << "The above jet satisfies reco selection for purity" << std::endl;
        //     std::cout << "dR=" << jets_.jtdyn_deltaR[jets_.nref] << " kt=" << jets_.jtdyn_kt[jets_.nref] << " truth jtpt=" << jets_.refpt[jets_.nref] << std::endl;
        //   }
        // }
        // std::cout << "Constituent vector lengths " << jets_.jtJetConstituent.size() << " " << jets_.refJetConstituent.size() << std::endl;
        jets_.jtJetConstituent = {};
        jets_.refJetConstituent = {};
         // std::cout << "Constituent vector lengths after clear " << jets_.jtJetConstituent.size() << " " << jets_.refJetConstituent.size() << std::endl;
        
      }
      if(doSubjetPurity){
        jets_.refsub11[jets_.nref] = sqrt(pow((sub1Gen->rap()-sub1Hyb->rap()),2)+pow((sub1Gen->phi()-sub1Hyb->phi()),2));
        jets_.refsub12[jets_.nref] = sqrt(pow((sub1Gen->rap()-sub2Hyb->rap()),2)+pow((sub1Gen->phi()-sub2Hyb->phi()),2));
        jets_.refsub21[jets_.nref] = sqrt(pow((sub2Gen->rap()-sub1Hyb->rap()),2)+pow((sub2Gen->phi()-sub1Hyb->phi()),2));
        jets_.refsub22[jets_.nref] = sqrt(pow((sub2Gen->rap()-sub2Hyb->rap()),2)+pow((sub2Gen->phi()-sub2Hyb->phi()),2));
      }
    }

    delete sub1Gen;
    delete sub2Gen;
    delete sub1Hyb;
    delete sub2Hyb;       

    jets_.nref++;
  } 
  // jets_.triggerJetInAcceptance = trigger_jet_in_acceptance;
  t->Fill();
  memset(&jets_,0,sizeof jets_);

}
//reco::Jet& jet and const reco::GenJet& jet - can replace the two IterDec with one using a template? 
void HiInclusiveJetSubstructure::IterativeDeclusteringRec(double groom_type, double groom_combine, const reco::Jet& jet, fastjet::PseudoJet *sub1, fastjet::PseudoJet *sub2)
{
  std::cout << "new jet" << std::endl;
  Double_t map_Corrected_jtpt = jets_.rawpt[jets_.nref];
  Double_t assymDelta_jetpt = 0;
  //take the bin limits of the jet axis in the asymmetry map for neutrals 
  std::vector <float> jetBoundsInMap = {};
  if(doPeripheralNeuPFScaling_ || doCompensatoryNeuPFScaling_){
    jetBoundsInMap = BinBoundsAsymmMap(jet.eta(), jet.phi(), *Asymm_map_);
    map_Corrected_jtpt = jets_.rawpt[jets_.nref]*(1+ReadJetAsymmMap(jet.eta(), jet.phi(), *Asymm_map_));
    assymDelta_jetpt = map_Corrected_jtpt - jets_.rawpt[jets_.nref];
    std::cout << assymDelta_jetpt << " difference in jet pT" << std::endl;
  }

  // random number generators for different smearing scenarios (systematics)
  TRandom *rand_track_sel = new TRandom3(0);
  TRandom *rand_charge_smear = new TRandom3(0);
  TRandom *rand_hcal_y_smear = new TRandom3(0);
  TRandom *rand_hcal_phi_smear = new TRandom3(0);
  std::vector<Int_t> posDonor = {};
  std::vector<Int_t> posAcceptor = {};
  // generalised angularities
  Int_t intjet_multi = 0;
  float jet_girth = 0;
  float jet_girth_new = 0;
  float jet_thrust = 0;
  float jet_LHA = 0;
  float jet_pTD = 0;
  std::vector<float> jet_PLJPkT = {};
  std::vector<float> jet_PLJPdR = {};
  std::vector<float> jet_PLJPeta = {};
  std::vector<float> jet_PLJPphi = {};
  Int_t nsplit = 0;
  double dyn_kt = std::numeric_limits<double>::min();
  Int_t dyn_split = 0;
  double z = 0;
  double dyn_eta = 0;
  double dyn_phi = 0;
  double dyn_deltaR = 0;
  // double dyn_var = std::numeric_limits<double>::min();
  double dyn_z = 0;
  double jet_radius_ca = 1.0;
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);
  fastjet::PseudoJet myjet;
  fastjet::PseudoJet mypart;
  myjet.reset(jet.p4().px(),jet.p4().py(),jet.p4().pz(),jet.p4().e());
  // Reclustering jet constituents with new algorithm

  try{
    std::vector<fastjet::PseudoJet> particles = {};
    auto daughters = jet.getJetConstituents();
    Int_t count_it = 0;
    double JetNeutralEnergy = 0.;

    //scaling of neutral candidates using an (eta,phi) map from Mikko 
    // loop over the neutral candidates and scale the core ones only (within the jet axis HCal granularity)
    float assymDelta_corePtChange = 0;
    float assymSum_peripheralChange = 0;
    float assymSum_allNeuChange = 0;
    if(isMC_ && (doPeripheralNeuPFScaling_ or doCompensatoryNeuPFScaling_)){
      fastjet::PseudoJet coreSum;
      fastjet::PseudoJet coreSum_modified;
      fastjet::PseudoJet peripheralSum;
      fastjet::PseudoJet peripheralSum_modified;
      fastjet::PseudoJet totalSum;
      fastjet::PseudoJet totalSum_modified;

      coreSum.reset(0,0,0,0);
      coreSum_modified.reset(0,0,0,0);
      // std::cout << "before total modification calculation" << std::endl;
      for(auto it = daughters.begin(); it!=daughters.end(); ++it){
        Int_t abspdgId = abs((**it).pdgId());
        //if neutral
        if(abspdgId!=130) continue;
        // for(auto it2 = it; it2!=daughters.end(); ++it2){
        //   if(it2!=it && (**it2).pdgId()==130){
        //     fastjet::PseudoJet one;
        //     one.reset((**it).px(), (**it).py(), (**it).pz(), (**it).energy());
        //     fastjet::PseudoJet two;
        //     two.reset((**it2).px(), (**it2).py(), (**it2).pz(), (**it2).energy());
        //     float dRonetwo = one.delta_R(two);
        //     // if(dRonetwo<0.087){
        //       std::cout << "Below 0.087 between neutrals: " << dRonetwo << std::endl;
        //     // }
        //   }
        // }
        fastjet::PseudoJet pseudoNeutral;
        fastjet::PseudoJet pseudoNeutral_modified;
        float neutralMapValue = ReadJetAsymmMap((**it).eta(), (**it).phi(), *Asymm_map_);
        //if the modification is negative, set to 0.01 so we don't flip the momentum of the particle
        // if(neutralMapValue < 0.) neutralMapValue = 0.01;
        pseudoNeutral.reset((**it).px(), (**it).py(), (**it).pz(), (**it).energy());
        pseudoNeutral_modified.reset((**it).px()*(1.+neutralMapValue), (**it).py()*(1.+neutralMapValue), (**it).pz()*(1.+neutralMapValue), (**it).energy()*(1.+neutralMapValue));
        totalSum = totalSum + pseudoNeutral;
        totalSum_modified = totalSum_modified + pseudoNeutral_modified;
        //if particle is in the same tower as the jet axis, add to the core, if not, add to sum of peripherals
        //don't modify the particles here, do afterwards, just record the modification amplitude
        if((**it).eta() > jetBoundsInMap.at(0) && (**it).eta() > jetBoundsInMap.at(0) && (**it).eta() < jetBoundsInMap.at(1) && (**it).phi() > jetBoundsInMap.at(2) && (**it).phi() < jetBoundsInMap.at(3)){
          coreSum = coreSum + pseudoNeutral;
          coreSum_modified = coreSum_modified + pseudoNeutral_modified;
          // std::cout << "Core neutral " << std::endl;
        }
        else{ //if not in the same tower as jet axis, add the change in pT to sum of peripherals
          peripheralSum = peripheralSum + pseudoNeutral;
          peripheralSum_modified = peripheralSum_modified + pseudoNeutral_modified;
          // std::cout << "Peripheral neutral " << std::endl;
        }
      }

      // std::cout << "after total modification calculation" << std::endl;
      // totalSum = coreSum + peripheralSum;
      // totalSum_modified = coreSum_modified + peripheralSum_modified;
      assymDelta_corePtChange = coreSum_modified.pt() - coreSum.pt();
      assymSum_peripheralChange = peripheralSum_modified.pt() - peripheralSum.pt();
      // std::cout << "Energy of neutral changed 4 mom " << (totalSum_modified - totalSum).e() << " and pT=" << (totalSum_modified - totalSum).pt() << std::endl;
      assymSum_allNeuChange = (totalSum_modified - totalSum).pt();
      //negate total momentum change if the energy of the 4-vector difference is negative
      if((totalSum_modified - totalSum).e() < 0) assymSum_allNeuChange = - assymSum_allNeuChange;
    }
    //the change in the peripheral neutrals must compensate for the difference coming from core
    float peripheralNeutralChange = assymDelta_jetpt - assymDelta_corePtChange; 
    double jetEnergyDifferenceOutsideCore = 0;
    double JetNeutralEnergyFraction = JetNeutralEnergy/jets_.jtrawE[jets_.nref];
    // std::cout << jet.pt() << " " << jets_.jtpt[jets_.nref] << " jet pT " << jets_.jtrawE[jets_.nref] << " raw jet E" << std::endl;
    // std::cout << " ratio " << JetNeutralEnergyFraction << std::endl;
    // std::cout << "PFE shift below" << std::endl;
    for(auto it = daughters.begin(); it!=daughters.end(); ++it){
      //if we want only charged constituents and the daughter charge is 0, skip it
      //hide all this jetID stuff in a function (hooraaay?)
      if(doHiJetID_){
        incrementJetID(**it);
      }
      if(doChargedConstOnly_ && (**it).charge()==0) continue;
      double PFE_scale = 1.;
      double charge_track_smear = 1.;
      double Hcal_y_smear = 1.;
      double Hcal_phi_smear = 1.;
      //if it is MC, rescale the 4-momentum of the charged particles (we accept only them above) by pfCCES(+-1%)
      if(isMC_){
        if((**it).charge()!=0)
          PFE_scale = pfChargedCandidateEnergyScale_;
        else if((**it).pdgId()==22){
          PFE_scale = pfGammaCandidateEnergyScale_;
          // std::cout << "Photon found! " << (**it).mass() << std::endl;
        }
        else if((**it).pdgId()==130){
          PFE_scale = pfNeutralCandidateEnergyScale_;
          // if(doCustomNeuPFScaling_){
          //   double map_value = ReadJetAsymmMap((**it).eta(), (**it).phi(), *Asymm_map_);
          // }
          //treat all neutrals equally in order to compensate for jet asymmetry (no matter if in core or not)
          if(doCompensatoryNeuPFScaling_){
            double map_value = ReadJetAsymmMap((**it).eta(), (**it).phi(), *Asymm_map_);
            double particle_dp = map_value*(**it).pt();
            std::cout << "Fraction of excess momentum " << assymDelta_jetpt << " going into particle " << particle_dp/assymSum_allNeuChange << std::endl;
            std::cout << "Total sum in changed neutrals: " << assymSum_allNeuChange << std::endl;
            map_value = map_value*assymDelta_jetpt/assymSum_allNeuChange;
            // if(map_value < -1.) map_value = -0.99;
            PFE_scale = 1 + map_value;
            std::cout << "Neutral scaling value " << PFE_scale << " for particle of pT=" << (**it).pt() << std::endl;
          }
          if(doPeripheralNeuPFScaling_){
            double map_value = ReadJetAsymmMap((**it).eta(), (**it).phi(), *Asymm_map_);
            //factor of 10 due to jet core tower ratio to whole jet (Mikko suggested, but is this not effectively energy fraction as below?)
            map_value = 10*map_value*peripheralNeutralChange/assymSum_peripheralChange;
            // if(map_value < -1.) map_value = -0.99;
            // std::cout << map_value << " map value for neutral at eta-phi = " << (**it).eta() << " " << (**it).phi() << std::endl;
            PFE_scale = 1 + map_value;
          }
          if(doRatioNeuPFScaling_){
            // std::cout << "Distance of neutral to jet core:" << 
            double map_value = ReadJetAsymmMap((**it).eta(), (**it).phi(), *Asymm_map_);
            //factor of 10 due to jet core tower ratio to whole jet (Mikko suggested, but is this not effectively energy fraction as below?)
            map_value = map_value*10;
            // std::cout << map_value << " map value for neutral at eta-phi = " << (**it).eta() << " " << (**it).phi() << std::endl;
            PFE_scale = 1 + map_value/JetNeutralEnergyFraction;
          }
          else if(doNaiveNeuPFScaling_){
            double map_value = ReadJetAsymmMap((**it).eta(), (**it).phi(), *Asymm_map_);
            //factor of 10 due to jet core tower ratio to whole jet (Mikko suggested, but is this not effectively energy fraction as below?)
            map_value = map_value*10;
            PFE_scale = 1 + map_value;
          }
          // std::cout << PFE_scale << " Changed to " << std::endl;
          if(pfNeutralSmear_){
            while(abs(Hcal_y_smear) > 0.087){
              Hcal_y_smear = rand_hcal_y_smear->Gaus(0, 0.087/2);
            }
            while(abs(Hcal_phi_smear) > 0.087){
              Hcal_phi_smear = rand_hcal_phi_smear->Gaus(0, 0.087/2);
            }
          }
          // std::cout << " Smearing amount y: " << Hcal_y_smear << " phi: " << Hcal_phi_smear << std::endl;
          // std::cout << "Neutral found with rapidity: " << (**it).rap() << " and pseudorap: " << (**it).eta() << std::endl;
        }
        else{
          std::cout << "Not supposed to be here!!!!! What is this particle: " << (**it).pdgId() << std::endl;
        }
      }
      //vary tracking efficiency - drop TrackVariation_% of particles within the jet if using only charged particles
      if(isMC_ && TrackVariation_ != 0. && doChargedConstOnly_){
        // std::cout << "doing the track variation" << std::endl;
        if(rand_track_sel->Uniform(0,1) < TrackVariation_) continue;
      }
      //vary tracking efficiency - smear by 10% TrackVariation_% of charged particles within the jet if using inclusive
      else if(isMC_ && TrackVariation_ != 0. && !doChargedConstOnly_ && (**it).charge()!=0 && rand_track_sel->Uniform(0,1) < TrackVariation_ ){
        double closest_to_particle = 0.087;
        fastjet::PseudoJet temp_part;
        temp_part.reset((**it).px(), (**it).py(), (**it).pz(), (**it).energy());
        Bool_t acceptorCharge = false;
        // int acceptorPos = 0;
        Int_t  count_it2 = 0;
        Int_t selectedAcceptor = 0;
        for(auto it2 = daughters.begin(); it2!=daughters.end(); ++it2){
          fastjet::PseudoJet temp_part_acceptor;
          temp_part_acceptor.reset((**it2).px(), (**it2).py(), (**it2).pz(), (**it2).energy());
          Double_t dR_DonorAcceptor = temp_part_acceptor.delta_R(temp_part);
          // std::cout << dR_DonorAcceptor << " distance between charged marked for removal and acceptors " << count_it << " and " << count_it2 << " particle numbers" << std::endl;
          // if( dR_DonorAcceptor < 0.1 && count_it2!=count_it && (**it2).charge()==0 && (**it).charge()==0 && (**it2).pdgId()!=22 && (**it).pdgId()!=22  ){
            // std::cout << "Neutrals within distance " << dR_DonorAcceptor << " " << count_it2 << " " << count_it << std::endl;
          // }
          if( dR_DonorAcceptor < closest_to_particle && count_it2!=count_it && (**it2).pdgId()!=22 ){
            selectedAcceptor = count_it2;
            closest_to_particle = dR_DonorAcceptor;
            // std::cout << " Acceptor particle ID = " << (**it2).pdgId() << std::endl;
          }
          count_it2++;
        }
        if(closest_to_particle < 0.087){
          posDonor.push_back(count_it);
          posAcceptor.push_back(selectedAcceptor);
          // std::cout << "Will move particle " << count_it << " to particle " << selectedAcceptor << " due to dR =" << closest_to_particle << std::endl; 
        }
        double charge_track_shift = rand_charge_smear->Gaus(0, (**it).energy()*0.1);
        charge_track_smear = ((**it).energy()-charge_track_shift)/(**it).energy();
        // std::cout << charge_track_smear << " smear factor for particle with energy " << (**it).energy() << std::endl;
      }
      // std::cout << "Rescaling charged pfCand energy by " << PFE_scale << std::endl;
      mypart.reset((**it).px(), (**it).py(), (**it).pz(), (**it).energy());
      // using the jet tower term in addition to the neutral particle energy fraction scaling can give to negative factor for neutrals
      // safety to keep particles positive (negative scaling reverses px,py,pz)
      if(PFE_scale < 0.){
        // std::cout << "Neutral scale negative, reset " << PFE_scale << std::endl;
        PFE_scale = 0.01;
      }
      fastjet::PseudoJet copy_mp = mypart;
      if( isMC_ ) mypart.reset((**it).px()*PFE_scale*charge_track_smear, (**it).py()*PFE_scale*charge_track_smear, (**it).pz()*PFE_scale*charge_track_smear, (**it).energy()*PFE_scale*charge_track_smear);
      if( isMC_ && pfNeutralSmear_ && (**it).pdgId()==130 ){
        // std::cout << "Old y = " << mypart.rap() << " and eta =" << mypart.eta() << std::endl;
        mypart.reset_PtYPhiM(mypart.perp(), mypart.rap()+Hcal_y_smear, mypart.phi()+Hcal_phi_smear, mypart.m());
        // mypart.reset_PtYPhiM(p, y, phi, m);
        // std::cout << "new y and phi " << mypart.rap() << " and " << mypart.phi() << " perp and m " << mypart.perp() << " " << mypart.m() << std::endl;

      }
      std::cout << "Particle pT change " << mypart.pt() - copy_mp.pt() << " 4-vec diff " << (mypart-copy_mp).pt() << std::endl;

      // if(mypart.delta_R(myjet) < 0.1){
      //   std::cout << mypart.delta_R(myjet) << " Distance between particle and jet" << std::endl;
      //   std::cout << "Jet eta phi " << jet.p4().eta() << " " << jet.p4().phi() << std::endl;
      //   std::cout << "Eta phi for weird neutrals: " << mypart.pseudorapidity() << " " << mypart.phi() << " with ID " << (**it).pdgId() << std::endl;
      // }

      particles.push_back(mypart);
      double frac_dR = mypart.delta_R(myjet)/rParam;
      double frac_pt = mypart.perp()/myjet.perp();
      intjet_multi++;
      jet_girth     += mypart.perp()*mypart.delta_R(myjet)/myjet.perp();
      jet_girth_new += frac_pt*frac_dR;
      jet_thrust    += frac_pt*frac_dR*frac_dR;
      jet_LHA       += frac_pt*sqrt(frac_dR);
      jet_pTD       += frac_pt*frac_pt;
      count_it++;
    }
   // std::cout << " N/CH/C/P = " << jets_.neutralN[jets_.nref] << "/" << jets_.chargedN[jets_.nref] << "/" << jets_.eN[jets_.nref]+jets_.muN[jets_.nref] << "/" <<  jets_.photonN[jets_.nref] << std::endl;
    // std::cout << " total particles " << particles.size() << std::endl;
    // std::cout << "Particle container has " << particles.size() << " reco particles" << std::endl;
    // std::cout << "acceptors" << std::endl;
    // for(size_t acc{0}; acc<posAcceptor.size(); ++acc){
    //   std::cout << posAcceptor.at(acc) << " ";
    // }
    // std::cout << std::endl;
    // std::cout << "donors" << std::endl;
    // for(size_t don{0}; don<posDonor.size(); ++don){
    //   std::cout << posDonor.at(don) << " ";
    // }
    // std::cout << std::endl;
    // merge the particle donors/acceptors, zero and remove donors from list
    //if two charged happen to be close and both get dropped, they swap possitions and one still remains 0 vector
    if(isMC_ && TrackVariation_ != 0. && !doChargedConstOnly_){
      for(size_t pos{0}; pos < posDonor.size(); ++pos){
        // std::cout << "before merger pt=" << particles.at(posAcceptor.at(pos)).perp() << std::endl;
        particles.at(posAcceptor.at(pos)) = particles.at(posAcceptor.at(pos)) + particles.at(posDonor.at(pos));
        // std::cout << "after merger pt=" << particles.at(posAcceptor.at(pos)).perp() << std::endl;
        // std::cout << "donor before merger pt=" << particles.at(posDonor.at(pos)).perp() << std::endl;
        particles.at(posDonor.at(pos)).reset(0., 0., 0., 0.);
        // std::cout << "donor after merger pt=" << particles.at(posDonor.at(pos)).perp() << std::endl;
      }
    }
    // remove the particles just in case
    // for(size_t p{0}; p < particles.size(); ++p){
    //   std::cout << particles.at(p).perp() << " ";
    // }
    // std::cout << std::endl;
    Int_t itp_idx = 0;
    // std::cout << jets_.eN[jets_.nref] + jets_.photonN[jets_.nref] + jets_.muN[jets_.nref] + jets_.chargedN[jets_.nref] + jets_.neutralN[jets_.nref] << " constituents recorded out of " << particles.size() << std::endl;
    // std::cout << (jets_.eg_HFSum[jets_.nref] + jets_.h_HFSum[jets_.nref] + jets_.eSum[jets_.nref] + jets_.muSum[jets_.nref] + jets_.photonSum[jets_.nref] + jets_.chargedSum[jets_.nref] + jets_.neutralSum[jets_.nref])/jets_.jtrawE[jets_.nref] << " fractional energy" << std::endl;
    // std::cout << "something else above? jet eta " << jets_.jteta[jets_.nref] << std::endl;

    // loop again and remove all particles with pT=0 if TrackVariation is being used
    if(isMC_ and TrackVariation_){
      for(auto itp = particles.begin(); itp!=particles.end(); ++itp){
        // std::cout << itp_idx << " " << itp->perp() << " in checks" << std::endl;
        if(itp->perp()==0 and itp!=particles.end()){
          // std::cout << "removing " << itp_idx << " particle" << std::endl;
          particles.erase(itp);
          itp--;
        }
      // itp_idx++;
      }
    }
    
    // std::cout << "Final length of particle vector " << particles.size() << std::endl;
    // for(size_t p{0}; p < particles.size(); ++p){
    //   std::cout << particles.at(p).perp() << " ";
    // }
    // std::cout << std::endl;
    if(particles.empty()){
      // jets_.jtdyn_var[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jtdyn_eta[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_phi[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_girth_new[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_thrust[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_LHA[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_pTD[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_PLJPkT.push_back(jet_PLJPkT);
      jets_.jt_PLJPdR.push_back(jet_PLJPdR);
      jets_.jt_PLJPeta.push_back(jet_PLJPeta);
      jets_.jt_PLJPphi.push_back(jet_PLJPphi);
      throw(123);
    }

    fastjet::PseudoJet sumParticles;
    sumParticles.reset(0,0,0,0);
    for(size_t i{0}; i < particles.size(); ++i){
      sumParticles = sumParticles + particles.at(i);
    }
    std::cout << "Modified jet pT vs constituent pT " << map_Corrected_jtpt - sumParticles.pt() << " = " << map_Corrected_jtpt << " - " << sumParticles.pt() << " for jet of eta phi " << jet.eta() << " " << jet.phi() << std::endl;
    // std::cout << "Clustering " << particles.size() << " number of reco particles" << std::endl;
    fastjet::ClusterSequence csiter(particles, jet_def);
    std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
    // std::cout << output_jets.size() << " size of output jets" << std::endl;

    // for(size_t h{0}; h < output_jets.size(); ++h){
      // std::cout << output_jets.at(h).perp() << " pT of output jet " << h << std::endl;
    // }
    output_jets = sorted_by_pt(output_jets);
    fastjet::PseudoJet jj = output_jets[0];
    // LookThroughJetSplits(jj,1);
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    fastjet::PseudoJet highest_splitting;
    if(!jj.has_parents(j1,j2)){
      jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jtdyn_eta[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_phi[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_girth_new[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_thrust[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_LHA[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_pTD[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_PLJPkT.push_back(jet_PLJPkT);
      jets_.jt_PLJPdR.push_back(jet_PLJPdR);
      jets_.jt_PLJPeta.push_back(jet_PLJPeta);
      jets_.jt_PLJPphi.push_back(jet_PLJPphi);
      throw(124);
    }
    while(jj.has_parents(j1,j2)){
      if(j1.perp() < j2.perp()) std::swap(j1,j2);
      double delta_R = j1.delta_R(j2);
      if(doHardestSplitMatching_ && isMC_) jets_.jtJetConstituent.push_back(j2);
      double k_t = j2.perp()*delta_R;
      z = j2.perp()/(j1.perp()+j2.perp());
      // double dyn = z*(1-z)*j2.perp()*pow(delta_R/rParam,mydynktcut);
      // double dyn = 1./output_jets[0].perp()*z*(1-z)*jj.perp()*pow(delta_R/rParam,mydynktcut);
      // std::cout << "Reco split " << nsplit << " with k_T=" << k_t << " z=" << z << " eta " << j2.eta() << " phi " << j2.phi() <<  std::endl;
      if(doPrimaryLJPReco_){
        jet_PLJPkT.push_back(k_t);
        jet_PLJPdR.push_back(delta_R);
        jet_PLJPeta.push_back(j2.eta());
        jet_PLJPphi.push_back(j2.phi());
      }
      if(k_t > dyn_kt){
        // dyn_var = dyn;
        highest_splitting = j2;
        dyn_kt = k_t;
        dyn_split = nsplit;
        dyn_deltaR = delta_R;
        dyn_z = z;
        dyn_eta = j2.eta();
        dyn_phi = j2.phi();
      }
      jj = j1;
      nsplit = nsplit+1;
    }
    // std::cout << eSum/jet.pt() << " vs given " << jet.chargedEmEnergyFraction() << std::endl; 
    // std::cout << highest_splitting.eta() << " " << highest_splitting.phi() << " highest reco splitting eta phi at " << dyn_split << std::endl;
    // jets_.jtdyn_var[jets_.nref] = dyn_var;
    jets_.jtdyn_split[jets_.nref] = dyn_split;
    jets_.jtdyn_eta[jets_.nref] = dyn_eta;
    jets_.jtdyn_phi[jets_.nref] = dyn_phi;
    jets_.jtdyn_deltaR[jets_.nref] = dyn_deltaR;
    jets_.jtdyn_kt[jets_.nref] = dyn_kt;
    jets_.jtdyn_z[jets_.nref] = dyn_z;
    jets_.jt_intjet_multi[jets_.nref] = intjet_multi;
    jets_.jt_girth[jets_.nref] = jet_girth;
    jets_.jt_girth_new[jets_.nref] = jet_girth_new;
    jets_.jt_thrust[jets_.nref] = jet_thrust;
    jets_.jt_LHA[jets_.nref] = jet_LHA;
    jets_.jt_pTD[jets_.nref] = jet_pTD;
    jets_.jt_PLJPkT.push_back(jet_PLJPkT);
    jets_.jt_PLJPdR.push_back(jet_PLJPdR);
    jets_.jt_PLJPeta.push_back(jet_PLJPeta);
    jets_.jt_PLJPphi.push_back(jet_PLJPphi);
  } 
  catch (fastjet::Error) { /*return -1;*/ }
  catch (Int_t MyNum){
    if(MyNum == 123)
      std::cout << "Whoops, seems the number of charged jet constituents is 0! Setting all reco jet split variables to numeric min." << std::endl;
    if(MyNum == 124)
      std::cout << "Jet does not have any parents, out of the loop!" << std::endl;
  }
}

void HiInclusiveJetSubstructure::IterativeDeclusteringGen(double groom_type, double groom_combine, const reco::GenJet& jet,fastjet::PseudoJet *sub1,fastjet::PseudoJet *sub2)
{
  Int_t intjet_multi = 0;
  float jet_girth = 0;
  float jet_girth_new = 0;
  float jet_thrust = 0;
  float jet_LHA = 0;
  float jet_pTD = 0;
  std::vector<float> jet_PLJPkT = {};
  std::vector<float> jet_PLJPdR = {};
  std::vector<float> jet_PLJPeta = {};
  std::vector<float> jet_PLJPphi = {};
  double nsplit = 0;
  // double dyn_theta = 0;
  double dyn_kt = std::numeric_limits<double>::min();
  double dyn_eta = 0;
  double dyn_phi = 0;

  Int_t dyn_split = 0;
  double dyn_deltaR = 0;
  double z = 0;
  double dyn_z = 0;
  double jet_radius_ca = 1.0;
  fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);
  fastjet::PseudoJet myjet;
  fastjet::PseudoJet mypart;
  myjet.reset(jet.p4().px(),jet.p4().py(),jet.p4().pz(),jet.p4().e());  
    // Reclustering jet constituents with new algorithm
  try{
    std::vector<fastjet::PseudoJet> particles = {};                         
    auto daughters = jet.getJetConstituents();
    for(auto it = daughters.begin(); it!=daughters.end(); ++it){
      //if we want only charged constituents and the daughter charge is 0, skip it
      if(doChargedConstOnly_ && (**it).charge()==0) continue;
      particles.push_back(fastjet::PseudoJet((**it).px(), (**it).py(), (**it).pz(), (**it).energy()));
      mypart.reset((**it).px(), (**it).py(), (**it).pz(), (**it).energy());
      double frac_dR = mypart.delta_R(myjet)/rParam;
      double frac_pt = mypart.perp()/myjet.perp();
      intjet_multi++;
      jet_girth     += mypart.perp()*mypart.delta_R(myjet)/myjet.perp();
      jet_girth_new += frac_pt*frac_dR;
      jet_thrust    += frac_pt*frac_dR*frac_dR;
      jet_LHA       += frac_pt*sqrt(frac_dR);
      jet_pTD       += frac_pt*frac_pt;
    }
    // std::cout << "Particle container has " << particles.size() << " reco particles" << std::endl;
    if(particles.empty()){
      // jets_.jtdyn_var[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.refdyn_eta[jets_.nref] = -std::numeric_limits<double>::max();
      jets_.refdyn_phi[jets_.nref] = -std::numeric_limits<double>::max();
      jets_.refdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.ref_girth[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_girth_new[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_thrust[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_LHA[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_pTD[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_PLJPkT.push_back(jet_PLJPkT);
      jets_.ref_PLJPdR.push_back(jet_PLJPdR);
      jets_.ref_PLJPeta.push_back(jet_PLJPeta);
      jets_.ref_PLJPphi.push_back(jet_PLJPphi);
      throw(123);
    }
    // std::cout << "Clustering " << particles.size() << " number of truth particles" << std::endl;
    fastjet::ClusterSequence csiter(particles, jet_def);
    std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
    output_jets = sorted_by_pt(output_jets);

    fastjet::PseudoJet jj = output_jets[0];
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    fastjet::PseudoJet highest_splitting;
    if(!jj.has_parents(j1,j2)){
      jets_.refdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.refdyn_eta[jets_.nref] = -std::numeric_limits<double>::max();
      jets_.refdyn_phi[jets_.nref] = -std::numeric_limits<double>::max();
      jets_.refdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.ref_girth[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_girth_new[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_thrust[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_LHA[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_pTD[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_PLJPkT.push_back(jet_PLJPkT);
      jets_.ref_PLJPdR.push_back(jet_PLJPdR);
      jets_.ref_PLJPeta.push_back(jet_PLJPeta);
      jets_.ref_PLJPphi.push_back(jet_PLJPphi);
      throw(124);
    }
    while(jj.has_parents(j1,j2)){
      if(j1.perp() < j2.perp()) std::swap(j1,j2);
      double delta_R = j1.delta_R(j2);
      if(doHardestSplitMatching_ && isMC_) jets_.refJetConstituent.push_back(j2);
      double k_t = j2.perp()*delta_R;
      z = j2.perp()/(j1.perp()+j2.perp());
      if(doPrimaryLJPTruth_){
        jet_PLJPkT.push_back(k_t);
        jet_PLJPdR.push_back(delta_R);
        jet_PLJPeta.push_back(j2.eta());
        jet_PLJPphi.push_back(j2.phi());
      }
      
      // std::cout << "Truth split " << nsplit << " with k_T=" << k_t << " z=" << z <<  " eta " << j2.eta() << " phi " << j2.phi() <<  std::endl;
      if(k_t > dyn_kt){
        highest_splitting = j2;
        dyn_kt = k_t;
        dyn_split = nsplit;
        dyn_deltaR = delta_R;
        dyn_z = z;
        dyn_eta = j2.eta();
        dyn_phi = j2.phi();
      }
      jj = j1;
      nsplit = nsplit+1;
    }
    // std::cout << highest_splitting.eta() << " " << highest_splitting.phi() << " highest truth splitting eta phi at " << dyn_split << std::endl;
    jets_.refdyn_split[jets_.nref] = dyn_split;
    jets_.refdyn_eta[jets_.nref] = dyn_eta;
    jets_.refdyn_phi[jets_.nref] = dyn_phi;
    jets_.refdyn_deltaR[jets_.nref] = dyn_deltaR;
    jets_.refdyn_kt[jets_.nref] = dyn_kt;
    jets_.refdyn_z[jets_.nref] = dyn_z;
    jets_.ref_intjet_multi[jets_.nref] = intjet_multi;
    jets_.ref_girth[jets_.nref] = jet_girth;
    jets_.ref_girth_new[jets_.nref] = jet_girth_new;
    jets_.ref_thrust[jets_.nref] = jet_thrust;
    jets_.ref_LHA[jets_.nref] = jet_LHA;
    jets_.ref_pTD[jets_.nref] = jet_pTD;
    jets_.ref_PLJPkT.push_back(jet_PLJPkT);
    jets_.ref_PLJPdR.push_back(jet_PLJPdR);
    jets_.ref_PLJPeta.push_back(jet_PLJPeta);
    jets_.ref_PLJPphi.push_back(jet_PLJPphi);
  } 
  catch (fastjet::Error) { /*return -1;*/ }
  catch (Int_t MyNum){
    if(MyNum == 123)
      std::cout << "Whoops, seems the number of charged jet constituents is 0! Setting all gen jet split variables to numeric min." << std::endl;
    if(MyNum == 124)
      std::cout << "Jet does not have any parents, out of the loop!" << std::endl;
  }
}
// My dream is indefinitely on hold :(
// template<typename T>
// void HiInclusiveJetSubstructure::TemplateDeclustering(Bool_t recoLevel, double groom_type, T& jet)
// {
//   TRandom *rand_track_sel = new TRandom3(0);
//   Int_t intjet_multi = 0;
//   float jet_girth = 0;
//   Int_t nsplit = 0;
//   double dyn_kt = std::numeric_limits<double>::min();
//   Int_t dyn_split = 0;
//   double z = 0;
//   double dyn_deltaR = 0;
//   // double dyn_var = std::numeric_limits<double>::min();
//   double dyn_z = 0;
//   fastjet::PseudoJet myjet;
//   fastjet::PseudoJet mypart;
//   myjet.reset(jet.p4().px(),jet.p4().py(),jet.p4().pz(),jet.p4().e());
//   // Reclustering jet constituents with new algorithm
//   try{
//     std::vector<fastjet::PseudoJet> particles = {};                         
//     auto daughters = jet.getJetConstituents();
//         // std::cout << "Number of pfCand " << pfCandidates.size() << std::endl;
//     // std::vector<int> vec_jet_consituent_charge;
//     for(auto it = daughters.begin(); it!=daughters.end(); ++it){
//       //if we want only charged constituents and the daughter charge is 0, skip it
//       if(doChargedConstOnly_ && (**it).charge()==0) continue;
//       double PFE_scale = 1.;
//       //if it is MC, rescale the 4-momentum of the charged particles (we accept only them above) by pfCCES(+-1%)
//       if(isMC_ && recoLevel && (**it).charge()!=0) PFE_scale = pfChargedCandidateEnergyScale_;
//       //shouldn't go into else in charged MC case 
//       else if(isMC_ && recoLevel && (**it).charge()==0) PFE_scale = pfNeutralCandidateEnergyScale_;
//       //vary tracking efficiency - drop ~4% of particles within the jet
//       if(isMC_ && recoLevel && doTrackVariation_){
//         // std::cout << "doing the track variation" << std::endl;
//         if(rand_track_sel->Uniform(0,1)<0.05) continue;
//       }
//       // std::cout << "Rescaling charged pfCand energy by " << PFE_scale << std::endl;
//       particles.push_back(fastjet::PseudoJet((**it).px()*PFE_scale, (**it).py()*PFE_scale, (**it).pz()*PFE_scale, (**it).energy()*PFE_scale));
//       mypart.reset((**it).px()*PFE_scale, (**it).py()*PFE_scale, (**it).pz()*PFE_scale, (**it).energy()*PFE_scale);
//       intjet_multi++;
//       jet_girth += mypart.perp()*mypart.delta_R(myjet)/myjet.perp();

//     }
//     //Call substructure loop function here

//     // std::cout << "Particle container has " << particles.size() << " reco particles" << std::endl;
//     if(particles.empty()){
//       // jets_.jtdyn_var[jets_.nref] = - std::numeric_limits<double>::max();
//       jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
//       jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
//       jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
//       jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
//       jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
//       jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
//       throw(123);
//     }
//     // std::cout << "Clustering " << particles.size() << " number of reco particles" << std::endl;
//     fastjet::ClusterSequence csiter(particles, jet_def);
//     std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
//     output_jets = sorted_by_pt(output_jets);
//     fastjet::PseudoJet jj = output_jets[0];
//     fastjet::PseudoJet j1;
//     fastjet::PseudoJet j2;
//     fastjet::PseudoJet highest_splitting;
//     if(!jj.has_parents(j1,j2)){
//       jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
//       jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
//       jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
//       jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
//       jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
//       jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
//       throw(124);
//     }
//     while(jj.has_parents(j1,j2)){
//       if(j1.perp() < j2.perp()) std::swap(j1,j2);
//       double delta_R = j1.delta_R(j2);
//       if(doHardestSplitMatching_ && isMC_) jets_.jtJetConstituent.push_back(j2);
//       double k_t = j2.perp()*delta_R;
//       z = j2.perp()/(j1.perp()+j2.perp());
//       // double dyn = z*(1-z)*j2.perp()*pow(delta_R/rParam,mydynktcut);
//       // double dyn = 1./output_jets[0].perp()*z*(1-z)*jj.perp()*pow(delta_R/rParam,mydynktcut);
//       // std::cout << "Reco split " << nsplit << " with k_T=" << k_t << " z=" << z << " eta " << j2.eta() << " phi " << j2.phi() <<  std::endl;
//       if(k_t > dyn_kt){
//         // dyn_var = dyn;
//         highest_splitting = j2;
//         dyn_kt = k_t;
//         dyn_split = nsplit;
//         dyn_deltaR = delta_R;
//         dyn_z = z;
//       }
//       jj = j1;
//       nsplit = nsplit+1;
//     }

//     // std::cout << highest_splitting.eta() << " " << highest_splitting.phi() << " highest reco splitting eta phi at " << dyn_split << std::endl;
//     // jets_.jtdyn_var[jets_.nref] = dyn_var;
//     jets_.jtdyn_split[jets_.nref] = dyn_split;
//     jets_.jtdyn_deltaR[jets_.nref] = dyn_deltaR;
//     jets_.jtdyn_kt[jets_.nref] = dyn_kt;
//     jets_.jtdyn_z[jets_.nref] = dyn_z;
//     jets_.jt_intjet_multi[jets_.nref] = intjet_multi;
//     jets_.jt_girth[jets_.nref] = jet_girth;
//   } 
//   catch (fastjet::Error) { /*return -1;*/ }
//   catch (Int_t MyNum){
//     if(MyNum == 123)
//       std::cout << "Whoops, seems the number of charged jet constituents is 0! Setting all reco jet split variables to numeric min." << std::endl;
//     if(MyNum == 124)
//       std::cout << "Jet does not have any parents, out of the loop!" << std::endl;
//   }
// }

// void HiInclusiveJetSubstructure::ClusterConstituents(std::vector<fastjet::PseudoJet> particles, Bool_t recoLevel){
//   double jet_radius_ca = 1.0;
//   fastjet::JetDefinition jet_def(fastjet::genkt_algorithm,jet_radius_ca,0,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);

// }

//maybe there is a more elegant way than the one below for matching...
void HiInclusiveJetSubstructure::RecoTruthSplitMatching(std::vector<fastjet::PseudoJet> &constituents_level1, fastjet::PseudoJet &hardest_level2, bool *bool_array, int *hardest_level1_split){
    //for now only include geometric matching, maybe consider pt/z, Lund plane location, etc...
  float min_dR = std::numeric_limits<float>::max();
  size_t closest_level1 = 0;
  // std::cout << "Starting loop over " << constituents_level1.size() << " particles" << std::endl;
  for(size_t i{0};i<constituents_level1.size();++i){
    float dR = constituents_level1.at(i).delta_R(hardest_level2);
    if(min_dR > dR){
      closest_level1 = i;
      min_dR = dR;
    }
  }
  // std::cout << "Compare particle " << static_cast<int>(closest_level1) << " with hardest split " << hardest_level1_split[jets_.nref] << std::endl;
  if(static_cast<int>(closest_level1) == hardest_level1_split[jets_.nref] ){
    bool_array[jets_.nref] = true;
  }
  else{
    // std::cout << "Sorry, closest pair is " << min_dR << " away with index " << static_cast<int>(closest_level1) << " as opposed to " << hardest_level1_split[jets_.nref] << std::endl;
    bool_array[jets_.nref] = false;
  }
}

float HiInclusiveJetSubstructure::ReadJetAsymmMap(float eta, float phi, TH2F Asymm_map){
  float result =  Asymm_map.GetBinContent(Asymm_map.GetXaxis()->FindBin(eta),Asymm_map.GetYaxis()->FindBin(phi));
  float high_bound = 0.15;
  float low_bound = -0.15;
  if( result < low_bound ){
    result = low_bound;
  }
  if( result > high_bound ){
    result = high_bound;
  }
  // std::cout << "For particle with eta-phi " << eta << " " << phi << " " << "result of " << result << std::endl;
  return result;
}

std::vector<float> HiInclusiveJetSubstructure::BinBoundsAsymmMap(float eta, float phi, TH2F Asymm_map){
  std::vector<float> result;
  float eta_low =  Asymm_map.GetXaxis()->GetBinLowEdge(Asymm_map.GetXaxis()->FindBin(eta));
  float eta_high =  Asymm_map.GetXaxis()->GetBinLowEdge(Asymm_map.GetXaxis()->FindBin(eta)+1);
  float phi_low = Asymm_map.GetYaxis()->GetBinLowEdge(Asymm_map.GetYaxis()->FindBin(phi));
  float phi_high = Asymm_map.GetYaxis()->GetBinLowEdge(Asymm_map.GetYaxis()->FindBin(phi)+1);
  return result = {eta_low, eta_high, phi_low, phi_high};
}

void HiInclusiveJetSubstructure::TruthRecoRecoTruthMatching(){
  // std::cout << jets_.jtdyn_split[jets_.nref] << " " << jets_.refdyn_split[jets_.nref] << " numbers of highest splits" << std::endl;
  if( jets_.jtdyn_split[jets_.nref] == std::numeric_limits<int>::min() || jets_.refdyn_split[jets_.nref] == std::numeric_limits<int>::min() || jets_.jtJetConstituent.size() == 0 || jets_.refJetConstituent.size() == 0 ){
    jets_.refdyn_isClosestToReco[jets_.nref] = false;
    jets_.jtdyn_isClosestToTruth[jets_.nref] = false;
    jets_.jtdyn_refdyn_dR[jets_.nref] = std::numeric_limits<float>::max();
    return;
  }
  //mind how the split number is defined in the reclustering
  fastjet::PseudoJet hardest_R_split = jets_.jtJetConstituent.at(jets_.jtdyn_split[jets_.nref]);
  fastjet::PseudoJet hardest_T_split = jets_.refJetConstituent.at(jets_.refdyn_split[jets_.nref]);
  // std::cout << hardest_R_split.eta() << " " << hardest_R_split.phi() << " hardest reco  splitting in matching" << std::endl;
  // std::cout << hardest_T_split.eta() << " " << hardest_T_split.phi() << " hardest truth splitting in matching" << std::endl;
  // std::cout << "Angle between hardest splits is dR = " << hardest_R_split.delta_R(hardest_T_split) << std::endl;
  jets_.jtdyn_refdyn_dR[jets_.nref] = hardest_R_split.delta_R(hardest_T_split);
  // std::cout << "truth loop" << std::endl;
  RecoTruthSplitMatching(jets_.refJetConstituent, hardest_R_split, jets_.refdyn_isClosestToReco, jets_.refdyn_split);
  // std::cout << "reco loop" << std::endl;
  RecoTruthSplitMatching(jets_.jtJetConstituent,  hardest_T_split, jets_.jtdyn_isClosestToTruth, jets_.jtdyn_split);
}

int HiInclusiveJetSubstructure::getPFJetMuon(const pat::Jet& pfJet, const reco::PFCandidateCollection *pfCandidateColl)
{

  int pfMuonIndex = -1;
  float ptMax = 0.;

  for(unsigned icand=0;icand<pfCandidateColl->size(); icand++) {
    const reco::PFCandidate& pfCandidate = pfCandidateColl->at(icand);
    int id = pfCandidate.particleId();
    if(abs(id) != 3) continue;
    if(reco::deltaR(pfJet,pfCandidate)>0.5) continue;

    double pt =  pfCandidate.pt();
    if(pt>ptMax){
      ptMax = pt;
      pfMuonIndex = (int) icand;
    }
  }

  return pfMuonIndex;
}

void HiInclusiveJetSubstructure::LookThroughJetSplits(fastjet::PseudoJet jj, int i=1){
  fastjet::PseudoJet j1;
  fastjet::PseudoJet j2;
  if(i==1)
    std::cout << "primary" << std::endl;
  else
    std::cout << "secondary" << std::endl;

  if(jj.has_parents(j1,j2)){
    // jj.has_parents(j1,j2);
    // std::cout << "here" << std::endl;

    LookThroughJetSplits(j2, 2);
    LookThroughJetSplits(j1, 1);
  }
}


double HiInclusiveJetSubstructure::getPtRel(const reco::PFCandidate& lep, const pat::Jet& jet )
{
  float lj_x = jet.p4().px();
  float lj_y = jet.p4().py();
  float lj_z = jet.p4().pz();

    // absolute values squared
  float lj2  = lj_x*lj_x+lj_y*lj_y+lj_z*lj_z;
  float lep2 = lep.px()*lep.px()+lep.py()*lep.py()+lep.pz()*lep.pz();

    // projection vec(mu) to lepjet axis
  float lepXlj = lep.px()*lj_x+lep.py()*lj_y+lep.pz()*lj_z;

    // absolute value squared and normalized
  float pLrel2 = lepXlj*lepXlj/lj2;

    // lep2 = pTrel2 + pLrel2
  float pTrel2 = lep2-pLrel2;

  return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
}

// Recursive function, but this version gets called only the first time

void HiInclusiveJetSubstructure::saveDaughters(const reco::GenParticle &gen){

  for(unsigned i = 0; i<gen.numberOfDaughters(); i++){
    const reco::Candidate & daughter = *gen.daughter(i);
    double daughterPt = daughter.pt();
    if(daughterPt<1.) continue;
    double daughterEta = daughter.eta();
    if(fabs(daughterEta)>3.) continue;
    int daughterPdgId = daughter.pdgId();
    int daughterStatus = daughter.status();
        // Special case when b->b+string, both b and string contain all daughters, so only take the string
    if(gen.pdgId()==daughterPdgId && gen.status()==3 && daughterStatus==2) continue;

        // cheesy way of finding strings which were already used
    if(daughter.pdgId()==92){
      for(unsigned ist = 0;ist<usedStringPts.size();ist++){
       if(fabs(daughter.pt() - usedStringPts[ist]) < 0.0001) return;
     }
     usedStringPts.push_back(daughter.pt());
   }
   jets_.bJetIndex[jets_.bMult] = jets_.nref;
   jets_.bStatus[jets_.bMult] = daughterStatus;
   jets_.bVx[jets_.bMult] = daughter.vx();
   jets_.bVy[jets_.bMult] = daughter.vy();
   jets_.bVz[jets_.bMult] = daughter.vz();
   jets_.bPt[jets_.bMult] = daughterPt;
   jets_.bEta[jets_.bMult] = daughterEta;
   jets_.bPhi[jets_.bMult] = daughter.phi();
   jets_.bPdg[jets_.bMult] = daughterPdgId;
   jets_.bChg[jets_.bMult] = daughter.charge();
   jets_.bMult++;
   saveDaughters(daughter);
 }
}

// This version called for all subsequent calls
void HiInclusiveJetSubstructure::saveDaughters(const reco::Candidate &gen){

  for(unsigned i = 0; i<gen.numberOfDaughters(); i++){
    const reco::Candidate & daughter = *gen.daughter(i);
    double daughterPt = daughter.pt();
    if(daughterPt<1.) continue;
    double daughterEta = daughter.eta();
    if(fabs(daughterEta)>3.) continue;
    int daughterPdgId = daughter.pdgId();
    int daughterStatus = daughter.status();
        // Special case when b->b+string, both b and string contain all daughters, so only take the string
    if(gen.pdgId()==daughterPdgId && gen.status()==3 && daughterStatus==2) continue;

        // cheesy way of finding strings which were already used
    if(daughter.pdgId()==92){
      for(unsigned ist=0;ist<usedStringPts.size();ist++){
        if(fabs(daughter.pt() - usedStringPts[ist]) < 0.0001) return;
      }
      usedStringPts.push_back(daughter.pt());
    }

    jets_.bJetIndex[jets_.bMult] = jets_.nref;
    jets_.bStatus[jets_.bMult] = daughterStatus;
    jets_.bVx[jets_.bMult] = daughter.vx();
    jets_.bVy[jets_.bMult] = daughter.vy();
    jets_.bVz[jets_.bMult] = daughter.vz();
    jets_.bPt[jets_.bMult] = daughterPt;
    jets_.bEta[jets_.bMult] = daughterEta;
    jets_.bPhi[jets_.bMult] = daughter.phi();
    jets_.bPdg[jets_.bMult] = daughterPdgId;
    jets_.bChg[jets_.bMult] = daughter.charge();
    jets_.bMult++;
    saveDaughters(daughter);
  }
}

//--------------------------------------------------------------------------------------------------
void HiInclusiveJetSubstructure::analyzeSubjets(const reco::Jet& jet) {

  std::vector<float> sjpt;
  std::vector<float> sjeta;
  std::vector<float> sjphi;
  std::vector<float> sjm;
  if(jet.numberOfDaughters()>0) {
    for (unsigned k = 0; k < jet.numberOfDaughters(); ++k) {
      const reco::Candidate & dp = *jet.daughter(k);
      sjpt.push_back(dp.pt());
      sjeta.push_back(dp.eta());
      sjphi.push_back(dp.phi());
      sjm.push_back(dp.mass());
    }
  } 
  else {
    sjpt.push_back(-999.);
    sjeta.push_back(-999.);
    sjphi.push_back(-999.);
    sjm.push_back(-999.);
  }
  jets_.jtSubJetPt.push_back(sjpt);
  jets_.jtSubJetEta.push_back(sjeta);
  jets_.jtSubJetPhi.push_back(sjphi);
  jets_.jtSubJetM.push_back(sjm);  
}

//--------------------------------------------------------------------------------------------------
int HiInclusiveJetSubstructure::getGroomedGenJetIndex(const reco::GenJet& jet) const {

    //Find closest soft-dropped gen jet
  double drMin = 100;
  int imatch = -1;
  for(unsigned int i = 0 ; i < gensubjets_->size(); ++i) {
    const reco::Jet& mjet = (*gensubjets_)[i];
    double dr = deltaR(jet,mjet);
    if(dr < drMin){
      imatch = i;
      drMin = dr;
    }
  }
  return imatch;
}

//--------------------------------------------------------------------------------------------------
void HiInclusiveJetSubstructure::analyzeRefSubjets(const reco::GenJet& jet) {

    //Find closest soft-dropped gen jet
  int imatch = getGroomedGenJetIndex(jet);

}

//--------------------------------------------------------------------------------------------------
void HiInclusiveJetSubstructure::analyzeGenSubjets(const reco::GenJet& jet) {
  //Find closest soft-dropped gen jet
  int imatch = getGroomedGenJetIndex(jet);
  double dr = 999.;

}

void HiInclusiveJetSubstructure::incrementJetID(const reco::Candidate& it){
  if(it.charge()!=0){
    // reco::Track const& track = (**it).pseudoTrack();
    TrackRef trk = it.get<TrackRef>();
    if(trk.isNonnull()){
    // if(!useQuality_ || (useQuality_ && reco::TrackBase::qualityByName(trackQuality_))){
      double ptcand = trk->pt();
      jets_.trackSum[jets_.nref] += ptcand;
      jets_.trackN[jets_.nref] += 1;

      if (ptcand > hardPtMin_) {
        jets_.trackHardSum[jets_.nref] += ptcand;
        jets_.trackHardN[jets_.nref] += 1;
      }
      if (ptcand > jets_.trackMax[jets_.nref])
        jets_.trackMax[jets_.nref] = ptcand;
    }
  }
  Double_t ecand = it.energy();
  // std::cout << (**it).pt() << " " << (**it).energy() << " pt and E" << std::endl;
  Int_t abspdgId = abs(it.pdgId());
  // if(abspdgId < 0) std::cout << "Negative pdg ID, what do?? " << (**it).pdgId() << std::endl;
  // if photon
  if( abspdgId == 22 ){
    jets_.photonSum[jets_.nref] += ecand;
    jets_.photonN[jets_.nref] += 1;
    if (ecand > hardPtMin_) {
      jets_.photonHardSum[jets_.nref] += ecand;
      jets_.photonHardN[jets_.nref] += 1;
    }
    if (ecand > jets_.photonMax[jets_.nref])
      jets_.photonMax[jets_.nref] = ecand;
  }
  // if neutral hadron
  else if( abspdgId == 130 ){
    jets_.neutralSum[jets_.nref] += ecand;
    jets_.neutralN[jets_.nref] += 1;
    if (ecand > jets_.neutralMax[jets_.nref])
      jets_.neutralMax[jets_.nref] = ecand;
  }
  // if charged hadron
  else if( abspdgId == 211 ){
    jets_.chargedSum[jets_.nref] += ecand;
    jets_.chargedN[jets_.nref] += 1;
    if (ecand > hardPtMin_) {
      jets_.chargedHardSum[jets_.nref] += ecand;
      jets_.chargedHardN[jets_.nref] += 1;
    }
    if (ecand > jets_.chargedMax[jets_.nref])
      jets_.chargedMax[jets_.nref] = ecand;
  }
  // electron
  else if( abspdgId == 11 ){
    jets_.eSum[jets_.nref] += ecand;
    jets_.eN[jets_.nref] += 1;
    if (ecand > jets_.eMax[jets_.nref])
      jets_.eMax[jets_.nref] = ecand;
  }
  // muon
  else if( abspdgId == 13 ){
    jets_.muSum[jets_.nref] += ecand;
    jets_.muN[jets_.nref] += 1;
    if (ecand > jets_.muMax[jets_.nref])
      jets_.muMax[jets_.nref] = ecand;
  }
  // hadronic forward
  else if( abspdgId == 1 ){
    jets_.h_HFSum[jets_.nref] += ecand;
    jets_.h_HFN[jets_.nref] += 1;
    if (ecand > jets_.h_HFMax[jets_.nref])
      jets_.h_HFMax[jets_.nref] = ecand;
  }
  // EG forward
  else if( abspdgId == 2 ){
    jets_.eg_HFSum[jets_.nref] += ecand;
    jets_.eg_HFN[jets_.nref] += 1;
    if (ecand > jets_.eg_HFMax[jets_.nref])
      jets_.eg_HFMax[jets_.nref] = ecand;
  }
  // what if else??
  else{
    std::cout << " Unknown particle! something else??? ID: " << abspdgId << std::endl;
  }
}

//--------------------------------------------------------------------------------------------------
// float HiInclusiveJetSubstructure::getAboveCharmThresh(TrackRefVector& selTracks, const TrackIPTagInfo& ipData, int sigOrVal)
// {

//     const double pdgCharmMass = 1.290;
//     btag::SortCriteria sc;
//     switch(sigOrVal){
//         case 1: //2d significance 
//         sc = reco::btag::IP2DSig;
//         case 2: //3d significance
//         sc = reco::btag::IP3DSig;
//         case 3:
//         sc = reco::btag::IP2DSig;  //values are not sortable!
//         case 4:
//         sc = reco::btag::IP3DSig;
// 	}
//     std::vector<std::size_t> indices = ipData.sortedIndexes(sc);
//     reco::TrackKinematics kin;
//     for(unsigned int i=0; i<indices.size(); i++){
//         size_t idx = indices[i];
//         const Track track = *(selTracks[idx]);
//         const btag::TrackIPData &data = ipData.impactParameterData()[idx];
//         kin.add(track);
//         if(kin.vectorSum().M() > pdgCharmMass){
//             switch(sigOrVal){
//                 case 1:
//                 return data.ip2d.significance();
//                 case 2:
//                 return data.ip3d.significance();
//                 case 3:
//                 return data.ip2d.value();
//                 case 4:
//                 return data.ip3d.value();
//             }
//         }
//     }
// 	return 0;
// }

DEFINE_FWK_MODULE(HiInclusiveJetSubstructure);
