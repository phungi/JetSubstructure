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
  if(isMC_){
    pfChargedCandidateEnergyScale_ = iConfig.getUntrackedParameter<double>("pfChargedEnergyScaleVar",1.);
    doTrackVariation_ = iConfig.getUntrackedParameter<bool>("doTrackVariation",false);
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
  std::cout << "Track efficiency var: " << doTrackVariation_ <<std::endl;
  std::cout << "Charged pfCand 4-mom var: " << pfChargedCandidateEnergyScale_ << std::endl;
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
  t->Branch("jtpt",jets_.jtpt,"jtpt[nref]/F");
  t->Branch("jteta",jets_.jteta,"jteta[nref]/F");
  t->Branch("jtphi",jets_.jtphi,"jtphi[nref]/F");
  // t->Branch("jtsym",jets_.jtsym,"jtsym[nref]/F");
  // t->Branch("jtrg",jets_.jtrg,"jtrg[nref]/F");
  // t->Branch("jtdynkt",jets_.jtdynkt,"jtdynkt[nref]/F");
  // t->Branch("jtdyn_pt1",jets_.jtdyn_pt1,"jtdyn_pt1[nref]/F");
  // t->Branch("jtdyn_var",jets_.jtdyn_var,"jtdyn_var[nref]/F");
  t->Branch("jtdyn_split",jets_.jtdyn_split,"jtdyn_split[nref]/I");
  // t->Branch("jtdyn_theta",jets_.jtdyn_theta,"jtdyn_theta[nref]/F");
  t->Branch("jtdyn_deltaR",jets_.jtdyn_deltaR,"jtdyn_deltaR[nref]/F");
  t->Branch("jtdyn_kt",jets_.jtdyn_kt,"jtdyn_kt[nref]/F");
  t->Branch("jtdyn_z",jets_.jtdyn_z,"jtdyn_z[nref]/F");
  t->Branch("jt_intjet_multi", jets_.jt_intjet_multi,"jt_intjet_multi[nref]/I");
  t->Branch("jt_girth", jets_.jt_girth, "jt_girth[nref]/F");
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
    // t->Branch("refdyn_theta",jets_.refdyn_theta,"refdyn_theta[nref]/F");
    t->Branch("refdyn_deltaR",jets_.refdyn_deltaR,"refdyn_deltaR[nref]/F");
    t->Branch("refdyn_kt",jets_.refdyn_kt,"refdyn_kt[nref]/F");
    t->Branch("refdyn_z",jets_.refdyn_z,"refdyn_z[nref]/F");
    t->Branch("ref_intjet_multi", jets_.ref_intjet_multi,"ref_intjet_multi[nref]/I");
    t->Branch("ref_girth", jets_.ref_girth, "ref_girth[nref]/F");
    t->Branch("jtdyn_isClosestToTruth", jets_.jtdyn_isClosestToTruth,"jtdyn_isClosestToTruth[nref]/O");
    t->Branch("refdyn_isClosestToReco", jets_.refdyn_isClosestToReco,"refdyn_isClosestToReco[nref]/O");
    t->Branch("jtdyn_refdyn_dR", jets_.jtdyn_refdyn_dR,"jtdyn_refdyn_dR[nref]/F");
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
}

void HiInclusiveJetSubstructure::analyze(const Event& iEvent, const EventSetup& iSetup) {
  int event = iEvent.id().event();
  int run = iEvent.id().run();
  int lumi = iEvent.id().luminosityBlock();

  jets_.run = run;
  jets_.evt = event;
  jets_.lumi = lumi;

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

  // auto const& pfCandidates = iEvent.get(pfCandidateToken_);
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
    if(std::abs(jet.eta()) > jetAbsEtaMax_-rParam) continue;
    //assume highest jet in event is also the trigger object, check if it's within the eta acceptance above
    if(j==0){ 
      jets_.triggerJetInAcceptance = true;
    }
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

    IterativeDeclusteringRec(groom_type, groom_combine, jet, sub1Hyb, sub2Hyb);

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
  TRandom *r1 = new TRandom3(0);
  Int_t intjet_multi = 0;
  float jet_girth = 0;
	Int_t nsplit = 0;
  double dyn_kt = std::numeric_limits<double>::min();
  Int_t dyn_split = 0;
  double z = 0;
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
        // std::cout << "Number of pfCand " << pfCandidates.size() << std::endl;
        // Geometrical PF Candidate x Jet Constituent Matching - Added by Bharadwaj - Apr 2023
        // poor man's matching, someone fix please
    // std::vector<int> vec_jet_consituent_charge;
    for(auto it = daughters.begin(); it!=daughters.end(); ++it){
      //if we want only charged constituents and the daughter charge is 0, skip it
      if(doChargedConstOnly_ && (**it).charge()==0) continue;
      double PFE_scale = 1.;
      //if it is MC, rescale the 4-momentum of the charged particles (we accept only them above) by pfCCES(+-1%)
      if(isMC_) PFE_scale = pfChargedCandidateEnergyScale_;
      //vary tracking efficiency - drop ~4% of particles within the jet
      if(isMC_ && doTrackVariation_){
        // std::cout << "doing the track variation" << std::endl;
        if(r1->Uniform(0,1)<0.05) continue;
      }
      // std::cout << "Rescaling charged pfCand energy by " << PFE_scale << std::endl;
      particles.push_back(fastjet::PseudoJet((**it).px()*PFE_scale, (**it).py()*PFE_scale, (**it).pz()*PFE_scale, (**it).energy()*PFE_scale));
      mypart.reset((**it).px()*PFE_scale, (**it).py()*PFE_scale, (**it).pz()*PFE_scale, (**it).energy()*PFE_scale);

      intjet_multi++;
      jet_girth += mypart.perp()*mypart.delta_R(myjet)/myjet.perp();

    }
    // std::cout << "Particle container has " << particles.size() << " reco particles" << std::endl;
    if(particles.empty()){
      // jets_.jtdyn_var[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
      throw(123);
    }
    // std::cout << "Clustering " << particles.size() << " number of reco particles" << std::endl;
    fastjet::ClusterSequence csiter(particles, jet_def);
    std::vector<fastjet::PseudoJet> output_jets = csiter.inclusive_jets(0);
    output_jets = sorted_by_pt(output_jets);
    fastjet::PseudoJet jj = output_jets[0];
    fastjet::PseudoJet j1;
    fastjet::PseudoJet j2;
    fastjet::PseudoJet highest_splitting;
    if(!jj.has_parents(j1,j2)){
      jets_.jtdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jtdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jtdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.jt_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.jt_girth[jets_.nref] = - std::numeric_limits<double>::max();
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
      if(k_t > dyn_kt){
        // dyn_var = dyn;
        highest_splitting = j2;
        dyn_kt = k_t;
        dyn_split = nsplit;
        dyn_deltaR = delta_R;
        dyn_z = z;
      }
      jj = j1;
      nsplit = nsplit+1;
    }

    // std::cout << highest_splitting.eta() << " " << highest_splitting.phi() << " highest reco splitting eta phi at " << dyn_split << std::endl;
    // jets_.jtdyn_var[jets_.nref] = dyn_var;
    jets_.jtdyn_split[jets_.nref] = dyn_split;
    jets_.jtdyn_deltaR[jets_.nref] = dyn_deltaR;
    jets_.jtdyn_kt[jets_.nref] = dyn_kt;
    jets_.jtdyn_z[jets_.nref] = dyn_z;
    jets_.jt_intjet_multi[jets_.nref] = intjet_multi;
    jets_.jt_girth[jets_.nref] = jet_girth;
  } 
  catch (fastjet::Error) { /*return -1;*/ }
  catch (Int_t MyNum){
    if(MyNum == 123)
      std::cout << "Whoops, seems the number of charged jet constituents is 0! Setting all reco jet split variables to numeric min." << std::endl;
    if(MyNum == 124)
      std::cout << "Jet does not have any parents, out of the loop!" << std::endl;
  }
}

void HiInclusiveJetSubstructure::IterativeDeclusteringGen(double groom_type, double groom_combine,const reco::GenJet& jet,fastjet::PseudoJet *sub1,fastjet::PseudoJet *sub2)
{
  Int_t intjet_multi = 0;
  float jet_girth = 0;
  double nsplit = 0;
  // double dyn_theta = 0;
  double dyn_kt = std::numeric_limits<double>::min();
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
      intjet_multi++;
      jet_girth += mypart.perp()*mypart.delta_R(myjet)/myjet.perp();
    }
    // std::cout << "Particle container has " << particles.size() << " reco particles" << std::endl;
    if(particles.empty()){
      // jets_.jtdyn_var[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_split[jets_.nref] = std::numeric_limits<int>::min();
      jets_.refdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.ref_girth[jets_.nref] = - std::numeric_limits<double>::max();
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
      jets_.refdyn_deltaR[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_kt[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.refdyn_z[jets_.nref] = - std::numeric_limits<double>::max();
      jets_.ref_intjet_multi[jets_.nref] = std::numeric_limits<int>::min();
      jets_.ref_girth[jets_.nref] = - std::numeric_limits<double>::max();
      throw(124);
    }
    while(jj.has_parents(j1,j2)){
      if(j1.perp() < j2.perp()) std::swap(j1,j2);
      double delta_R = j1.delta_R(j2);
      if(doHardestSplitMatching_ && isMC_) jets_.refJetConstituent.push_back(j2);
      double k_t = j2.perp()*delta_R;
      z = j2.perp()/(j1.perp()+j2.perp());
      // std::cout << "Truth split " << nsplit << " with k_T=" << k_t << " z=" << z <<  " eta " << j2.eta() << " phi " << j2.phi() <<  std::endl;
      if(k_t > dyn_kt){
        highest_splitting = j2;
        dyn_kt = k_t;
        dyn_split = nsplit;
        dyn_deltaR = delta_R;
        dyn_z = z;
      }
      jj = j1;
      nsplit = nsplit+1;
    }
    // std::cout << highest_splitting.eta() << " " << highest_splitting.phi() << " highest truth splitting eta phi at " << dyn_split << std::endl;
    jets_.refdyn_split[jets_.nref] = dyn_split;
    jets_.refdyn_deltaR[jets_.nref] = dyn_deltaR;
    jets_.refdyn_kt[jets_.nref] = dyn_kt;
    jets_.refdyn_z[jets_.nref] = dyn_z;
    jets_.ref_intjet_multi[jets_.nref] = intjet_multi;
    jets_.ref_girth[jets_.nref] = jet_girth;
  } 
  catch (fastjet::Error) { /*return -1;*/ }
  catch (Int_t MyNum){
    if(MyNum == 123)
      std::cout << "Whoops, seems the number of charged jet constituents is 0! Setting all gen jet split variables to numeric min." << std::endl;
    if(MyNum == 124)
      std::cout << "Jet does not have any parents, out of the loop!" << std::endl;
  }
}

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
