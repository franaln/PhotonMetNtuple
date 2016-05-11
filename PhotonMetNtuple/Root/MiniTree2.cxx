#include "PhotonMetNtuple/MiniTree2.h"
#include "PhotonMetNtuple/Utils.h"

#include "xAODTruth/xAODTruthHelpers.h"
#include "FourMomUtils/xAODP4Helpers.h"


#define IGEV 0.001

MiniTree2::MiniTree2(const std::string& name): asg::AsgMetadataTool( name )
{
  declareProperty("SystematicList", m_sysList); //pass here the list of systematics
  declareProperty("OutFile", m_outfile); //here we should pass *file = wk()->getOutputFile ("output");
  declareProperty("IsMC", m_ismc);

  tree = 0;
}

MiniTree2::~MiniTree2()
{
  
  for(const auto& sys : m_sysList) {
    
    if(sys.affectsKinematics){
      const CP::SystematicSet& sysSet = sys.systset;
      std::string sys_name = sysSet.name().c_str();
      
      delete ph_pt_map[sys_name];
      delete ph_eta_map[sys_name];
      delete ph_phi_map[sys_name];
      delete ph_iso_map[sys_name];
      delete ph_w_map[sys_name];

      delete ph_noniso_pt_map[sys_name];
      delete ph_noniso_eta_map[sys_name];
      delete ph_noniso_phi_map[sys_name];
      delete ph_noniso_iso_map[sys_name];
      delete ph_noniso_w_map[sys_name];
      
      delete jet_pt_map[sys_name];
      delete jet_eta_map[sys_name];
      delete jet_phi_map[sys_name];
      delete jet_e_map[sys_name];
      delete jet_isb_map[sys_name];
      delete jet_w_map[sys_name];
      
      delete el_pt_map[sys_name];
      delete el_eta_map[sys_name];
      delete el_phi_map[sys_name];
      delete el_ch_map[sys_name];
      delete el_w_map[sys_name];

      delete mu_pt_map[sys_name];
      delete mu_eta_map[sys_name];
      delete mu_phi_map[sys_name];
      delete mu_ch_map[sys_name];
      delete mu_w_map[sys_name];
    }
  }
     
  delete ph_pt_map["Nominal"];
  delete ph_eta_map["Nominal"];
  delete ph_phi_map["Nominal"];
  delete ph_iso_map["Nominal"];
  delete ph_w_map["Nominal"];

  delete ph_noniso_pt_map["Nominal"];
  delete ph_noniso_eta_map["Nominal"];
  delete ph_noniso_phi_map["Nominal"];
  delete ph_noniso_iso_map["Nominal"];
  delete ph_noniso_w_map["Nominal"];

  delete jet_pt_map["Nominal"];
  delete jet_eta_map["Nominal"];
  delete jet_phi_map["Nominal"];
  delete jet_e_map["Nominal"];
  delete jet_isb_map["Nominal"];
  delete jet_w_map["Nominal"];

  delete el_pt_map["Nominal"];
  delete el_eta_map["Nominal"];
  delete el_phi_map["Nominal"];
  delete el_ch_map["Nominal"];
  delete el_w_map["Nominal"];

  delete mu_pt_map["Nominal"];
  delete mu_eta_map["Nominal"];
  delete mu_phi_map["Nominal"];
  delete mu_ch_map["Nominal"];
  delete mu_w_map["Nominal"];
}

TString MiniTree2::BookName(TString branch, TString sys_name) 
{
  TString ret = branch + "_" + sys_name ;
  return ret;
}

StatusCode MiniTree2::initialize() 
{
  //Init the nominal tree
  ph_loose_pt = new std::vector<float>(); 
  ph_loose_eta = new std::vector<float>();
  ph_loose_phi = new std::vector<float>();
  ph_loose_isem = new std::vector<unsigned int>();
  ph_loose_iso20 = new std::vector<float>();
  ph_loose_iso40 = new std::vector<float>();

  ph_pt = new std::vector<float>(); 
  ph_eta = new std::vector<float>();
  ph_phi = new std::vector<float>();
  ph_iso = new std::vector<float>();
  ph_w = new std::vector<float>();

  ph_truth_pt = new std::vector<float>(); 
  ph_truth_eta = new std::vector<float>();
  ph_truth_phi = new std::vector<float>();
  ph_truth_id = new std::vector<int>();
  ph_truth_type = new std::vector<int>();
  ph_truth_origin = new std::vector<int>();

  ph_noniso_pt = new std::vector<float>(); 
  ph_noniso_eta = new std::vector<float>();
  ph_noniso_phi = new std::vector<float>();
  ph_noniso_iso = new std::vector<float>();
  ph_noniso_w = new std::vector<float>();

  jet_pt = new std::vector<float>(); 
  jet_eta = new std::vector<float>();
  jet_phi = new std::vector<float>();
  jet_e = new std::vector<float>();
  jet_isb = new std::vector<bool>();
  jet_w = new std::vector<float>();
  
  el_pt = new std::vector<float>(); 
  el_eta = new std::vector<float>();
  el_phi = new std::vector<float>();
  el_ch = new std::vector<int>();
  el_w = new std::vector<float>();
  
  mu_pt = new std::vector<float>(); 
  mu_eta = new std::vector<float>();
  mu_phi = new std::vector<float>();
  mu_ch = new std::vector<int>();
  mu_w = new std::vector<float>();

  // Tree  
  TString tree_name = "mini";

  tree = new TTree(tree_name, tree_name);
  tree->SetDirectory(m_outfile);	 

  // Nominal blocks
  tree->Branch("event_number", &event_number, "event_number/I");
  tree->Branch("avg_mu", &avg_mu, "avg_mu/I");

  tree->Branch("weight_mc", &weight_mc);
  tree->Branch("weight_pu", &weight_pu);
  tree->Branch("weight_sf", &weight_sf_map["Nominal"]);
  tree->Branch("weight_btag", &weight_btag_map["Nominal"]);
 
  // loose photon only in nominal
  tree->Branch("ph_loose_n", &ph_loose_n, "ph_loose_n/I");
  tree->Branch("ph_loose_eta", ph_loose_eta);
  tree->Branch("ph_loose_phi", ph_loose_phi);
  tree->Branch("ph_loose_pt",  ph_loose_pt);
  tree->Branch("ph_loose_iso20",  ph_loose_iso20);
  tree->Branch("ph_loose_iso40",  ph_loose_iso40);

  tree->Branch("ph_n", &ph_n_map["Nominal"], "ph_n/I");
  tree->Branch("ph_pt",  ph_pt);
  tree->Branch("ph_eta", ph_eta);
  tree->Branch("ph_phi", ph_phi);
  tree->Branch("ph_iso", ph_iso);
  tree->Branch("ph_w",   ph_w);

  if (m_ismc) {
    tree->Branch("ph_truth_pt",  ph_truth_pt);
    tree->Branch("ph_truth_eta", ph_truth_eta);
    tree->Branch("ph_truth_phi", ph_truth_phi);
    tree->Branch("ph_truth_id", ph_truth_id);
    tree->Branch("ph_truth_type", ph_truth_type);
    tree->Branch("ph_truth_origin", ph_truth_origin);
  }

  tree->Branch("ph_noniso_n", &ph_noniso_n_map["Nominal"], "ph_noniso_n/I");
  tree->Branch("ph_noniso_pt",  ph_noniso_pt);
  tree->Branch("ph_noniso_eta", ph_noniso_eta);
  tree->Branch("ph_noniso_phi", ph_noniso_phi);
  tree->Branch("ph_noniso_iso", ph_noniso_iso);
  tree->Branch("ph_noniso_w",   ph_noniso_w);

  tree->Branch("jet_n", &jet_n_map["Nominal"], "jet_n/I");
  tree->Branch("bjet_n", &bjet_n_map["Nominal"], "bjet_n/I");
  tree->Branch("jet_eta", jet_eta);
  tree->Branch("jet_phi", jet_phi);
  tree->Branch("jet_pt",  jet_pt);
  tree->Branch("jet_e",   jet_e);
  tree->Branch("jet_isb", jet_isb);
  tree->Branch("jet_w",  jet_w);
  
  tree->Branch("el_n", &el_n_map["Nominal"], "el_n/I");
  tree->Branch("el_eta", el_eta);
  tree->Branch("el_phi", el_phi);
  tree->Branch("el_pt",  el_pt);
  tree->Branch("el_ch",  el_ch);
  tree->Branch("el_w",   el_w);
  
  tree->Branch("mu_n", &mu_n_map["Nominal"], "mu_n/I");
  tree->Branch("mu_eta", mu_eta);
  tree->Branch("mu_phi", mu_phi);
  tree->Branch("mu_pt",  mu_pt);
  tree->Branch("mu_ch",  mu_ch);
  tree->Branch("mu_w",   mu_w);
  
  tree->Branch("met_et", &met_et_map["Nominal"]);
  tree->Branch("met_phi", &met_phi_map["Nominal"]);

  tree->Branch("ht", &ht_map["Nominal"]);				
  tree->Branch("meff", &meff_map["Nominal"]);				
  tree->Branch("rt2", &rt2_map["Nominal"]);				
  tree->Branch("rt4", &rt4_map["Nominal"]);				

  tree->Branch("dphi_gamjet", &dphi_gamjet_map["Nominal"]);
  tree->Branch("dphi_jetmet", &dphi_jetmet_map["Nominal"]);
  tree->Branch("dphi_gammet", &dphi_gammet_map["Nominal"]);

  std::string sys_name = "Nominal"; 
   
  weight_sf_map.insert(std::pair<std::string, float>(sys_name, 1.));
  weight_btag_map.insert(std::pair<std::string, float>(sys_name, 1.));
  
  ph_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  ph_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, ph_pt));
  ph_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_eta));
  ph_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_phi));
  ph_iso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_iso));
  ph_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_w));

  ph_noniso_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  ph_noniso_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_pt));
  ph_noniso_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_eta));
  ph_noniso_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_phi));
  ph_noniso_iso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_iso));
  ph_noniso_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_w));

  jet_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  bjet_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  jet_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, jet_pt));
  jet_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_eta));
  jet_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_phi));
  jet_e_map.insert  (std::pair<std::string, std::vector<float>*>(sys_name, jet_e));
  jet_isb_map.insert(std::pair<std::string, std::vector<bool>*> (sys_name, jet_isb));
  jet_w_map.insert  (std::pair<std::string, std::vector<float>*>(sys_name, jet_w));

  el_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  el_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, el_pt));
  el_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_eta));
  el_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_phi));
  el_ch_map.insert (std::pair<std::string, std::vector<int>*>  (sys_name, el_ch));
  el_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_w));

  mu_n_map.insert (std::pair<std::string, int>(sys_name, 0));  
  mu_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, mu_pt));
  mu_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_eta));
  mu_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_phi));
  mu_ch_map.insert (std::pair<std::string, std::vector<int>*>  (sys_name, mu_ch));
  mu_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_w));
  
  met_et_map.insert(std::pair<std::string, float>(sys_name, 0.));
  met_phi_map.insert(std::pair<std::string, float>(sys_name, 0.));

  ht_map.insert(std::pair<std::string, float>(sys_name, 0.));
  meff_map.insert(std::pair<std::string, float>(sys_name, 0.));
  rt2_map.insert(std::pair<std::string, float>(sys_name, 0.));
  rt4_map.insert(std::pair<std::string, float>(sys_name, 0.));

  dphi_jetmet_map.insert(std::pair<std::string, float>(sys_name, 0.));
  dphi_gamjet_map.insert(std::pair<std::string, float>(sys_name, 0.));
  dphi_gammet_map.insert(std::pair<std::string, float>(sys_name, 0.));

  // Systematics blocks
  for (const auto& sys : m_sysList) {

    if (sys.affectsKinematics || sys.affectsWeights ) {

      bool syst_affectsPhotons = ST::testAffectsObject(xAOD::Type::Photon, sys.affectsType);
      bool syst_affectsElectrons = ST::testAffectsObject(xAOD::Type::Electron, sys.affectsType);
      bool syst_affectsMuons = ST::testAffectsObject(xAOD::Type::Muon, sys.affectsType);
      bool syst_affectsJets = ST::testAffectsObject(xAOD::Type::Jet, sys.affectsType);
      bool syst_affectsTaus = ST::testAffectsObject(xAOD::Type::Tau, sys.affectsType);
      bool syst_affectsBTag = ST::testAffectsObject(xAOD::Type::BTag, sys.affectsType);

      //if (syst_affectsPhotons && !syst_affectsElectrons) continue;
      if (syst_affectsTaus) continue;
      
      ph_pt = new std::vector<float>(); 
      ph_eta = new std::vector<float>();
      ph_phi = new std::vector<float>();
      ph_iso = new std::vector<float>();
      ph_w = new std::vector<float>();

      ph_noniso_pt = new std::vector<float>(); 
      ph_noniso_eta = new std::vector<float>();
      ph_noniso_phi = new std::vector<float>();
      ph_noniso_iso = new std::vector<float>();
      ph_noniso_w = new std::vector<float>();

      jet_pt = new std::vector<float>(); 
      jet_eta = new std::vector<float>();
      jet_phi = new std::vector<float>();
      jet_e = new std::vector<float>();
      jet_isb = new std::vector<bool>();
      jet_w = new std::vector<float>();

      el_pt = new std::vector<float>(); 
      el_eta = new std::vector<float>();
      el_phi = new std::vector<float>();
      el_ch = new std::vector<int>();
      el_w = new std::vector<float>();

      mu_pt = new std::vector<float>(); 
      mu_eta = new std::vector<float>();
      mu_phi = new std::vector<float>();
      mu_ch = new std::vector<int>();
      mu_w = new std::vector<float>();

      const CP::SystematicSet& sysSet = sys.systset;
      sys_name = sysSet.name().c_str();

      TString tree_name = "mini_" + sys_name;

      // jets/btagging
      if (syst_affectsJets || syst_affectsBTag) {
        
        if (sys.affectsWeights) {
          tree->Branch(BookName("jet_w", sys_name), jet_w);
        }
        else {
          tree->Branch(BookName("jet_n", sys_name) , &jet_n_map[sys_name]);
          tree->Branch(BookName("jet_pt", sys_name) , jet_pt);
          tree->Branch(BookName("jet_eta", sys_name), jet_eta);
          tree->Branch(BookName("jet_phi", sys_name), jet_phi);
          tree->Branch(BookName("jet_e", sys_name)  , jet_e);
          tree->Branch(BookName("jet_isb", sys_name), jet_isb);
        }
      }
     
      // photons
      if (syst_affectsPhotons) {
        if (sys.affectsWeights) {
          tree->Branch(BookName("ph_w", sys_name), ph_w);
        }
        else {
          tree->Branch(BookName("ph_n", sys_name) , &ph_n_map[sys_name]);
          tree->Branch(BookName("ph_pt", sys_name) , ph_pt);
          tree->Branch(BookName("ph_eta", sys_name), ph_eta);
          tree->Branch(BookName("ph_phi", sys_name), ph_phi);
          tree->Branch(BookName("ph_iso", sys_name), ph_phi);
	  	}
      }
      
      // electrons
      if (syst_affectsElectrons) {
        if (sys.affectsWeights) {
          tree->Branch(BookName("el_w", sys_name), el_w);
        }
        else {
          tree->Branch(BookName("el_n", sys_name), &el_n_map[sys_name]);
          tree->Branch(BookName("el_pt", sys_name), el_pt);
          tree->Branch(BookName("el_eta", sys_name), el_eta);
          tree->Branch(BookName("el_phi", sys_name), el_phi);
          tree->Branch(BookName("el_ch", sys_name), el_ch);
	  	}
      }

      // muons
      if (syst_affectsMuons) {
        if (sys.affectsWeights) {
          tree->Branch(BookName("mu_w", sys_name), mu_w);
        }
        else {
          tree->Branch(BookName("mu_n", sys_name) , &mu_n_map[sys_name]);
          tree->Branch(BookName("mu_pt", sys_name) , mu_pt);
          tree->Branch(BookName("mu_eta", sys_name), mu_eta);
          tree->Branch(BookName("mu_phi", sys_name), mu_phi);
          tree->Branch(BookName("mu_ch", sys_name) , mu_ch);
	  	}
      }
      
      if (sys.affectsKinematics) {
        tree->Branch(BookName("met_et", sys_name), &met_et_map[sys_name]);
        tree->Branch(BookName("met_phi", sys_name), &met_phi_map[sys_name]);
      
        tree->Branch(BookName("ht", sys_name), &ht_map[sys_name]);
        tree->Branch(BookName("meff", sys_name), &meff_map[sys_name]);
        tree->Branch(BookName("rt2", sys_name), &rt2_map[sys_name]);
        tree->Branch(BookName("rt4", sys_name), &rt4_map[sys_name]);
      
        tree->Branch(BookName("dphi_gamjet", sys_name), &dphi_gamjet_map[sys_name]);
        tree->Branch(BookName("dphi_jetmet", sys_name), &dphi_jetmet_map[sys_name]);
        tree->Branch(BookName("dphi_gammet", sys_name), &dphi_gammet_map[sys_name]);
                    
      }
      
      // tree_map.insert(std::pair<std::string, TTree*>(sys_name, tree));
      ph_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, ph_pt));
      ph_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_eta));
      ph_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_phi));
      ph_iso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_iso));
      ph_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_w));

      // ph_noniso_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_pt));
      // ph_noniso_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_eta));
      // ph_noniso_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_phi));
      // ph_noniso_iso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_iso));
      // ph_noniso_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_w));

      jet_pt_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_pt));
      jet_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_eta));
      jet_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_phi));
      jet_e_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_e));
      jet_isb_map.insert(std::pair<std::string, std::vector<bool>*>(sys_name, jet_isb));
      jet_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_w));

      el_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, el_pt));
      el_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_eta));
      el_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_phi));
      el_ch_map.insert (std::pair<std::string, std::vector<int>*>  (sys_name, el_ch));
      el_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_w));

      mu_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, mu_pt));
      mu_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_eta));
      mu_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_phi));
      mu_ch_map.insert (std::pair<std::string, std::vector<int>*>  (sys_name, mu_ch));
      mu_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_w));

      met_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));

      ht_map.insert(std::pair<std::string, float>(sys_name, -99.));
      meff_map.insert(std::pair<std::string, float>(sys_name, -99.));
      rt2_map.insert(std::pair<std::string, float>(sys_name, -99.));
      rt4_map.insert(std::pair<std::string, float>(sys_name, -99.));

      dphi_gamjet_map.insert(std::pair<std::string, float>(sys_name, -99.));
      dphi_jetmet_map.insert(std::pair<std::string, float>(sys_name, -99.));
      dphi_gammet_map.insert(std::pair<std::string, float>(sys_name, -99.));

    }
  }

  return StatusCode::SUCCESS;
}

bool MiniTree2::process(AnalysisCollections collections, std::string sysname) 
{
  // Setting up a barebone output tree
  ph_pt_map[sysname]->clear();
  ph_eta_map[sysname]->clear();
  ph_phi_map[sysname]->clear();
  ph_iso_map[sysname]->clear();
  ph_w_map[sysname]->clear();

  if (m_ismc) {
    ph_truth_pt->clear();
    ph_truth_eta->clear();
    ph_truth_phi->clear();
    ph_truth_id->clear();
    ph_truth_type->clear();
    ph_truth_origin->clear();
  }

  if (sysname == "Nominal") {
    ph_noniso_pt_map[sysname]->clear();
    ph_noniso_eta_map[sysname]->clear();
    ph_noniso_phi_map[sysname]->clear();
    ph_noniso_iso_map[sysname]->clear();
    ph_noniso_w_map[sysname]->clear();
  }

  jet_pt_map[sysname]->clear();
  jet_eta_map[sysname]->clear();
  jet_phi_map[sysname]->clear();
  jet_e_map[sysname]->clear();
  jet_isb_map[sysname]->clear();
  jet_w_map[sysname]->clear();
  
  el_pt_map[sysname]->clear();
  el_eta_map[sysname]->clear();
  el_phi_map[sysname]->clear();
  el_ch_map[sysname]->clear();
  el_w_map[sysname]->clear();

  mu_pt_map[sysname]->clear();
  mu_eta_map[sysname]->clear();
  mu_phi_map[sysname]->clear();
  mu_ch_map[sysname]->clear();
  mu_w_map[sysname]->clear();
  
  // Setting up some basic filtering rule
  collections.photons->setStore(collections.photons_aux);
  collections.electrons->setStore(collections.electrons_aux);
  collections.muons->setStore(collections.muons_aux);
  collections.jets->setStore(collections.jets_aux);
  collections.met->setStore(collections.met_aux);

  // all loose photons
  if (sysname == "Nominal") {

    ph_loose_pt->clear();
    ph_loose_eta->clear();
    ph_loose_phi->clear();
    ph_loose_iso20->clear();
    ph_loose_iso40->clear();

    int loose_photons = 0;
    for (const auto& ph_itr : *collections.photons) {
      
      if (ph_itr->auxdata<char>("baseline") == 1  &&
          ph_itr->auxdata<char>("passOR") == 1  &&
          PassEtaCut(ph_itr, 2.37)) {
        
        loose_photons += 1;
        
        ph_loose_pt->push_back(ph_itr->pt()*IGEV);
        ph_loose_eta->push_back(ph_itr->eta());
        ph_loose_phi->push_back(ph_itr->phi());
        
        ph_loose_iso20->push_back(ph_itr->isolationValue(xAOD::Iso::topoetcone20)*IGEV);
        ph_loose_iso40->push_back(ph_itr->isolationValue(xAOD::Iso::topoetcone40)*IGEV);

        // bool is_signal = ph_itr->auxdata<char>("signal") && (ph_itr->pt() > 125000.) && (ph_itr->eta() < 2.37);
        //       // ph_loose_signal->push_back(is_signal);

        //       // unsigned int m_isem = ph_itr->selectionisEM(xAOD::EgammaParameters::SelectionisEM::isEMTight);
        //       // std::cout << m_isem << std::endl;
      }
    }
    ph_loose_n = loose_photons;
  }

  Double_t total_weight_sf = 1.;

  // electrons
  int el_n = 0;;
  for (const auto& el_itr : *collections.electrons) {
    
    if (el_itr->auxdata<char>("baseline") == 1 &&
        el_itr->auxdata<char>("passOR") == 1 &&
        el_itr->auxdata<char>("signal") == 1 && 
        PassEtaCut(el_itr)) {
      
      el_n += 1;
      el_pt_map[sysname]->push_back(el_itr->pt()*IGEV);
      el_eta_map[sysname]->push_back(el_itr->eta());
      el_phi_map[sysname]->push_back(el_itr->phi());
      el_ch_map[sysname]->push_back(el_itr->trackParticle()->charge());
      el_w_map[sysname]->push_back( el_itr->auxdata<double>("effscalefact") );
      total_weight_sf *= el_itr->auxdata<double>("effscalefact");
    }
  }
  el_n_map[sysname] = el_n;

  // muons
  int mu_n = 0;
  for (const auto& mu_itr : *collections.muons) {
    
    if (mu_itr->auxdata<char>("baseline") == 1 &&
        mu_itr->auxdata<char>("passOR") == 1 &&
        mu_itr->auxdata<char>("signal") == 1 &&
        PassEtaCut(mu_itr)) {

      mu_n += 1;      
      mu_pt_map[sysname] ->push_back(mu_itr->pt()*IGEV);
      mu_eta_map[sysname]->push_back(mu_itr->eta());
      mu_phi_map[sysname]->push_back(mu_itr->phi());
      mu_ch_map[sysname] ->push_back(mu_itr->primaryTrackParticle()->charge());
      mu_w_map[sysname]->push_back(mu_itr->auxdata<double>("effscalefact"));

      total_weight_sf *= mu_itr->auxdata<double>("effscalefact");
    }
  }
  mu_n_map[sysname] = mu_n;

  // jets
  int jet_n = 0;
  int bjet_n = 0;
  for (const auto& jet_itr : *collections.jets) {     

    if (jet_itr->auxdata<char>("baseline") == 1  &&
        jet_itr->auxdata<char>("passOR") == 1 &&
        jet_itr->auxdata<char>("signal") == 1 &&
        PassEtaCut(jet_itr, 2.5) && 
        jet_itr->pt()>30000.) {
      
      jet_n++;
      jet_pt_map[sysname]->push_back(jet_itr->pt()*IGEV);
      jet_eta_map[sysname]->push_back(jet_itr->eta());
      jet_phi_map[sysname]->push_back(jet_itr->phi());
      jet_e_map[sysname]->push_back(jet_itr->e()*IGEV);
      jet_w_map[sysname]->push_back(jet_itr->auxdata<double>("effscalefact"));

      total_weight_sf *= jet_itr->auxdata<double>("effscalefact");

      int isbjet = int(jet_itr->auxdata<char>("bjet"));
      if (isbjet)
        bjet_n++;
      
      jet_isb_map[sysname]->push_back(isbjet);
    }
  }
  jet_n_map[sysname] = jet_n;
  bjet_n_map[sysname] = bjet_n;

  // photons
  int ph_n = 0;
  int ph_noniso_n = 0;

  for (const auto& ph_itr : *collections.photons) {
   
    if (ph_itr->auxdata<char>("baseline") == 1 &&
        ph_itr->auxdata<char>("signal")   == 1 &&
        ph_itr->auxdata<char>("passOR")   == 1 &&
        PassEtaCut(ph_itr, 2.37) ) {
      
      float iso40 = ph_itr->isolationValue(xAOD::Iso::topoetcone40)*IGEV;
      float iso = iso40 - 0.022 * ph_itr->pt()*IGEV;

      if (iso < 2.45) {

        ph_n += 1;
      
        ph_pt_map[sysname] ->push_back(ph_itr->pt()*IGEV);
        ph_eta_map[sysname]->push_back(ph_itr->eta());
        ph_phi_map[sysname]->push_back(ph_itr->phi());
        
        ph_iso_map[sysname]->push_back(iso);

        double sf = ph_itr->auxdata<double>("effscalefact");
        ph_w_map[sysname]->push_back(sf);
        total_weight_sf *= sf;
        
        // truth info
        if (m_ismc) {
          const xAOD::TruthParticle *ph_truth = xAOD::TruthHelpers::getTruthParticle(*ph_itr);
          
          if (ph_truth) {
            ph_truth_pt->push_back(ph_truth->pt()*IGEV);
            ph_truth_eta->push_back(ph_truth->eta());
            ph_truth_phi->push_back(ph_truth->phi());
            ph_truth_id->push_back(ph_truth->pdgId());
          }
          else {
            ph_truth_pt->push_back(-99.);
            ph_truth_eta->push_back(-99.);
            ph_truth_phi->push_back(-99.);
            ph_truth_id->push_back(-99);
          }

          ph_truth_type->push_back(xAOD::TruthHelpers::getParticleTruthType(*ph_itr));
          ph_truth_origin->push_back(xAOD::TruthHelpers::getParticleTruthOrigin(*ph_itr));

        }

      }
      else if (sysname == "Nominal" && iso > 5.45 && iso < 29.45) {

        ph_noniso_n += 1;
      
        ph_noniso_pt_map[sysname] ->push_back(ph_itr->pt()*IGEV);
        ph_noniso_eta_map[sysname]->push_back(ph_itr->eta());
        ph_noniso_phi_map[sysname]->push_back(ph_itr->phi());
        ph_noniso_w_map[sysname]->push_back(ph_itr->auxdata<double>("effscalefact") );

        total_weight_sf *= ph_itr->auxdata<double>("effscalefact");
        
        ph_noniso_iso_map[sysname]->push_back(iso);
      }

    }
  }
  ph_n_map[sysname] = ph_n;
  ph_noniso_n_map[sysname] = ph_noniso_n;

  // met
  float etmiss_etx = 0.;
  float etmiss_ety = 0.;
  float etmiss_et = 0.;
  
  xAOD::MissingETContainer::const_iterator met_it = collections.met->find("Final");
  if (met_it == collections.met->end()) {
    Error("PhotonMetNtuple:MiniTree", "No RefFinal inside MET container");
  }
  else {
    etmiss_etx = (*met_it)->mpx() * IGEV;
    etmiss_ety = (*met_it)->mpy() * IGEV;
    etmiss_et = sqrt(etmiss_etx*etmiss_etx + etmiss_ety*etmiss_ety);
  }
  
  TLorentzVector met(etmiss_etx, etmiss_ety, 0., etmiss_et);
  met_et_map[sysname] = met.Perp();
  met_phi_map[sysname] = met.Phi();  
  
  // Extra variables
  Double_t sum_jet_pt = 0.0;
  Double_t sum_jet2_pt = 0.0;
  Double_t sum_jet4_pt = 0.0;

  for (int i=0; i<jet_n; i++) {

    Double_t pt = (*jet_pt_map[sysname])[i];

    if (jet_n >= 2 && i < 2) 
      sum_jet2_pt += pt;
    if (jet_n >= 4 && i < 4) 
      sum_jet4_pt += pt;
    
    sum_jet_pt += pt;
  }
  
  // Ht
  Double_t ht = sum_jet_pt;
  if (ph_n > 0)
    ht += (*ph_pt_map[sysname])[0];
  
  ht_map[sysname] = ht;

  // Meff
  meff_map[sysname] = ht + met_et_map[sysname];
  
  // Rt
  rt2_map[sysname] = sum_jet2_pt/sum_jet_pt;
  rt4_map[sysname] = sum_jet4_pt/sum_jet_pt;
  
  // min dphi between met and the first two jets
  Double_t dphi1 = 4.;
  Double_t dphi2 = 4.;
  if (jet_n > 0) dphi1 = get_dphi((*jet_phi_map[sysname])[0], met.Phi());
  if (jet_n > 1) dphi2 = get_dphi((*jet_phi_map[sysname])[1], met.Phi());
  
  dphi_jetmet_map[sysname] = TMath::Min(dphi1, dphi2);

  // dphi between leading photon and leading jet
  if (ph_n > 0 && jet_n > 0) 
    dphi_gamjet_map[sysname] = get_dphi((*ph_phi_map[sysname])[0], (*jet_phi_map[sysname])[0]);
  
  // dphi between leading photon and MET
  if (ph_n > 0)
    dphi_gammet_map[sysname] = get_dphi((*ph_phi_map[sysname])[0], met.Phi());

  // weigths
  //weight_mc = collections.weight_mc;
  //weight_pu = collections.weight_pu;
  weight_sf_map[sysname] = total_weight_sf;

  
  // Skim: at least one loose photon with pt>75 or an electron
  int photons_loose = 0;
  for (const auto& ph_itr : *collections.photons) {
    if (ph_itr->auxdata<char>("baseline") == 1  &&
        ph_itr->auxdata<char>("passOR") == 1  &&
        PassEtaCut(ph_itr, 2.37) &&
        ph_itr->pt()*IGEV > 75) {
      photons_loose += 1;
    }
  }
  return (photons_loose > 0 || el_n > 0);
  //return true;
}

// Call after all the syst have been processed and ONLY IF one of them passed the event selection criteria
StatusCode MiniTree2::FillTree()
{
  tree->Fill();
  return StatusCode::SUCCESS;
}

bool MiniTree2::PassEtaCut(const xAOD::IParticle *part, Double_t maxeta) 
{
  Double_t eta = fabs(part->eta());

  if (eta > maxeta)
    return false;

  if (eta >= 1.37 && eta < 1.52)
    return false;

  return true;
}
