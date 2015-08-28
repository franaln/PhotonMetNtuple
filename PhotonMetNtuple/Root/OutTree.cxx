#include "PhotonMetNtuple/OutTree.h"
#include "PhotonMetNtuple/Utils.h"

#define IGEV 0.001

OutTree::OutTree(const std::string& name): asg::AsgMetadataTool( name )
{
  
  declareProperty("SystematicList", m_sysList); //pass here the list of systematics
  declareProperty("OutFile", m_outfile); //here we should pass *file = wk()->getOutputFile ("output");
  
  tree = 0;
}

OutTree::~OutTree()
{
  
  for(const auto& sys : m_sysList) {
    
    if(sys.affectsKinematics){
      const CP::SystematicSet& sysSet = sys.systset;
      std::string sys_name = sysSet.name().c_str();
      
      delete ph_pt_map[sys_name];
      delete ph_eta_map[sys_name];
      delete ph_phi_map[sys_name];
      delete ph_w_map[sys_name];
      
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
  delete ph_w_map["Nominal"];

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

TString OutTree::BookName(TString branch, TString sys_name) 
{
  TString ret = branch + "_" + sys_name ;
  return ret;
}

StatusCode OutTree::initialize() 
{
  //Init the nominal tree
  ph_pt = new std::vector<float>(); 
  ph_eta = new std::vector<float>();
  ph_phi = new std::vector<float>();
  ph_w = new std::vector<float>();

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
  
  TString tree_name = "mini";

  //Get ready to leak memory all over the place
  tree = new TTree(tree_name, tree_name);
  tree->SetDirectory(m_outfile);	 

  tree->Branch("event_number", &event_number, "event_number/I");
  tree->Branch("weight_mc", &weight_mc);
  tree->Branch("weight_pu", &weight_pu);
           
  tree->Branch("ph_n", &ph_n_map["Nominal"], "ph_n/I");
  tree->Branch("ph_eta", ph_eta);
  tree->Branch("ph_phi", ph_phi);
  tree->Branch("ph_pt",  ph_pt);
  tree->Branch("ph_w",   ph_w);

  tree->Branch("jet_n", &jet_n_map["Nominal"], "jet_n/I");
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
  tree->Branch("rt2", &rt2_map["Nominal"]);				
  tree->Branch("rt4", &rt4_map["Nominal"]);				

  tree->Branch("dphi_gamjet", &dphi_gamjet_map["Nominal"]);
  tree->Branch("dphi_jetmet", &dphi_jetmet_map["Nominal"]);
  tree->Branch("dphi_gammet", &dphi_gammet_map["Nominal"]);


  std::string sys_name = "Nominal"; 
   
  //  SFwei_map.insert (std::pair<std::string, float>(sys_name, 1.));
  //bTagSF_map.insert (std::pair<std::string, float>(sys_name, 1.));
  
  ph_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  ph_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, ph_pt));
  ph_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_eta));
  ph_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_phi));
  ph_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_w));

  jet_n_map.insert (std::pair<std::string, int>(sys_name, 0));
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
  rt2_map.insert(std::pair<std::string, float>(sys_name, 0.));
  rt4_map.insert(std::pair<std::string, float>(sys_name, 0.));

  dphi_jetmet_map.insert(std::pair<std::string, float>(sys_name, 0.));
  dphi_gamjet_map.insert(std::pair<std::string, float>(sys_name, 0.));
  dphi_gammet_map.insert(std::pair<std::string, float>(sys_name, 0.));
   
  //Here I have to initialize all the trees
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
      ph_w = new std::vector<float>();

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

      //Get ready to leak memory all over the place
      if (syst_affectsJets || syst_affectsBTag) {
        if (sys.affectsWeights) {
          tree->Branch(BookName("jet_w", sys_name), jet_w);
        }
        else {
          tree->Branch(BookName("jet_pt", sys_name) , jet_pt);
          tree->Branch(BookName("jet_eta", sys_name), jet_eta);
          tree->Branch(BookName("jet_phi", sys_name), jet_phi);
          tree->Branch(BookName("jet_e", sys_name)  , jet_e);
          tree->Branch(BookName("jet_isb", sys_name), jet_isb);
        }
      }
      

      if (syst_affectsPhotons) {
        if (sys.affectsWeights) {
          tree->Branch(BookName("ph_w", sys_name), ph_w);
        }
        else {
          tree->Branch(BookName("ph_pt", sys_name) , ph_pt);
          tree->Branch(BookName("ph_eta", sys_name), ph_eta);
          tree->Branch(BookName("ph_phi", sys_name), ph_phi);
	  	}
      }
      
      if (syst_affectsElectrons) {
        if (sys.affectsWeights) {
          tree->Branch(BookName("el_w", sys_name), el_w);
        }
        else {
          tree->Branch(BookName("el_pt", sys_name), el_pt);
          tree->Branch(BookName("el_eta", sys_name), el_eta);
          tree->Branch(BookName("el_phi", sys_name), el_phi);
          tree->Branch(BookName("el_ch", sys_name), el_ch);
	  	}
      }

      if (syst_affectsMuons) {
        if (sys.affectsWeights) {
          tree->Branch(BookName("mu_w", sys_name), mu_w);
        }
        else {
          tree->Branch(BookName("mu_pt", sys_name) , mu_pt);
          tree->Branch(BookName("mu_eta", sys_name), mu_eta);
          tree->Branch(BookName("mu_phi", sys_name), mu_phi);
          tree->Branch(BookName("mu_ch", sys_name) , mu_ch);
	  	}
      }
      
      if (sys.affectsKinematics) {
        tree->Branch(BookName("met_et", sys_name), &met_et_map[sys_name]);
        tree->Branch(BookName("met_phi", sys_name), &met_phi_map[sys_name]);
      }
      
      // tree_map.insert(std::pair<std::string, TTree*>(sys_name, tree));
      ph_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, ph_pt));
      ph_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_eta));
      ph_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_phi));
      ph_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_w));

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

      met_et_map.insert(std::pair<std::string, float>(sys_name, 0.));
      met_phi_map.insert(std::pair<std::string, float>(sys_name, 0.));
    }
  }

  return StatusCode::SUCCESS;
}


bool OutTree::process(AnalysisCollections collections, std::string sysname) 
{
  const char *APP_NAME = "process()";
  
  //Setting up a barebone output tree
  ph_pt_map[sysname]->clear();
  ph_eta_map[sysname]->clear();
  ph_phi_map[sysname]->clear();
  ph_w_map[sysname]->clear();
  
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
  
  //Setting up some basic filtering rule
  event_number = collections.event_number;
  collections.photons->setStore(collections.photons_aux);
  collections.electrons->setStore(collections.electrons_aux);
  collections.muons->setStore(collections.muons_aux);
  collections.jets->setStore(collections.jets_aux);
  collections.met->setStore(collections.met_aux);
  
  // electrons
  int el_n = 0;;
  for (const auto& el_itr : *collections.electrons) {
    
    if (el_itr->auxdata<char>("passOR") == 1 &&
        el_itr->auxdata<char>("signal") == 1 &&
        el_itr->pt() > 25000. && 
        (fabs( el_itr->eta()) < 2.47) ){
      
      el_n += 1;
      el_pt_map[sysname]->push_back(el_itr->pt()*IGEV);
      el_eta_map[sysname]->push_back(el_itr->eta());
      el_phi_map[sysname]->push_back(el_itr->phi());
      el_ch_map[sysname]->push_back(el_itr->trackParticle()->charge());
      el_w_map[sysname]->push_back( el_itr->auxdata<double>("effscalefact") );
    }
  }
  el_n_map[sysname] = el_n;

  // muons
  int mu_n = 0;
  for (const auto& mu_itr : *collections.muons) {
    
    if (mu_itr->auxdata<char>("passOR") == 1 &&
        mu_itr->auxdata<char>("signal") == 1 &&
        mu_itr->pt() > 25000. && 
        (fabs(mu_itr->eta()) < 2.4)) {

      mu_n += 1;      
      mu_pt_map[sysname] ->push_back(mu_itr->pt()*IGEV);
      mu_eta_map[sysname]->push_back(mu_itr->eta());
      mu_phi_map[sysname]->push_back(mu_itr->phi());
      mu_ch_map[sysname] ->push_back(mu_itr->primaryTrackParticle()->charge());
      mu_w_map[sysname]->push_back(mu_itr->auxdata<double>("effscalefact"));
    }
  }
  mu_n_map[sysname] = mu_n;

  // jets
  int jet_n = 0;
  for (const auto& jet_itr : *collections.jets) {     
    if (jet_itr->auxdata<char>("baseline") == 1  &&
        jet_itr->auxdata<char>("passOR") == 1  &&
        jet_itr->pt() > 20000. && 
        fabs(jet_itr->eta()) < 2.5) {  

      jet_n += 1;      
      jet_pt_map[sysname]->push_back(jet_itr->pt()*IGEV);
      jet_eta_map[sysname]->push_back(jet_itr->eta());
      jet_phi_map[sysname]->push_back(jet_itr->phi());
      jet_e_map[sysname]->push_back(jet_itr->e()*IGEV);
      jet_w_map[sysname]->push_back(jet_itr->auxdata<double>("effscalefact"));
      jet_isb_map[sysname]->push_back(jet_itr->auxdata<char>("bjet")==1);
    }
  }
  jet_n_map[sysname] = jet_n;

  // photons
  int ph_n = 0;
  for (const auto& ph_itr : *collections.photons) {
   
    if (ph_itr->auxdata<char>("signal") == 1  &&
        ph_itr->auxdata<char>("passOR") == 1  &&
        ph_itr->pt() > 75000. && fabs(ph_itr->eta()) < 2.37) {  
       
      ph_n += 1;

      ph_pt_map[sysname] ->push_back(ph_itr->pt()*IGEV);
      ph_eta_map[sysname]->push_back(ph_itr->eta());
      ph_phi_map[sysname]->push_back(ph_itr->phi());
      
      ph_w_map[sysname]->push_back(ph_itr->auxdata<double>("effscalefact") );
    }
  }
  ph_n_map[sysname] = ph_n;
  
  
  float etmiss_etx = 0.;
  float etmiss_ety = 0.;
  float etmiss_et = 0.;
  
  xAOD::MissingETContainer::const_iterator met_it = collections.met->find("Final");
  if (met_it == collections.met->end()) {
    Error(APP_NAME, "No RefFinal inside MET container");
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
  
  /// Ht
  Double_t ht = sum_jet_pt;
  if (ph_n > 0)
    ht += (*ph_pt_map[sysname])[0];
  
  ht_map[sysname] = ht;
  
  /// Rt
  rt2_map[sysname] = sum_jet2_pt/sum_jet_pt;
  rt4_map[sysname] = sum_jet4_pt/sum_jet_pt;
  
  // Delta phi between met and closest jet
  Double_t dphi1 = 4.;
  Double_t dphi2 = 4.;
  
  if (jet_n > 0) dphi1 = get_dphi((*jet_phi_map[sysname])[0], met.Phi());
  if (jet_n > 1) dphi2 = get_dphi((*jet_phi_map[sysname])[1], met.Phi());
  
  dphi_jetmet_map[sysname] = TMath::Min(dphi1, dphi2);
  
  // dphi beteen leading photon and leading jet
  if (ph_n > 0 && jet_n > 0) 
    dphi_gamjet_map[sysname] = get_dphi((*ph_phi_map[sysname])[0], (*jet_phi_map[sysname])[0]);
  
  // dphi beteen leading photon and MET
  if (ph_n > 0)
    dphi_gammet_map[sysname] = get_dphi((*ph_phi_map[sysname])[0], met.Phi());

  // weigths
  weight_mc = collections.weight_mc;
  weight_pu = collections.weight_pu;
  
  //Apply your custom filter here
  return true; 
}

// Call after all the syst have been processed and ONLY IF one of them passed the event selection criteria
StatusCode OutTree::FillTree()
{
  tree->Fill();
  return StatusCode::SUCCESS;
}
