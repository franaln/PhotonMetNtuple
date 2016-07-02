#include "PhotonMetNtuple/BaselineTree.h"
#include "PhotonMetNtuple/Utils.h"

#include "xAODTruth/xAODTruthHelpers.h"
#include "FourMomUtils/xAODP4Helpers.h"

#define IGEV 0.001


BaselineTree::BaselineTree(const std::string& name): 
  asg::AsgMetadataTool(name)
{
  declareProperty("OutFile", m_outfile); //here we should pass *file = wk()->getOutputFile ("output");
  declareProperty("IsMC", m_ismc);

  tree = 0;
}

BaselineTree::~BaselineTree()
{
  delete ph_pt;
  delete ph_eta;
  delete ph_etas2;
  delete ph_phi;
  delete ph_iso;

  delete jet_pt;
  delete jet_eta;
  delete jet_phi;
  delete jet_e;
  delete jet_isb;

  delete el_pt;
  delete el_eta;
  delete el_etas2;
  delete el_phi;
  delete el_ch;

  delete mu_pt;
  delete mu_eta;
  delete mu_phi;
  delete mu_ch;

}

StatusCode BaselineTree::initialize() 
{
  //Init the nominal tree
  ph_pt = new std::vector<float>(); 
  ph_eta = new std::vector<float>();
  ph_etas2 = new std::vector<float>();
  ph_phi = new std::vector<float>();
  ph_iso = new std::vector<float>();
  ph_passOR = new std::vector<int>();
  ph_signal = new std::vector<int>();
  ph_isol = new std::vector<int>();

  jet_pt = new std::vector<float>(); 
  jet_eta = new std::vector<float>();
  jet_phi = new std::vector<float>();
  jet_e = new std::vector<float>();
  jet_isb = new std::vector<bool>();
  jet_passOR = new std::vector<int>();
  jet_signal = new std::vector<int>();
  
  el_pt = new std::vector<float>(); 
  el_eta = new std::vector<float>();
  el_etas2 = new std::vector<float>();
  el_phi = new std::vector<float>();
  el_ch = new std::vector<int>();
  el_passOR = new std::vector<int>();
  el_signal = new std::vector<int>();
  el_isol = new std::vector<int>();
  
  mu_pt = new std::vector<float>(); 
  mu_eta = new std::vector<float>();
  mu_phi = new std::vector<float>();
  mu_ch = new std::vector<int>();
  mu_passOR = new std::vector<int>();
  mu_signal = new std::vector<int>();
  mu_isol = new std::vector<int>();

  // Tree  
  TString tree_name = "baseline";

  tree = new TTree(tree_name, tree_name);
  tree->SetDirectory(m_outfile);	 

  // Nominal blocks
  tree->Branch("event", &event_number, "event/I");
  tree->Branch("run", &run_number, "run/I");
  
  tree->Branch("ph_n", &ph_n, "ph_n/I");
  tree->Branch("ph_pt",  ph_pt);
  tree->Branch("ph_eta", ph_eta);
  tree->Branch("ph_etas2", ph_etas2);
  tree->Branch("ph_phi", ph_phi);
  tree->Branch("ph_iso", ph_iso);
  tree->Branch("ph_passOR",   ph_passOR);
  tree->Branch("ph_signal",   ph_signal);
  tree->Branch("ph_isol",   ph_isol);

  tree->Branch("jet_n", &jet_n, "jet_n/I");
  tree->Branch("bjet_n", &bjet_n, "bjet_n/I");
  tree->Branch("jet_eta", jet_eta);
  tree->Branch("jet_phi", jet_phi);
  tree->Branch("jet_pt",  jet_pt);
  tree->Branch("jet_e",   jet_e);
  tree->Branch("jet_isb", jet_isb);
  tree->Branch("jet_passOR",  jet_passOR);
  tree->Branch("jet_signal",  jet_signal);
  
  tree->Branch("el_n", &el_n, "el_n/I");
  tree->Branch("el_eta", el_eta);
  tree->Branch("el_etas2", el_etas2);
  tree->Branch("el_phi", el_phi);
  tree->Branch("el_pt",  el_pt);
  tree->Branch("el_ch",  el_ch);
  tree->Branch("el_passOR",   el_passOR);
  tree->Branch("el_signal",   el_signal);
  tree->Branch("el_isol",   el_isol);

  tree->Branch("mu_n", &mu_n, "mu_n/I");
  tree->Branch("mu_eta", mu_eta);
  tree->Branch("mu_phi", mu_phi);
  tree->Branch("mu_pt",  mu_pt);
  tree->Branch("mu_ch",  mu_ch);
  tree->Branch("mu_passOR",   mu_passOR);
  tree->Branch("mu_signal",   mu_signal);
  tree->Branch("mu_isol",   mu_isol);
  
  tree->Branch("met_et", &met_et);
  tree->Branch("met_phi", &met_phi);

  return StatusCode::SUCCESS;
}

void BaselineTree::clear()
{
  event_number = 0;
  run_number = 0;
}

bool BaselineTree::process(AnalysisBaselineCollections collections) 
{
  // clear
  ph_pt->clear();
  ph_eta->clear();
  ph_etas2->clear();
  ph_phi->clear();
  ph_iso->clear();

  ph_passOR->clear();
  ph_signal->clear();
  ph_isol->clear();
  
  jet_pt->clear();
  jet_eta->clear();
  jet_phi->clear();
  jet_e->clear();
  jet_isb->clear();

  jet_passOR->clear();
  jet_signal->clear();
  
  el_pt->clear();
  el_eta->clear();
  el_etas2->clear();
  el_phi->clear();
  el_ch->clear();

  el_passOR->clear();
  el_signal->clear();
  el_isol->clear();

  mu_pt->clear();
  mu_eta->clear();
  mu_phi->clear();
  mu_ch->clear();

  mu_passOR->clear();
  mu_signal->clear();
  mu_isol->clear();


  // Setting up some basic filtering rule
  collections.photons->setStore(collections.photons_aux);
  collections.electrons->setStore(collections.electrons_aux);
  collections.muons->setStore(collections.muons_aux);
  collections.jets->setStore(collections.jets_aux);
  collections.met->setStore(collections.met_aux);


  // electrons
  el_n = 0;;
  for (const auto& el_itr : *collections.electrons) {

    if (el_itr->auxdata<char>("baseline") != 1)
      continue;

      el_n += 1;
      el_pt->push_back(el_itr->pt()*IGEV);
      el_eta->push_back(el_itr->eta());
      el_etas2->push_back(el_itr->caloCluster()->etaBE(2));
      el_phi->push_back(el_itr->phi());
      el_ch->push_back(el_itr->trackParticle()->charge());

      el_passOR->push_back(el_itr->auxdata<char>("passOR"));
      el_signal->push_back(el_itr->auxdata<char>("signal"));
      el_isol->push_back(el_itr->auxdata<char>("isol"));

  }

  // muons
  mu_n = 0;
  for (const auto& mu_itr : *collections.muons) {
    
    if (mu_itr->auxdata<char>("baseline") != 1)
      continue;

    mu_n += 1;      
    mu_pt->push_back(mu_itr->pt()*IGEV);
    mu_eta->push_back(mu_itr->eta());
    mu_phi->push_back(mu_itr->phi());
    mu_ch ->push_back(mu_itr->primaryTrackParticle()->charge());

    mu_passOR->push_back(mu_itr->auxdata<char>("passOR"));
    mu_signal->push_back(mu_itr->auxdata<char>("signal"));
    mu_isol->push_back(mu_itr->auxdata<char>("isol"));
  }

  // jets
  jet_n = 0;
  bjet_n = 0;
  for (const auto& jet_itr : *collections.jets) {     

    if (jet_itr->auxdata<char>("baseline") != 1)
      continue;

    jet_n++;
    jet_pt->push_back(jet_itr->pt()*IGEV);
    jet_eta->push_back(jet_itr->eta());
    jet_phi->push_back(jet_itr->phi());
    jet_e->push_back(jet_itr->e()*IGEV);
    
    int isbjet = int(jet_itr->auxdata<char>("bjet"));
    if (isbjet)
      bjet_n++;
      
    jet_isb->push_back(isbjet);

    jet_passOR->push_back(jet_itr->auxdata<char>("passOR"));
    jet_signal->push_back(jet_itr->auxdata<char>("signal"));
  }

  // photons
  ph_n = 0;
  for (const auto& ph_itr : *collections.photons) {

    if (ph_itr->auxdata<char>("baseline") != 1)
      continue;
      
    // separate iso and noniso photons
    float iso40 = ph_itr->isolationValue(xAOD::Iso::topoetcone40)*IGEV;
    float iso = iso40 - 0.022 * ph_itr->pt()*IGEV;

    ph_n += 1;
    ph_pt ->push_back(ph_itr->pt()*IGEV);
    ph_eta->push_back(ph_itr->eta());
    ph_etas2->push_back(ph_itr->caloCluster()->etaBE(2));
    ph_phi->push_back(ph_itr->phi());
    ph_iso->push_back(iso);
    
    ph_passOR->push_back(ph_itr->auxdata<char>("passOR"));
    ph_signal->push_back(ph_itr->auxdata<char>("signal"));
    ph_isol->push_back(ph_itr->auxdata<char>("isol"));
    
  }

  // met
  xAOD::MissingETContainer::const_iterator met_it = collections.met->find("Final");
  if (met_it == collections.met->end()) {
    Error("PhotonMetNtuple:BaselineTree", "No RefFinal inside MET container");
  }

  met_et = (*met_it)->met() * IGEV; 
  met_phi = (*met_it)->phi(); 

  return true;
}

// Call after all the syst have been processed and ONLY IF one of them passed the event selection criteria
StatusCode BaselineTree::FillTree()
{
  tree->Fill();
  return StatusCode::SUCCESS;
}

