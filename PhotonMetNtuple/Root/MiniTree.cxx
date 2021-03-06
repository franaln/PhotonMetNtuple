#include "PhotonMetNtuple/MiniTree.h"
#include "PhotonMetNtuple/Utils.h"

#include "xAODTruth/xAODTruthHelpers.h"
#include "FourMomUtils/xAODP4Helpers.h"

#define IGEV 0.001


MiniTree::MiniTree(const std::string& name): 
  asg::AsgMetadataTool(name)
{

  declareProperty("SystematicList", m_sys_list); //pass here the list of systematics
  declareProperty("OutFile", m_outfile); //here we should pass *file = wk()->getOutputFile ("output");
  declareProperty("IsMC", m_ismc=false);
  declareProperty("SaveMediumElectrons", m_save_medium_electrons=true);

  tree = 0;
}

MiniTree::~MiniTree()
{
  
  for (const auto& sys : m_sys_list) {
    
    if (sys.affectsKinematics) {
      const CP::SystematicSet& sysSet = sys.systset;
      std::string sys_name = sysSet.name().c_str();
      
      delete ph_pt_map[sys_name];
      delete ph_eta_map[sys_name];
      delete ph_etas2_map[sys_name];
      delete ph_phi_map[sys_name];
      delete ph_etcone40_map[sys_name];
      delete ph_ptcone20_map[sys_name];
      delete ph_iso_map[sys_name];
      delete ph_trackiso_map[sys_name];
      delete ph_conv_map[sys_name];
      delete ph_w_map[sys_name];

      delete ph_noniso_pt_map[sys_name];
      delete ph_noniso_eta_map[sys_name];
      delete ph_noniso_etas2_map[sys_name];
      delete ph_noniso_phi_map[sys_name];
      delete ph_noniso_etcone40_map[sys_name];
      delete ph_noniso_ptcone20_map[sys_name];
      delete ph_noniso_iso_map[sys_name];
      delete ph_noniso_trackiso_map[sys_name];
      delete ph_noniso_conv_map[sys_name];
      delete ph_noniso_w_map[sys_name];
      
      delete jet_pt_map[sys_name];
      delete jet_eta_map[sys_name];
      delete jet_phi_map[sys_name];
      delete jet_e_map[sys_name];
      delete jet_isb_map[sys_name];
      delete jet_w_map[sys_name];
      
      delete el_pt_map[sys_name];
      delete el_eta_map[sys_name];
      delete el_etas2_map[sys_name];
      delete el_phi_map[sys_name];
      delete el_ch_map[sys_name];
      delete el_w_map[sys_name];

      delete el_medium_pt_map[sys_name];
      delete el_medium_eta_map[sys_name];
      delete el_medium_etas2_map[sys_name];
      delete el_medium_phi_map[sys_name];
      delete el_medium_ch_map[sys_name];

      delete mu_pt_map[sys_name];
      delete mu_eta_map[sys_name];
      delete mu_phi_map[sys_name];
      delete mu_ch_map[sys_name];
      delete mu_w_map[sys_name];
    }
  }
     
  delete ph_pt_map["Nominal"];
  delete ph_eta_map["Nominal"];
  delete ph_etas2_map["Nominal"];
  delete ph_phi_map["Nominal"];
  delete ph_etcone40_map["Nominal"];
  delete ph_ptcone20_map["Nominal"];
  delete ph_iso_map["Nominal"];
  delete ph_trackiso_map["Nominal"];
  delete ph_conv_map["Nominal"];
  delete ph_w_map["Nominal"];

  delete ph_noniso_pt_map["Nominal"];
  delete ph_noniso_eta_map["Nominal"];
  delete ph_noniso_etas2_map["Nominal"];
  delete ph_noniso_phi_map["Nominal"];
  delete ph_noniso_etcone40_map["Nominal"];  
  delete ph_noniso_ptcone20_map["Nominal"];
  delete ph_noniso_iso_map["Nominal"];
  delete ph_noniso_trackiso_map["Nominal"];
  delete ph_noniso_conv_map["Nominal"];
  delete ph_noniso_w_map["Nominal"];

  delete jet_pt_map["Nominal"];
  delete jet_eta_map["Nominal"];
  delete jet_phi_map["Nominal"];
  delete jet_e_map["Nominal"];
  delete jet_isb_map["Nominal"];
  delete jet_w_map["Nominal"];

  delete el_pt_map["Nominal"];
  delete el_eta_map["Nominal"];
  delete el_etas2_map["Nominal"];
  delete el_phi_map["Nominal"];
  delete el_ch_map["Nominal"];
  delete el_w_map["Nominal"];

  delete el_medium_pt_map ["Nominal"];
  delete el_medium_eta_map["Nominal"];
  delete el_medium_etas2_map["Nominal"];
  delete el_medium_phi_map["Nominal"];
  delete el_medium_ch_map ["Nominal"];

  delete mu_pt_map["Nominal"];
  delete mu_eta_map["Nominal"];
  delete mu_phi_map["Nominal"];
  delete mu_ch_map["Nominal"];
  delete mu_w_map["Nominal"];
}

TString MiniTree::BookName(TString branch, TString sys_name) 
{
  TString ret = branch + "_" + sys_name;
  return ret;
}

StatusCode MiniTree::initialize() 
{
  //Init the nominal tree
  ph_pt = new std::vector<float>(); 
  ph_eta = new std::vector<float>();
  ph_etas2 = new std::vector<float>();
  ph_phi = new std::vector<float>();
  ph_etcone40 = new std::vector<float>();
  ph_ptcone20 = new std::vector<float>();
  ph_iso = new std::vector<float>();
  ph_trackiso = new std::vector<float>();
  ph_conv = new std::vector<int>();
  ph_w = new std::vector<float>();

  ph_truth_pt = new std::vector<float>(); 
  ph_truth_eta = new std::vector<float>();
  ph_truth_phi = new std::vector<float>();
  ph_truth_id = new std::vector<int>();
  ph_truth_type = new std::vector<int>();
  ph_truth_origin = new std::vector<int>();

  ph_noniso_pt = new std::vector<float>(); 
  ph_noniso_eta = new std::vector<float>();
  ph_noniso_etas2 = new std::vector<float>();
  ph_noniso_phi = new std::vector<float>();
  ph_noniso_etcone40 = new std::vector<float>();
  ph_noniso_ptcone20 = new std::vector<float>();
  ph_noniso_iso = new std::vector<float>();
  ph_noniso_trackiso = new std::vector<float>();
  ph_noniso_conv = new std::vector<int>();
  ph_noniso_w = new std::vector<float>();

  jet_pt = new std::vector<float>(); 
  jet_eta = new std::vector<float>();
  jet_phi = new std::vector<float>();
  jet_e = new std::vector<float>();
  jet_isb = new std::vector<bool>();
  jet_w = new std::vector<float>();
  
  el_pt = new std::vector<float>(); 
  el_eta = new std::vector<float>();
  el_etas2 = new std::vector<float>();
  el_phi = new std::vector<float>();
  el_ch = new std::vector<int>();
  el_w = new std::vector<float>();

  el_medium_pt = new std::vector<float>(); 
  el_medium_eta = new std::vector<float>();
  el_medium_etas2 = new std::vector<float>();
  el_medium_phi = new std::vector<float>();
  el_medium_ch = new std::vector<int>();
  
  mu_pt = new std::vector<float>(); 
  mu_eta = new std::vector<float>();
  mu_phi = new std::vector<float>();
  mu_ch = new std::vector<int>();
  mu_w = new std::vector<float>();

  // Tree  
  TString tree_name = "mini";

  tree = new TTree(tree_name, tree_name);
  tree->SetDirectory(m_outfile);	 

  // Event variables
  tree->Branch("run", &run_number, "run/I");
  tree->Branch("lb", &lumi_block, "lb/I");
  tree->Branch("event", &event_number, "event/l");
  tree->Branch("avgmu", &avg_mu, "avgmu/F");

  if (m_ismc) {
    tree->Branch("fs", &final_state, "fs/i");
    tree->Branch("mcveto",    &mcveto   , "mcveto/i");
  }
  tree->Branch("year", &year, "year/i");

  tree->Branch("pass_g120", &pass_g120, "pass_g120/i");
  tree->Branch("pass_g140", &pass_g140, "pass_g140/i");
  tree->Branch("pass_g70_xe70", &pass_g70_xe70, "pass_g70_xe70/i");

  // Weights
  tree->Branch("weight_mc", &weight_mc, "weight_mc/F"); // no syst
  tree->Branch("weight_sf", &weight_sf_map["Nominal"], "weight_sf/F"); // one for each syst
  tree->Branch("weight_pu", &weight_pu, "weight_pu/F");
  tree->Branch("weight_pu_down", &weight_pu_down, "weight_pu_down/F");
  tree->Branch("weight_pu_up", &weight_pu_up, "weight_pu_up/F");
  tree->Branch("PRWHash", &PRWHash, "PRWHash/l"); // for PURW


  // Nominal blocks
  tree->Branch("ph_n", &ph_n_map["Nominal"], "ph_n/I");
  tree->Branch("ph_pt",  ph_pt);
  tree->Branch("ph_eta", ph_eta);
  tree->Branch("ph_etas2", ph_etas2);
  tree->Branch("ph_phi", ph_phi);
  tree->Branch("ph_etcone40", ph_etcone40);
  tree->Branch("ph_ptcone20", ph_ptcone20);
  tree->Branch("ph_iso", ph_iso);
  tree->Branch("ph_trackiso", ph_trackiso);
  tree->Branch("ph_conv", ph_conv);
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
  tree->Branch("ph_noniso_etas2", ph_noniso_etas2);
  tree->Branch("ph_noniso_phi", ph_noniso_phi);
  tree->Branch("ph_noniso_etcone40", ph_noniso_etcone40);
  tree->Branch("ph_noniso_ptcone20", ph_noniso_ptcone20);
  tree->Branch("ph_noniso_iso", ph_noniso_iso);
  tree->Branch("ph_noniso_trackiso", ph_noniso_trackiso);
  tree->Branch("ph_noniso_conv", ph_noniso_conv);
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
  tree->Branch("el_etas2", el_etas2);
  tree->Branch("el_phi", el_phi);
  tree->Branch("el_pt",  el_pt);
  tree->Branch("el_ch",  el_ch);
  tree->Branch("el_w",   el_w);

  if (m_save_medium_electrons) {
    tree->Branch("el_medium_n", &el_medium_n_map["Nominal"], "el_medium_n/I");
    tree->Branch("el_medium_eta", el_medium_eta);
    tree->Branch("el_medium_etas2", el_medium_etas2);
    tree->Branch("el_medium_phi", el_medium_phi);
    tree->Branch("el_medium_pt",  el_medium_pt);
    tree->Branch("el_medium_ch",  el_medium_ch);
  }

  tree->Branch("mu_n", &mu_n_map["Nominal"], "mu_n/I");
  tree->Branch("mu_eta", mu_eta);
  tree->Branch("mu_phi", mu_phi);
  tree->Branch("mu_pt",  mu_pt);
  tree->Branch("mu_ch",  mu_ch);
  tree->Branch("mu_w",   mu_w);
  
  tree->Branch("met_et",    &met_et_map["Nominal"]);
  tree->Branch("met_phi",   &met_phi_map["Nominal"]);
  tree->Branch("met_sumet", &met_sumet_map["Nominal"]);
  tree->Branch("met_sig",   &met_sig_map["Nominal"]);

  tree->Branch("met_track_et",    &met_track_et_map["Nominal"]);
  tree->Branch("met_track_phi",   &met_track_phi_map["Nominal"]);
  tree->Branch("met_soft_et",  &met_soft_et_map["Nominal"]);
  tree->Branch("met_soft_phi", &met_soft_phi_map["Nominal"]);
  tree->Branch("met_ele_et",   &met_ele_et_map["Nominal"]);
  tree->Branch("met_ele_phi",  &met_ele_phi_map["Nominal"]);
  tree->Branch("met_gam_et",   &met_gam_et_map["Nominal"]);
  tree->Branch("met_gam_phi",  &met_gam_phi_map["Nominal"]);
  tree->Branch("met_muon_et",  &met_muon_et_map["Nominal"]);
  tree->Branch("met_muon_phi", &met_muon_phi_map["Nominal"]);
  tree->Branch("met_jet_et",   &met_jet_et_map["Nominal"]);
  tree->Branch("met_jet_phi",  &met_jet_phi_map["Nominal"]);

  tree->Branch("ht0",   &ht0_map["Nominal"]);				
  tree->Branch("ht",   &ht_map["Nominal"]);				
  tree->Branch("meff", &meff_map["Nominal"]);				
  tree->Branch("rt1",  &rt1_map["Nominal"]);				
  tree->Branch("rt2",  &rt2_map["Nominal"]);				
  tree->Branch("rt3",  &rt3_map["Nominal"]);				
  tree->Branch("rt4",  &rt4_map["Nominal"]);				
  tree->Branch("mt_gam",  &mt_gam_map["Nominal"]);				

  tree->Branch("dphi_gamjet", &dphi_gamjet_map["Nominal"]);
  tree->Branch("dphi_jet1met", &dphi_jet1met_map["Nominal"]);
  tree->Branch("dphi_jetmet", &dphi_jetmet_map["Nominal"]);
  tree->Branch("dphi_gammet", &dphi_gammet_map["Nominal"]);

  tree->Branch("dphi_met_trackmet", &dphi_met_trackmet_map["Nominal"]);

  std::string sys_name = "Nominal"; 
   
  weight_sf_map.insert(std::pair<std::string, float>(sys_name, 1.));
  
  ph_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  ph_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, ph_pt));
  ph_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_eta));
  ph_etas2_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_etas2));
  ph_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_phi));
  ph_etcone40_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_etcone40));
  ph_ptcone20_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_ptcone20));
  ph_iso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_iso));
  ph_trackiso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_trackiso));
  ph_conv_map.insert(std::pair<std::string, std::vector<int>*>(sys_name, ph_conv));
  ph_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_w));

  ph_noniso_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  ph_noniso_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_pt));
  ph_noniso_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_eta));
  ph_noniso_etas2_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_etas2));
  ph_noniso_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_phi));
  ph_noniso_etcone40_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_etcone40));
  ph_noniso_ptcone20_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_ptcone20));
  ph_noniso_iso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_iso));
  ph_noniso_trackiso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_trackiso));
  ph_noniso_conv_map.insert(std::pair<std::string, std::vector<int>*>(sys_name, ph_noniso_conv));
  ph_noniso_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_w));

  jet_n_map.insert  (std::pair<std::string, int>(sys_name, 0));
  bjet_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  jet_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, jet_pt));
  jet_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_eta));
  jet_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_phi));
  jet_e_map.insert  (std::pair<std::string, std::vector<float>*>(sys_name, jet_e));
  jet_isb_map.insert(std::pair<std::string, std::vector<bool>*> (sys_name, jet_isb));
  jet_w_map.insert  (std::pair<std::string, std::vector<float>*>(sys_name, jet_w));

  el_n_map.insert(std::pair<std::string, int>(sys_name, 0));
  el_pt_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_pt));
  el_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_eta));
  el_etas2_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_etas2));
  el_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_phi));
  el_ch_map.insert(std::pair<std::string, std::vector<int>*>  (sys_name, el_ch));
  el_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_w));

  el_medium_n_map.insert (std::pair<std::string, int>(sys_name, 0));
  el_medium_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, el_medium_pt));
  el_medium_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_medium_eta));
  el_medium_etas2_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_medium_etas2));
  el_medium_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_medium_phi));
  el_medium_ch_map.insert (std::pair<std::string, std::vector<int>*>  (sys_name, el_medium_ch));

  mu_n_map.insert (std::pair<std::string, int>(sys_name, 0));  
  mu_pt_map.insert (std::pair<std::string, std::vector<float>*>(sys_name, mu_pt));
  mu_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_eta));
  mu_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_phi));
  mu_ch_map.insert (std::pair<std::string, std::vector<int>*>  (sys_name, mu_ch));
  mu_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_w));
  
  met_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_sumet_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_sig_map.insert(std::pair<std::string, float>(sys_name, -99.));

  met_track_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_track_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_soft_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_soft_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_ele_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_ele_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_gam_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_gam_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_muon_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_muon_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_jet_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
  met_jet_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));

  ht0_map.insert(std::pair<std::string, float>(sys_name, -99.));
  ht_map.insert(std::pair<std::string, float>(sys_name, -99.));
  meff_map.insert(std::pair<std::string, float>(sys_name, -99.));
  rt1_map.insert(std::pair<std::string, float>(sys_name, -99.));
  rt2_map.insert(std::pair<std::string, float>(sys_name, -99.));
  rt3_map.insert(std::pair<std::string, float>(sys_name, -99.));
  rt4_map.insert(std::pair<std::string, float>(sys_name, -99.));
  mt_gam_map.insert(std::pair<std::string, float>(sys_name, -99.));

  dphi_jet1met_map.insert(std::pair<std::string, float>(sys_name, -99.));
  dphi_jetmet_map.insert(std::pair<std::string, float>(sys_name, -99.));
  dphi_gamjet_map.insert(std::pair<std::string, float>(sys_name, -99.));
  dphi_gammet_map.insert(std::pair<std::string, float>(sys_name, -99.));

  dphi_met_trackmet_map.insert(std::pair<std::string, float>(sys_name, -99.));

  // Systematics blocks
  for (const auto& sys : m_sys_list) {

    if (sys.affectsKinematics || sys.affectsWeights ) {

      bool syst_affectsPhotons   = ST::testAffectsObject(xAOD::Type::Photon,   sys.affectsType);
      bool syst_affectsElectrons = ST::testAffectsObject(xAOD::Type::Electron, sys.affectsType);
      bool syst_affectsMuons     = ST::testAffectsObject(xAOD::Type::Muon,     sys.affectsType);
      bool syst_affectsJets      = ST::testAffectsObject(xAOD::Type::Jet,      sys.affectsType);
      bool syst_affectsTaus      = ST::testAffectsObject(xAOD::Type::Tau,      sys.affectsType);
      bool syst_affectsBTag      = ST::testAffectsObject(xAOD::Type::BTag,     sys.affectsType);

      if (syst_affectsTaus)
        continue;
      
      ph_pt = new std::vector<float>(); 
      ph_eta = new std::vector<float>();
      ph_etas2 = new std::vector<float>();
      ph_phi = new std::vector<float>();
      ph_iso = new std::vector<float>();
      ph_trackiso = new std::vector<float>();
      ph_etcone40 = new std::vector<float>();
      ph_ptcone20 = new std::vector<float>();
      ph_conv = new std::vector<int>();
      ph_w = new std::vector<float>();

      ph_noniso_pt = new std::vector<float>(); 
      ph_noniso_eta = new std::vector<float>();
      ph_noniso_etas2 = new std::vector<float>();
      ph_noniso_phi = new std::vector<float>();
      ph_noniso_etcone40 = new std::vector<float>();
      ph_noniso_ptcone20 = new std::vector<float>();
      ph_noniso_iso = new std::vector<float>();
      ph_noniso_trackiso = new std::vector<float>();
      ph_noniso_conv = new std::vector<int>();
      ph_noniso_w = new std::vector<float>();

      jet_pt = new std::vector<float>(); 
      jet_eta = new std::vector<float>();
      jet_phi = new std::vector<float>();
      jet_e = new std::vector<float>();
      jet_isb = new std::vector<bool>();
      jet_w = new std::vector<float>();

      el_pt = new std::vector<float>(); 
      el_eta = new std::vector<float>();
      el_etas2 = new std::vector<float>();
      el_phi = new std::vector<float>();
      el_ch = new std::vector<int>();
      el_w = new std::vector<float>();

      el_medium_pt = new std::vector<float>(); 
      el_medium_eta = new std::vector<float>();
      el_medium_etas2 = new std::vector<float>();
      el_medium_phi = new std::vector<float>();
      el_medium_ch = new std::vector<int>();

      mu_pt = new std::vector<float>(); 
      mu_eta = new std::vector<float>();
      mu_phi = new std::vector<float>();
      mu_ch = new std::vector<int>();
      mu_w = new std::vector<float>();

      const CP::SystematicSet& sysSet = sys.systset;
      sys_name = sysSet.name().c_str();

      // total sf weight
      if (sys.affectsWeights) {
        tree->Branch(BookName("weight_sf", sys_name), &weight_sf_map[sys_name], BookName("weight_sf", sys_name)+"/F");
      }

      // jets/btagging
      if (syst_affectsJets || syst_affectsBTag) {
        
        if (sys.affectsWeights) {
          tree->Branch(BookName("jet_w", sys_name), jet_w);
        }
        else {
          tree->Branch(BookName("jet_n", sys_name), &jet_n_map[sys_name], BookName("jet_n", sys_name)+"/I");
          tree->Branch(BookName("bjet_n", sys_name) , &bjet_n_map[sys_name]);
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
          tree->Branch(BookName("ph_noniso_w", sys_name), ph_noniso_w);
        }
        else {
          tree->Branch(BookName("ph_n", sys_name) , &ph_n_map[sys_name], BookName("ph_n", sys_name)+"/I");
          tree->Branch(BookName("ph_pt", sys_name) , ph_pt);
          tree->Branch(BookName("ph_eta", sys_name), ph_eta);
          tree->Branch(BookName("ph_etas2", sys_name), ph_etas2);
          tree->Branch(BookName("ph_phi", sys_name), ph_phi);
          tree->Branch(BookName("ph_etcone40", sys_name), ph_etcone40);
          tree->Branch(BookName("ph_ptcone20", sys_name), ph_ptcone20);
          tree->Branch(BookName("ph_iso", sys_name), ph_iso);
          tree->Branch(BookName("ph_trackiso", sys_name), ph_trackiso);
          tree->Branch(BookName("ph_conv", sys_name), ph_conv);

          tree->Branch(BookName("ph_noniso_n", sys_name), &ph_noniso_n_map[sys_name], BookName("ph_noniso_n", sys_name)+"/I");
          tree->Branch(BookName("ph_noniso_pt", sys_name) , ph_noniso_pt);
          tree->Branch(BookName("ph_noniso_eta", sys_name), ph_noniso_eta);
          tree->Branch(BookName("ph_noniso_etas2", sys_name), ph_noniso_etas2);
          tree->Branch(BookName("ph_noniso_phi", sys_name), ph_noniso_phi);
          tree->Branch(BookName("ph_noniso_etcone40", sys_name), ph_noniso_etcone40);
          tree->Branch(BookName("ph_noniso_ptcone20", sys_name), ph_noniso_ptcone20);
          tree->Branch(BookName("ph_noniso_iso", sys_name), ph_noniso_iso);
          tree->Branch(BookName("ph_noniso_trackiso", sys_name), ph_noniso_trackiso);
          tree->Branch(BookName("ph_noniso_conv", sys_name), ph_noniso_conv);
	  	}

      }
      
      // electrons
      if (syst_affectsElectrons) {
        if (sys.affectsWeights) {
          tree->Branch(BookName("el_w", sys_name), el_w);
        }
        else {
          tree->Branch(BookName("el_n", sys_name), &el_n_map[sys_name], BookName("el_n", sys_name)+"/I");
          tree->Branch(BookName("el_pt", sys_name), el_pt);
          tree->Branch(BookName("el_eta", sys_name), el_eta);
          tree->Branch(BookName("el_etas2", sys_name), el_etas2);
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
          tree->Branch(BookName("mu_n", sys_name) , &mu_n_map[sys_name], BookName("mu_n", sys_name)+"/I");
          tree->Branch(BookName("mu_pt", sys_name) , mu_pt);
          tree->Branch(BookName("mu_eta", sys_name), mu_eta);
          tree->Branch(BookName("mu_phi", sys_name), mu_phi);
          tree->Branch(BookName("mu_ch", sys_name) , mu_ch);
	  	}
      }
      
      if (sys.affectsKinematics) {
        tree->Branch(BookName("met_et", sys_name), &met_et_map[sys_name]);
        tree->Branch(BookName("met_phi", sys_name), &met_phi_map[sys_name]);
        tree->Branch(BookName("met_sumet", sys_name), &met_sumet_map[sys_name]);
        tree->Branch(BookName("met_sig", sys_name), &met_sig_map[sys_name]);

        tree->Branch(BookName("met_track_et", sys_name),    &met_track_et_map[sys_name]);
        tree->Branch(BookName("met_track_phi", sys_name),   &met_track_phi_map[sys_name]);
        
        tree->Branch(BookName("met_soft_et", sys_name),  &met_soft_et_map[sys_name]);
        tree->Branch(BookName("met_soft_phi", sys_name), &met_soft_phi_map[sys_name]);
        tree->Branch(BookName("met_ele_et", sys_name),   &met_ele_et_map[sys_name]);
        tree->Branch(BookName("met_ele_phi", sys_name),  &met_ele_phi_map[sys_name]);
        tree->Branch(BookName("met_gam_et", sys_name),   &met_gam_et_map[sys_name]);
        tree->Branch(BookName("met_gam_phi", sys_name),  &met_gam_phi_map[sys_name]);
        tree->Branch(BookName("met_muon_et", sys_name),  &met_muon_et_map[sys_name]);
        tree->Branch(BookName("met_muon_phi", sys_name), &met_muon_phi_map[sys_name]);
        tree->Branch(BookName("met_jet_et", sys_name),   &met_jet_et_map[sys_name]);
        tree->Branch(BookName("met_jet_phi", sys_name),  &met_jet_phi_map[sys_name]);
      
        tree->Branch(BookName("ht0", sys_name), &ht0_map[sys_name]);
        tree->Branch(BookName("ht", sys_name), &ht_map[sys_name]);
        tree->Branch(BookName("meff", sys_name), &meff_map[sys_name]);
        tree->Branch(BookName("rt1", sys_name), &rt1_map[sys_name]);
        tree->Branch(BookName("rt2", sys_name), &rt2_map[sys_name]);
        tree->Branch(BookName("rt3", sys_name), &rt3_map[sys_name]);
        tree->Branch(BookName("rt4", sys_name), &rt4_map[sys_name]);
        tree->Branch(BookName("mt_gam", sys_name), &mt_gam_map[sys_name]);
      
        tree->Branch(BookName("dphi_gamjet", sys_name), &dphi_gamjet_map[sys_name]);
        tree->Branch(BookName("dphi_jet1met", sys_name), &dphi_jet1met_map[sys_name]);
        tree->Branch(BookName("dphi_jetmet", sys_name), &dphi_jetmet_map[sys_name]);
        tree->Branch(BookName("dphi_gammet", sys_name), &dphi_gammet_map[sys_name]);
        tree->Branch(BookName("dphi_met_trackmet", sys_name), &dphi_met_trackmet_map[sys_name]);
      }
      
      // tree_map.insert(std::pair<std::string, TTree*>(sys_name, tree));
      ph_n_map.insert(std::pair<std::string, int>(sys_name, 0));
      ph_pt_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_pt));
      ph_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_eta));
      ph_etas2_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_etas2));
      ph_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_phi));
      ph_etcone40_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_etcone40));
      ph_ptcone20_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_ptcone20));
      ph_iso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_iso));
      ph_trackiso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_trackiso));
      ph_conv_map.insert(std::pair<std::string, std::vector<int>*>(sys_name, ph_conv));
      ph_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_w));

      ph_noniso_n_map.insert(std::pair<std::string, int>(sys_name, 0));
      ph_noniso_pt_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_pt));
      ph_noniso_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_eta));
      ph_noniso_etas2_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_etas2));
      ph_noniso_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_phi));
      ph_noniso_etcone40_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_etcone40));
      ph_noniso_ptcone20_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_ptcone20));
      ph_noniso_iso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_iso));
      ph_noniso_trackiso_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_trackiso));
      ph_noniso_conv_map.insert(std::pair<std::string, std::vector<int>*>(sys_name, ph_noniso_conv));
      ph_noniso_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, ph_noniso_w));

      jet_n_map.insert(std::pair<std::string, int>(sys_name, 0));
      bjet_n_map.insert(std::pair<std::string, int>(sys_name, 0));
      jet_pt_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_pt));
      jet_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_eta));
      jet_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_phi));
      jet_e_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_e));
      jet_isb_map.insert(std::pair<std::string, std::vector<bool>*>(sys_name, jet_isb));
      jet_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, jet_w));

      el_n_map.insert(std::pair<std::string, int>(sys_name, 0));
      el_pt_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_pt));
      el_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_eta));
      el_etas2_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_etas2));
      el_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_phi));
      el_ch_map.insert (std::pair<std::string, std::vector<int>*>  (sys_name, el_ch));
      el_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_w));

      el_medium_n_map.insert(std::pair<std::string, int>(sys_name, 0));
      el_medium_pt_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_medium_pt));
      el_medium_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_medium_eta));
      el_medium_etas2_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_medium_etas2));
      el_medium_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, el_medium_phi));
      el_medium_ch_map.insert(std::pair<std::string, std::vector<int>*>  (sys_name, el_medium_ch));

      mu_n_map.insert(std::pair<std::string, int>(sys_name, 0));  
      mu_pt_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_pt));
      mu_eta_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_eta));
      mu_phi_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_phi));
      mu_ch_map.insert(std::pair<std::string, std::vector<int>*>  (sys_name, mu_ch));
      mu_w_map.insert(std::pair<std::string, std::vector<float>*>(sys_name, mu_w));

      met_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_sumet_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_sig_map.insert(std::pair<std::string, float>(sys_name, -99.));

      met_track_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_track_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_soft_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_soft_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_ele_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_ele_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_gam_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_gam_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_muon_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_muon_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_jet_et_map.insert(std::pair<std::string, float>(sys_name, -99.));
      met_jet_phi_map.insert(std::pair<std::string, float>(sys_name, -99.));
      
      ht0_map.insert(std::pair<std::string, float>(sys_name, -99.));
      ht_map.insert(std::pair<std::string, float>(sys_name, -99.));
      meff_map.insert(std::pair<std::string, float>(sys_name, -99.));
      mt_gam_map.insert(std::pair<std::string, float>(sys_name, -99.));

      rt1_map.insert(std::pair<std::string, float>(sys_name, -99.));
      rt2_map.insert(std::pair<std::string, float>(sys_name, -99.));
      rt3_map.insert(std::pair<std::string, float>(sys_name, -99.));
      rt4_map.insert(std::pair<std::string, float>(sys_name, -99.));

      dphi_gamjet_map.insert(std::pair<std::string, float>(sys_name, -99.));
      dphi_jet1met_map.insert(std::pair<std::string, float>(sys_name, -99.));
      dphi_jetmet_map.insert(std::pair<std::string, float>(sys_name, -99.));
      dphi_gammet_map.insert(std::pair<std::string, float>(sys_name, -99.));

      dphi_met_trackmet_map.insert(std::pair<std::string, float>(sys_name, -99.));

      weight_sf_map.insert(std::pair<std::string, float>(sys_name, 1.));
    }
  }

  return StatusCode::SUCCESS;
}

void MiniTree::clear()
{
  run_number = 0;
  lumi_block = 0;
  event_number = 0;
  avg_mu = 0;

  mcveto = 0;
  final_state = 0;
  year = 0;

  pass_g120 = 0;
  pass_g140 = 0;
  pass_g70_xe70 = 0;

  weight_mc = 1.;
  weight_pu = 1.;
  weight_pu_down = 1.;
  weight_pu_up = 1.;

  PRWHash = 0;
}

bool MiniTree::process(AnalysisCollections collections, std::string sysname) 
{
  // clear vectors
  ph_pt_map[sysname]->clear();
  ph_eta_map[sysname]->clear();
  ph_etas2_map[sysname]->clear();
  ph_phi_map[sysname]->clear();
  ph_etcone40_map[sysname]->clear();
  ph_ptcone20_map[sysname]->clear();
  ph_iso_map[sysname]->clear();
  ph_trackiso_map[sysname]->clear();
  ph_conv_map[sysname]->clear();
  ph_w_map[sysname]->clear();

  if (m_ismc) {
    ph_truth_pt->clear();
    ph_truth_eta->clear();
    ph_truth_phi->clear();
    ph_truth_id->clear();
    ph_truth_type->clear();
    ph_truth_origin->clear();
  }

  ph_noniso_pt_map[sysname]->clear();
  ph_noniso_eta_map[sysname]->clear();
  ph_noniso_etas2_map[sysname]->clear();
  ph_noniso_phi_map[sysname]->clear();
  ph_noniso_etcone40_map[sysname]->clear();
  ph_noniso_ptcone20_map[sysname]->clear();
  ph_noniso_iso_map[sysname]->clear();
  ph_noniso_trackiso_map[sysname]->clear();
  ph_noniso_conv_map[sysname]->clear();
  ph_noniso_w_map[sysname]->clear();
  
  jet_pt_map[sysname]->clear();
  jet_eta_map[sysname]->clear();
  jet_phi_map[sysname]->clear();
  jet_e_map[sysname]->clear();
  jet_isb_map[sysname]->clear();
  jet_w_map[sysname]->clear();
  
  el_pt_map[sysname]->clear();
  el_eta_map[sysname]->clear();
  el_etas2_map[sysname]->clear();
  el_phi_map[sysname]->clear();
  el_ch_map[sysname]->clear();
  el_w_map[sysname]->clear();

  if (m_save_medium_electrons && sysname == "Nominal") {
    el_medium_pt_map[sysname]->clear();
    el_medium_eta_map[sysname]->clear();
    el_medium_etas2_map[sysname]->clear();
    el_medium_phi_map[sysname]->clear();
    el_medium_ch_map[sysname]->clear();
  }
  
  mu_pt_map[sysname]->clear();
  mu_eta_map[sysname]->clear();
  mu_phi_map[sysname]->clear();
  mu_ch_map[sysname]->clear();
  mu_w_map[sysname]->clear();

  // clear event variables
  met_phi_map[sysname]   = -99.;
  met_et_map[sysname]     = -99.;
  met_sumet_map[sysname] = -99.;
  met_sig_map[sysname]   = -99.;

  met_track_et_map[sysname] = -99.;
  met_track_phi_map[sysname] = -99.;
  met_soft_et_map[sysname] = -99.;
  met_soft_phi_map[sysname] = -99.;
  met_ele_et_map[sysname] = -99.;
  met_ele_phi_map[sysname] = -99.;
  met_gam_et_map[sysname] = -99.;
  met_gam_phi_map[sysname] = -99.;
  met_muon_et_map[sysname] = -99.;
  met_muon_phi_map[sysname] = -99.;
  met_jet_et_map[sysname] = -99.;
  met_jet_phi_map[sysname] = -99.;

  ht0_map[sysname] = -99.;
  ht_map[sysname] = -99.;
  meff_map[sysname] = -99.;

  rt1_map[sysname] = -99.;
  rt2_map[sysname] = -99.;
  rt3_map[sysname] = -99.;
  rt4_map[sysname] = -99.;

  mt_gam_map[sysname] = -99.;

  dphi_jet1met_map[sysname] = -99.;
  dphi_jetmet_map[sysname] = -99.;
  dphi_gamjet_map[sysname] = -99.;
  dphi_gammet_map[sysname] = -99.;

  dphi_met_trackmet_map[sysname] = -99.;

  // Setting up some basic filtering rule
  collections.photons->setStore(collections.photons_aux);
  collections.electrons->setStore(collections.electrons_aux);
  collections.muons->setStore(collections.muons_aux);
  collections.jets->setStore(collections.jets_aux);
  collections.met->setStore(collections.met_aux);
  collections.met_track->setStore(collections.met_track_aux);


  Double_t total_weight_sf = 1.;

  // electrons
  int el_n = 0;
  for (const auto& el_itr : *collections.electrons) {
    
    if (el_itr->auxdata<char>("baseline") == 1 &&
        el_itr->auxdata<char>("passOR") == 1 &&
        el_itr->auxdata<char>("signal") == 1 && 
        el_itr->auxdata<char>("isol") == 1 && 
        PassEtaCut(el_itr)) {
      
      el_n += 1;
      el_pt_map[sysname]->push_back(el_itr->pt()*IGEV);
      el_eta_map[sysname]->push_back(el_itr->eta());
      el_etas2_map[sysname]->push_back(el_itr->caloCluster()->etaBE(2));
      el_phi_map[sysname]->push_back(el_itr->phi());
      el_ch_map[sysname]->push_back(el_itr->trackParticle()->charge());
      el_w_map[sysname]->push_back( el_itr->auxdata<double>("effscalefact") );

      if (m_ismc)
        total_weight_sf *= el_itr->auxdata<double>("effscalefact");
    }
  }
  el_n_map[sysname] = el_n;

  // medium electrons
  int el_medium_n = 0;
  if (m_save_medium_electrons && sysname == "Nominal") {
    for (const auto& el_itr : *collections.electrons) {

      if (el_itr->auxdata<char>("baseline") == 1 &&
          el_itr->auxdata<char>("passOR") == 1 &&
          el_itr->auxdata<char>("medium") == 1 && 
          el_itr->auxdata<char>("isol") == 1 && 
          PassEtaCut(el_itr)) {
        
        el_medium_n += 1;
        el_medium_pt_map[sysname]->push_back(el_itr->pt()*IGEV);
        el_medium_eta_map[sysname]->push_back(el_itr->eta());
        el_medium_etas2_map[sysname]->push_back(el_itr->caloCluster()->etaBE(2));
        el_medium_phi_map[sysname]->push_back(el_itr->phi());
        el_medium_ch_map[sysname]->push_back(el_itr->trackParticle()->charge());
      }
    }
    el_medium_n_map[sysname] = el_medium_n;
  }

  // muons
  int mu_n = 0;
  for (const auto& mu_itr : *collections.muons) {
    
    if (mu_itr->auxdata<char>("baseline") == 1 &&
        mu_itr->auxdata<char>("passOR") == 1 &&
        mu_itr->auxdata<char>("signal") == 1 &&
        mu_itr->auxdata<char>("isol") == 1 &&
        PassEtaCut(mu_itr)) {

      mu_n += 1;      
      mu_pt_map[sysname] ->push_back(mu_itr->pt()*IGEV);
      mu_eta_map[sysname]->push_back(mu_itr->eta());
      mu_phi_map[sysname]->push_back(mu_itr->phi());
      mu_ch_map[sysname] ->push_back(mu_itr->primaryTrackParticle()->charge());

      float sf = 1.;
      if (m_ismc)
        sf = mu_itr->auxdata<double>("effscalefact");

      mu_w_map[sysname]->push_back(sf);
      total_weight_sf *= sf;

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

      float sf = 1.;
      if (m_ismc)
        sf = jet_itr->auxdata<double>("effscalefact");
      
      jet_w_map[sysname]->push_back(sf);
      total_weight_sf *= sf;

      int isbjet = (jet_itr->auxdata<char>("bjet") == 1);
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
        PassEtaCut(ph_itr, 2.37)) {

      // separate iso and noniso photons
      float etcone40 = ph_itr->isolationValue(xAOD::Iso::topoetcone40)*IGEV;
      float ptcone20 = ph_itr->isolationValue(xAOD::Iso::ptcone20)*IGEV;

      float pt_gev = ph_itr->pt()*IGEV;

      float calo_iso  = etcone40 - 0.022 * pt_gev ; // calo iso without pt dependence
      float track_iso = ptcone20 / pt_gev;          // track iso without pt dependence

      // isolated (only calo iso < 2.45)
      if (ph_itr->auxdata<char>("isol") == 1) {

        ph_n += 1;
        ph_pt_map[sysname] ->push_back(ph_itr->pt()*IGEV);
        ph_eta_map[sysname]->push_back(ph_itr->eta());
        ph_etas2_map[sysname]->push_back(ph_itr->caloCluster()->etaBE(2));
        ph_phi_map[sysname]->push_back(ph_itr->phi());

        ph_etcone40_map[sysname]->push_back(etcone40);
        ph_ptcone20_map[sysname]->push_back(ptcone20);

        ph_iso_map[sysname]      ->push_back(calo_iso);
        ph_trackiso_map[sysname] ->push_back(track_iso);

        // unconverted = 0 : unconverted photon
        // singleSi = 1 : one track only, with Si hits
        // singleTRT = 2 : one track only, no Si hits (TRT only)
        // doubleSi = 3 : two tracks, both with Si hits
        // doubleTRT = 4 : two tracks, none with Si hits (TRT only)
        // doubleSiTRT = 5 : two tracks, only one with Si hits
        ph_conv_map[sysname] ->push_back(ph_itr->conversionType());

        double sf = 1.;
        if (m_ismc)
          sf = ph_itr->auxdata<double>("effscalefact");
        
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

          ph_truth_type  ->push_back(xAOD::TruthHelpers::getParticleTruthType(*ph_itr));
          ph_truth_origin->push_back(xAOD::TruthHelpers::getParticleTruthOrigin(*ph_itr));
        }
      }

      // noniso photons
      else {

        ph_noniso_n += 1;

        ph_noniso_pt_map[sysname]   ->push_back(ph_itr->pt()*IGEV);
        ph_noniso_eta_map[sysname]  ->push_back(ph_itr->eta());
        ph_noniso_etas2_map[sysname]->push_back(ph_itr->caloCluster()->etaBE(2));
        ph_noniso_phi_map[sysname]  ->push_back(ph_itr->phi());

        ph_noniso_etcone40_map[sysname]->push_back(etcone40);
        ph_noniso_ptcone20_map[sysname]->push_back(ptcone20);

        ph_noniso_iso_map[sysname]      ->push_back(calo_iso);
        ph_noniso_trackiso_map[sysname] ->push_back(track_iso);
        
        ph_noniso_conv_map[sysname] ->push_back(ph_itr->conversionType());
        
        double sf = 1.;
        if (m_ismc)
          sf = ph_itr->auxdata<double>("effscalefact");

        ph_noniso_w_map[sysname]->push_back(sf);
      }
      
    }
  }
  ph_n_map[sysname] = ph_n;
  ph_noniso_n_map[sysname] = ph_noniso_n;

  int ph_skim = ph_n + ph_noniso_n;

  // SF weigth
  weight_sf_map[sysname] = total_weight_sf;

  // met
  xAOD::MissingETContainer::const_iterator met_it = collections.met->find("Final");
  if (met_it == collections.met->end()) {
    Error("PhotonMetNtuple:MiniTree", "No RefFinal inside MET container");
  }
  met_et_map[sysname]    = (*met_it)->met() * IGEV; 
  met_phi_map[sysname]   = (*met_it)->phi(); 
  met_sumet_map[sysname] = (*met_it)->sumet() * IGEV; 
  met_sig_map[sysname]   = ((*met_it)->met() * IGEV) / sqrt((*met_it)->sumet() * IGEV); 

  xAOD::MissingETContainer::const_iterator tst_it = collections.met->find("PVSoftTrk");
  if (tst_it == collections.met->end()) {
    Error("PhotonMetNtuple:MiniTree", "No PVSoftTrk inside MET container");
  }
  met_soft_et_map[sysname]    = (*tst_it)->met() * IGEV;
  met_soft_phi_map[sysname]   = (*tst_it)->phi();

  xAOD::MissingETContainer::const_iterator refele_it = collections.met->find("RefEle");
  xAOD::MissingETContainer::const_iterator refgam_it = collections.met->find("RefGamma");
  xAOD::MissingETContainer::const_iterator refmuon_it = collections.met->find("Muons");
  xAOD::MissingETContainer::const_iterator refjet_it = collections.met->find("RefJet");

  met_ele_et_map[sysname]    = (*refele_it)->met() * IGEV;
  met_ele_phi_map[sysname]   = (*refele_it)->phi();
  met_gam_et_map[sysname]    = (*refgam_it)->met() * IGEV;
  met_gam_phi_map[sysname]   = (*refgam_it)->phi();
  met_muon_et_map[sysname]   = (*refmuon_it)->met() * IGEV;
  met_muon_phi_map[sysname]  = (*refmuon_it)->phi();
  met_jet_et_map[sysname]    = (*refjet_it)->met() * IGEV;
  met_jet_phi_map[sysname]   = (*refjet_it)->phi();

  xAOD::MissingETContainer::const_iterator met_track_it = collections.met_track->find("Track");
  if (met_track_it == collections.met_track->end()) {
    Error("PhotonMetNtuple:MiniTree", "No RefFinal inside TrackMET container");
  }
  else {
    met_track_et_map[sysname]    = (*met_track_it)->met() * IGEV; 
    met_track_phi_map[sysname]   = (*met_track_it)->phi(); 
  }


  // Compute extra variables
  Double_t sum_jet_pt = 0.0;
  Double_t sum_jet1_pt = 0.0;
  Double_t sum_jet2_pt = 0.0;
  Double_t sum_jet3_pt = 0.0;
  Double_t sum_jet4_pt = 0.0;

  for (int i=0; i<jet_n; i++) {

    Double_t jetpt = (*jet_pt_map[sysname])[i];

    sum_jet_pt += jetpt;

    if (i < 1) 
      sum_jet1_pt += jetpt;
    if (i < 2) 
      sum_jet2_pt += jetpt;
    if (i < 3) 
      sum_jet3_pt += jetpt;
    if (i < 4) 
      sum_jet4_pt += jetpt;
  }
  
  // Ht
  Double_t ht = sum_jet_pt;
  if (ph_n > 0)
    ht += (*ph_pt_map[sysname])[0];

  ht0_map[sysname] = sum_jet_pt;
  ht_map[sysname] = ht;

  // Meff
  meff_map[sysname] = ht + met_et_map[sysname];
  
  // Rt
  if (jet_n > 0) {
    rt1_map[sysname] = sum_jet1_pt/sum_jet_pt;
    rt2_map[sysname] = sum_jet2_pt/sum_jet_pt;
    rt3_map[sysname] = sum_jet3_pt/sum_jet_pt;
    rt4_map[sysname] = sum_jet4_pt/sum_jet_pt;
  }

  // min dphi between met and the first two jets
  Double_t dphi1 = 4.;
  Double_t dphi2 = 4.;
  if (jet_n > 0) dphi1 = get_dphi((*jet_phi_map[sysname])[0], met_phi_map[sysname]);
  if (jet_n > 1) dphi2 = get_dphi((*jet_phi_map[sysname])[1], met_phi_map[sysname]);
  
  if (jet_n > 0) {
    dphi_jet1met_map[sysname] = dphi1;
    dphi_jetmet_map[sysname] = TMath::Min(dphi1, dphi2);
  }

  // dphi between leading photon and leading jet
  if (ph_n > 0 && jet_n > 0) 
    dphi_gamjet_map[sysname] = get_dphi((*ph_phi_map[sysname])[0], (*jet_phi_map[sysname])[0]);
  
  // dphi between leading photon and MET
  if (ph_n > 0) {
    dphi_gammet_map[sysname] = get_dphi((*ph_phi_map[sysname])[0], met_phi_map[sysname]);

    // MT
    Float_t mt2_gam = 2 * met_et_map[sysname] * (*ph_pt_map[sysname])[0] * (1 - TMath::Cos(dphi_gammet_map[sysname]));
    mt_gam_map[sysname] = TMath::Sqrt(mt2_gam);
  }

  // Add dphi met track
  dphi_met_trackmet_map[sysname] = get_dphi(met_phi_map[sysname], met_track_phi_map[sysname]);

  // Skim: at least one signal photon (no iso cut) with pt>75 (or a medium electron for Nominal)
  bool pass_skim = false;
  if (sysname == "Nominal")
    pass_skim = (ph_skim > 0 || el_medium_n > 0);
  else
    pass_skim = (ph_skim > 0);

  return pass_skim;
}

// Call after all the syst have been processed and ONLY IF one of them passed the event selection criteria
StatusCode MiniTree::FillTree()
{
  tree->Fill();
  return StatusCode::SUCCESS;
}

bool MiniTree::PassEtaCut(const xAOD::IParticle *part, Float_t maxeta)
{
  Double_t eta = fabs(part->eta());

  if (eta > maxeta)
    return false;

  return true;
}

bool MiniTree::PassEtaCut(const xAOD::Photon *part, Float_t maxeta) 
{
  Double_t eta = fabs(part->caloCluster()->etaBE(2));

  if (eta > maxeta)
    return false;

  if (eta > 1.37 && eta < 1.52)
    return false;

  return true;
}

bool MiniTree::PassEtaCut(const xAOD::Electron *part, Float_t maxeta) 
{
  Double_t eta = fabs(part->caloCluster()->etaBE(2));

  if (eta > maxeta)
    return false;

  if (eta > 1.37 && eta < 1.52)
    return false;

  return true;
}

