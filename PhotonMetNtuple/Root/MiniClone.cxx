#define mini_cxx

#include <PhotonMetNtuple/MiniClone.h>

#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

MiniClone::MiniClone(TString orig_path, TString clone_path) : 
  orig_tree(0),
  clone_tree(0),
  m_ismc(false),
  m_efake_sample(false),
  m_jfake_sample(false),
  m_is_old_version(false)
{
  orig_file = new TFile(orig_path);
  orig_file->GetObject("mini", orig_tree);

  clone_file = new TFile(clone_path, "recreate");
  clone_tree = new TTree("mini", "mini");
}

MiniClone::~MiniClone()
{
  if (orig_tree) delete orig_tree->GetCurrentFile();
  if (clone_tree) delete clone_tree->GetCurrentFile();

}

void MiniClone::Init()
{
  InitOriginalTree();
  CreateCloneTree();
}

Int_t MiniClone::GetEntry(Long64_t entry)
{
  // Read contents of entry.
   if (!orig_tree) return 0;
   return orig_tree->GetEntry(entry);
}

void MiniClone::InitOriginalTree()
{
   // Set object pointer
   ph_pt = 0;
   ph_eta = 0;
   ph_etas2 = 0;
   ph_phi = 0;
   ph_etcone40 = 0;
   ph_ptcone20 = 0;
   ph_iso = 0;
   ph_trackiso = 0;
   ph_conv = 0;
   ph_w = 0;
   ph_truth_pt = 0;
   ph_truth_eta = 0;
   ph_truth_phi = 0;
   ph_truth_id = 0;
   ph_truth_type = 0;
   ph_truth_origin = 0;
   ph_noniso_pt = 0;
   ph_noniso_eta = 0;
   ph_noniso_etas2 = 0;
   ph_noniso_phi = 0;
   ph_noniso_etcone40 = 0;
   ph_noniso_ptcone20 = 0;
   ph_noniso_iso = 0;
   ph_noniso_trackiso = 0;
   ph_noniso_conv = 0;
   ph_noniso_w = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_pt = 0;
   jet_e = 0;
   jet_isb = 0;
   jet_w = 0;
   el_eta = 0;
   el_etas2 = 0;
   el_phi = 0;
   el_pt = 0;
   el_ch = 0;
   el_w = 0;
   el_medium_eta = 0;
   el_medium_etas2 = 0;
   el_medium_phi = 0;
   el_medium_pt = 0;
   el_medium_ch = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_pt = 0;
   mu_ch = 0;
   mu_w = 0;

   // Set branch addresses and branch pointers
   orig_tree->SetBranchAddress("run", &run, &b_run);
   orig_tree->SetBranchAddress("lb", &lb, &b_lb);
   orig_tree->SetBranchAddress("event", &event, &b_event);
   orig_tree->SetBranchAddress("avgmu", &avgmu, &b_avgmu);

   orig_tree->SetBranchAddress("year", &year, &b_year);
   orig_tree->SetBranchAddress("pass_g120", &pass_g120, &b_pass_g120);
   orig_tree->SetBranchAddress("pass_g140", &pass_g140, &b_pass_g140);

   if (orig_tree->GetListOfBranches()->FindObject("pass_g0_xe70"))
     orig_tree->SetBranchAddress("pass_g0_xe70", &pass_g70_xe70, &b_pass_g70_xe70);
   else
     orig_tree->SetBranchAddress("pass_g70_xe70", &pass_g70_xe70, &b_pass_g70_xe70);

   if (m_ismc) orig_tree->SetBranchAddress("mcveto", &mcveto, &b_mcveto);
   if (m_ismc) orig_tree->SetBranchAddress("fs", &fs, &b_final_state);

   orig_tree->SetBranchAddress("weight_mc", &weight_mc, &b_weight_mc);
   orig_tree->SetBranchAddress("weight_sf", &weight_sf, &b_weight_sf);
   orig_tree->SetBranchAddress("weight_pu", &weight_pu, &b_weight_pu);
   orig_tree->SetBranchAddress("weight_pu_down", &weight_pu_down, &b_weight_pu_down);
   orig_tree->SetBranchAddress("weight_pu_up", &weight_pu_up, &b_weight_pu_up);
   orig_tree->SetBranchAddress("PRWHash", &PRWHash, &b_PRWHash);

   if (m_efake_sample) {
     orig_tree->SetBranchAddress("weight_feg", &weight_feg, &b_weight_feg);
     orig_tree->SetBranchAddress("weight_feg_dn", &weight_feg_dn, &b_weight_feg_dn);
     orig_tree->SetBranchAddress("weight_feg_up", &weight_feg_up, &b_weight_feg_up);
   }
   if (m_jfake_sample) {
     orig_tree->SetBranchAddress("weight_fjg", &weight_fjg, &b_weight_fjg);
     orig_tree->SetBranchAddress("weight_fjg_dn", &weight_fjg_dn, &b_weight_fjg_dn);
     orig_tree->SetBranchAddress("weight_fjg_up", &weight_fjg_up, &b_weight_fjg_up);
   }

   orig_tree->SetBranchAddress("ph_n", &ph_n, &b_ph_n);
   orig_tree->SetBranchAddress("ph_pt", &ph_pt, &b_ph_pt);
   orig_tree->SetBranchAddress("ph_eta", &ph_eta, &b_ph_eta);
   orig_tree->SetBranchAddress("ph_etas2", &ph_etas2, &b_ph_etas2);
   orig_tree->SetBranchAddress("ph_phi", &ph_phi, &b_ph_phi);
   orig_tree->SetBranchAddress("ph_etcone40", &ph_etcone40, &b_ph_etcone40);
   orig_tree->SetBranchAddress("ph_ptcone20", &ph_ptcone20, &b_ph_ptcone20);
   orig_tree->SetBranchAddress("ph_iso", &ph_iso, &b_ph_iso);
   orig_tree->SetBranchAddress("ph_trackiso", &ph_trackiso, &b_ph_trackiso);
   orig_tree->SetBranchAddress("ph_conv", &ph_conv, &b_ph_conv);
   orig_tree->SetBranchAddress("ph_w", &ph_w, &b_ph_w);

   if (m_ismc) {
     orig_tree->SetBranchAddress("ph_truth_pt", &ph_truth_pt, &b_ph_truth_pt);
     orig_tree->SetBranchAddress("ph_truth_eta", &ph_truth_eta, &b_ph_truth_eta);
     orig_tree->SetBranchAddress("ph_truth_phi", &ph_truth_phi, &b_ph_truth_phi);
     orig_tree->SetBranchAddress("ph_truth_id", &ph_truth_id, &b_ph_truth_id);
     orig_tree->SetBranchAddress("ph_truth_type", &ph_truth_type, &b_ph_truth_type);
     orig_tree->SetBranchAddress("ph_truth_origin", &ph_truth_origin, &b_ph_truth_origin);
   }

   orig_tree->SetBranchAddress("ph_noniso_n", &ph_noniso_n, &b_ph_noniso_n);
   orig_tree->SetBranchAddress("ph_noniso_pt", &ph_noniso_pt, &b_ph_noniso_pt);
   orig_tree->SetBranchAddress("ph_noniso_eta", &ph_noniso_eta, &b_ph_noniso_eta);
   orig_tree->SetBranchAddress("ph_noniso_etas2", &ph_noniso_etas2, &b_ph_noniso_etas2);
   orig_tree->SetBranchAddress("ph_noniso_phi", &ph_noniso_phi, &b_ph_noniso_phi);
   orig_tree->SetBranchAddress("ph_noniso_etcone40", &ph_noniso_etcone40, &b_ph_noniso_etcone40);
   orig_tree->SetBranchAddress("ph_noniso_ptcone20", &ph_noniso_ptcone20, &b_ph_noniso_ptcone20);
   orig_tree->SetBranchAddress("ph_noniso_iso", &ph_noniso_iso, &b_ph_noniso_iso);
   orig_tree->SetBranchAddress("ph_noniso_trackiso", &ph_noniso_trackiso, &b_ph_noniso_trackiso);
   orig_tree->SetBranchAddress("ph_noniso_conv", &ph_noniso_conv, &b_ph_noniso_conv);
   orig_tree->SetBranchAddress("ph_noniso_w", &ph_noniso_w, &b_ph_noniso_w);

   orig_tree->SetBranchAddress("jet_n", &jet_n, &b_jet_n);
   orig_tree->SetBranchAddress("bjet_n", &bjet_n, &b_bjet_n);
   orig_tree->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   orig_tree->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   orig_tree->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   orig_tree->SetBranchAddress("jet_e", &jet_e, &b_jet_e);
   orig_tree->SetBranchAddress("jet_isb", &jet_isb, &b_jet_isb);
   orig_tree->SetBranchAddress("jet_w", &jet_w, &b_jet_w);

   orig_tree->SetBranchAddress("el_n", &el_n, &b_el_n);
   orig_tree->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   orig_tree->SetBranchAddress("el_etas2", &el_etas2, &b_el_etas2);
   orig_tree->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   orig_tree->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
   orig_tree->SetBranchAddress("el_ch", &el_ch, &b_el_ch);
   orig_tree->SetBranchAddress("el_w", &el_w, &b_el_w);

   if (orig_tree->GetListOfBranches()->FindObject("el_medium_n")) {
     orig_tree->SetBranchAddress("el_medium_n", &el_medium_n, &b_el_medium_n);
     orig_tree->SetBranchAddress("el_medium_eta", &el_medium_eta, &b_el_medium_eta);
     orig_tree->SetBranchAddress("el_medium_etas2", &el_medium_etas2, &b_el_medium_etas2);
     orig_tree->SetBranchAddress("el_medium_phi", &el_medium_phi, &b_el_medium_phi);
     orig_tree->SetBranchAddress("el_medium_pt", &el_medium_pt, &b_el_medium_pt);
     orig_tree->SetBranchAddress("el_medium_ch", &el_medium_ch, &b_el_medium_ch);
   }

   orig_tree->SetBranchAddress("mu_n", &mu_n, &b_mu_n);
   orig_tree->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   orig_tree->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   orig_tree->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
   orig_tree->SetBranchAddress("mu_ch", &mu_ch, &b_mu_ch);
   orig_tree->SetBranchAddress("mu_w", &mu_w, &b_mu_w);

   orig_tree->SetBranchAddress("met_et", &met_et, &b_met_et);
   orig_tree->SetBranchAddress("met_phi", &met_phi, &b_met_phi);

   if (orig_tree->GetListOfBranches()->FindObject("met_sumet")) {
     orig_tree->SetBranchAddress("met_sumet", &met_sumet, &b_met_sumet);
     orig_tree->SetBranchAddress("met_sig", &met_sig, &b_met_sig);
   }
   
   if (orig_tree->GetListOfBranches()->FindObject("met_track_et")) {
     orig_tree->SetBranchAddress("met_track_et", &met_track_et, &b_met_track_et);
     orig_tree->SetBranchAddress("met_track_phi", &met_track_phi, &b_met_track_phi);

     orig_tree->SetBranchAddress("met_soft_et", &met_soft_et, &b_met_soft_et);
     orig_tree->SetBranchAddress("met_soft_phi", &met_soft_phi, &b_met_soft_phi);
     orig_tree->SetBranchAddress("met_ele_et", &met_ele_et, &b_met_ele_et);
     orig_tree->SetBranchAddress("met_ele_phi", &met_ele_phi, &b_met_ele_phi);
     orig_tree->SetBranchAddress("met_gam_et", &met_gam_et, &b_met_gam_et);
     orig_tree->SetBranchAddress("met_gam_phi", &met_gam_phi, &b_met_gam_phi);
     orig_tree->SetBranchAddress("met_muon_et", &met_muon_et, &b_met_muon_et);
     orig_tree->SetBranchAddress("met_muon_phi", &met_muon_phi, &b_met_muon_phi);
     orig_tree->SetBranchAddress("met_jet_et", &met_jet_et, &b_met_jet_et);
     orig_tree->SetBranchAddress("met_jet_phi", &met_jet_phi, &b_met_jet_phi);
   }
   else
     m_is_old_version = true;

   orig_tree->SetBranchAddress("ht0", &ht0, &b_ht0);
   orig_tree->SetBranchAddress("ht", &ht, &b_ht);
   orig_tree->SetBranchAddress("meff", &meff, &b_meff);
   orig_tree->SetBranchAddress("rt1", &rt1, &b_rt1);
   orig_tree->SetBranchAddress("rt2", &rt2, &b_rt2);
   orig_tree->SetBranchAddress("rt3", &rt3, &b_rt3);
   orig_tree->SetBranchAddress("rt4", &rt4, &b_rt4);
   orig_tree->SetBranchAddress("mt_gam", &mt_gam, &b_mt_gam);
 
   orig_tree->SetBranchAddress("dphi_gamjet", &dphi_gamjet, &b_dphi_gamjet);
   orig_tree->SetBranchAddress("dphi_jetmet", &dphi_jetmet, &b_dphi_jetmet);
   orig_tree->SetBranchAddress("dphi_jet1met", &dphi_jet1met, &b_dphi_jet1met);
   orig_tree->SetBranchAddress("dphi_gammet", &dphi_gammet, &b_dphi_gammet);

   orig_tree->SetBranchAddress("dphi_met_trackmet", &dphi_met_trackmet, &b_dphi_met_trackmet);
}

void MiniClone::CreateCloneTree()
{
  new_ph_pt    = new std::vector<float>(); 
  new_ph_eta   = new std::vector<float>();
  new_ph_etas2 = new std::vector<float>();
  new_ph_phi   = new std::vector<float>();
  new_ph_etcone40   = new std::vector<float>();
  new_ph_ptcone20 = new std::vector<float>();
  new_ph_iso = new std::vector<float>();
  new_ph_trackiso = new std::vector<float>();
  new_ph_conv = new std::vector<int>();
  new_ph_w     = new std::vector<float>();

  new_ph_truth_pt    = new std::vector<float>(); 
  new_ph_truth_eta   = new std::vector<float>();
  new_ph_truth_phi   = new std::vector<float>();
  new_ph_truth_id   = new std::vector<int>();
  new_ph_truth_type = new std::vector<int>();
  new_ph_truth_origin = new std::vector<int>();

  new_ph_noniso_pt    = new std::vector<float>(); 
  new_ph_noniso_eta   = new std::vector<float>();
  new_ph_noniso_etas2   = new std::vector<float>();
  new_ph_noniso_phi   = new std::vector<float>();
  new_ph_noniso_etcone40   = new std::vector<float>();
  new_ph_noniso_ptcone20   = new std::vector<float>();
  new_ph_noniso_iso   = new std::vector<float>();
  new_ph_noniso_trackiso   = new std::vector<float>();
  new_ph_noniso_conv = new std::vector<int>();
  new_ph_noniso_w     = new std::vector<float>();

  new_jet_pt  = new std::vector<float>(); 
  new_jet_eta = new std::vector<float>();
  new_jet_phi = new std::vector<float>();
  new_jet_e   = new std::vector<float>();
  new_jet_isb = new std::vector<bool>();
  new_jet_w   = new std::vector<float>();
  
  new_el_pt   = new std::vector<float>(); 
  new_el_eta  = new std::vector<float>();
  new_el_etas2 = new std::vector<float>();
  new_el_phi  = new std::vector<float>();
  new_el_ch   = new std::vector<int>();
  new_el_w    = new std::vector<float>();

  new_el_medium_pt   = new std::vector<float>(); 
  new_el_medium_eta  = new std::vector<float>();
  new_el_medium_etas2 = new std::vector<float>();
  new_el_medium_phi  = new std::vector<float>();
  new_el_medium_ch   = new std::vector<int>();
  
  new_mu_pt  = new std::vector<float>(); 
  new_mu_eta = new std::vector<float>();
  new_mu_phi = new std::vector<float>();
  new_mu_ch  = new std::vector<int>();
  new_mu_w   = new std::vector<float>();

  clone_tree->Branch("run", &new_run, "run/I");
  clone_tree->Branch("lb", &new_lb, "lb/I");
  clone_tree->Branch("event", &new_event, "event/l");
  clone_tree->Branch("avgmu", &new_avgmu, "avgmu/F");

  if (m_ismc) clone_tree->Branch("fs", &new_fs, "fs/i");

  //clone_tree->Branch("pass_tst_cleaning", &new_pass_tst_cleaning, "pass_tst_cleaning/i");
  clone_tree->Branch("year", &new_year, "year/i");
  if (m_ismc) clone_tree->Branch("mcveto", &new_mcveto, "mcveto/i");
  clone_tree->Branch("pass_g120", &new_pass_g120, "pass_g120/i");
  clone_tree->Branch("pass_g140", &new_pass_g140, "pass_g140/i");
  clone_tree->Branch("pass_g70_xe70", &new_pass_g70_xe70, "pass_g70_xe70/i");

  clone_tree->Branch("weight_mc", &new_weight_mc, "weight_mc/F");
  clone_tree->Branch("weight_pu", &new_weight_pu, "weight_pu/F");
  clone_tree->Branch("weight_sf", &new_weight_sf, "weight_sf/F");
  clone_tree->Branch("weight_pu_down", &new_weight_pu_down, "weight_pu_down/F");
  clone_tree->Branch("weight_pu_up", &new_weight_pu_up, "weight_pu_up/F");
  clone_tree->Branch("PRWHash", &new_PRWHash, "PRWHash/l");

  if (m_efake_sample) {
    clone_tree->Branch("weight_feg",    &new_weight_feg,    "weight_mc/F");
    clone_tree->Branch("weight_feg_dn", &new_weight_feg_dn, "weight_mc_dn/F");
    clone_tree->Branch("weight_feg_up", &new_weight_feg_up, "weight_mc_up/F");    
  }

  if (m_jfake_sample) {
    clone_tree->Branch("weight_fjg",    &new_weight_fjg,    "weight_mc/F");
    clone_tree->Branch("weight_fjg_dn", &new_weight_fjg_dn, "weight_mc_dn/F");
    clone_tree->Branch("weight_fjg_up", &new_weight_fjg_up, "weight_mc_up/F");    
  }

  clone_tree->Branch("ph_n", &new_ph_n, "ph_n/I");
  clone_tree->Branch("ph_pt",  new_ph_pt);
  clone_tree->Branch("ph_eta", new_ph_eta);
  clone_tree->Branch("ph_etas2", new_ph_etas2);
  clone_tree->Branch("ph_phi", new_ph_phi);
  clone_tree->Branch("ph_etcone40", new_ph_etcone40);
  clone_tree->Branch("ph_ptcone20", new_ph_ptcone20);
  clone_tree->Branch("ph_iso", new_ph_iso);
  clone_tree->Branch("ph_trackiso", new_ph_trackiso);
  clone_tree->Branch("ph_conv", new_ph_conv);
  clone_tree->Branch("ph_w",   new_ph_w);

  if (m_ismc) {
    clone_tree->Branch("ph_truth_pt",  new_ph_truth_pt);
    clone_tree->Branch("ph_truth_eta", new_ph_truth_eta);
    clone_tree->Branch("ph_truth_phi", new_ph_truth_phi);
    clone_tree->Branch("ph_truth_id",  new_ph_truth_id);
    clone_tree->Branch("ph_truth_type", new_ph_truth_type);
    clone_tree->Branch("ph_truth_origin", new_ph_truth_origin);
  }

  clone_tree->Branch("ph_noniso_n", &new_ph_noniso_n, "ph_noniso_n/I");
  clone_tree->Branch("ph_noniso_pt",  new_ph_noniso_pt);
  clone_tree->Branch("ph_noniso_eta", new_ph_noniso_eta);
  clone_tree->Branch("ph_noniso_etas2", new_ph_noniso_etas2);
  clone_tree->Branch("ph_noniso_phi", new_ph_noniso_phi);
  clone_tree->Branch("ph_noniso_etcone40", new_ph_noniso_etcone40);
  clone_tree->Branch("ph_noniso_ptcone20", new_ph_noniso_ptcone20);
  clone_tree->Branch("ph_noniso_iso", new_ph_noniso_iso);
  clone_tree->Branch("ph_noniso_trackiso", new_ph_noniso_trackiso);
  clone_tree->Branch("ph_noniso_conv", new_ph_noniso_conv);
  clone_tree->Branch("ph_noniso_w",   new_ph_noniso_w);

  clone_tree->Branch("jet_n", &new_jet_n, "jet_n/I");
  clone_tree->Branch("bjet_n", &new_bjet_n, "bjet_n/I");
  clone_tree->Branch("jet_eta", new_jet_eta);
  clone_tree->Branch("jet_phi", new_jet_phi);
  clone_tree->Branch("jet_pt",  new_jet_pt);
  clone_tree->Branch("jet_e",   new_jet_e);
  clone_tree->Branch("jet_isb", new_jet_isb);
  clone_tree->Branch("jet_w",  new_jet_w);
  
  clone_tree->Branch("el_n", &new_el_n, "el_n/I");
  clone_tree->Branch("el_eta", new_el_eta);
  clone_tree->Branch("el_etas2", new_el_etas2);
  clone_tree->Branch("el_phi", new_el_phi);
  clone_tree->Branch("el_pt",  new_el_pt);
  clone_tree->Branch("el_ch",  new_el_ch);
  clone_tree->Branch("el_w",   new_el_w);

  clone_tree->Branch("el_medium_n", &new_el_medium_n, "el_medium_n/I");
  clone_tree->Branch("el_medium_eta", new_el_medium_eta);
  clone_tree->Branch("el_medium_etas2", new_el_medium_etas2);
  clone_tree->Branch("el_medium_phi", new_el_medium_phi);
  clone_tree->Branch("el_medium_pt",  new_el_medium_pt);
  clone_tree->Branch("el_medium_ch",  new_el_medium_ch);
  
  clone_tree->Branch("mu_n", &new_mu_n, "mu_n/I");
  clone_tree->Branch("mu_eta", new_mu_eta);
  clone_tree->Branch("mu_phi", new_mu_phi);
  clone_tree->Branch("mu_pt",  new_mu_pt);
  clone_tree->Branch("mu_ch",  new_mu_ch);
  clone_tree->Branch("mu_w",   new_mu_w);
  
  clone_tree->Branch("met_et",    &new_met_et);
  clone_tree->Branch("met_phi",   &new_met_phi);
  clone_tree->Branch("met_sumet", &new_met_sumet);
  clone_tree->Branch("met_sig",   &new_met_sig);

  clone_tree->Branch("met_track_et",   &new_met_track_et);
  clone_tree->Branch("met_track_phi",   &new_met_track_phi);
  clone_tree->Branch("met_soft_et",   &new_met_soft_et);
  clone_tree->Branch("met_soft_phi",   &new_met_soft_phi);
  clone_tree->Branch("met_ele_et",   &new_met_ele_et);
  clone_tree->Branch("met_ele_phi",   &new_met_ele_phi);
  clone_tree->Branch("met_gam_et",   &new_met_gam_et);
  clone_tree->Branch("met_gam_phi",   &new_met_gam_phi);
  clone_tree->Branch("met_muon_et",   &new_met_muon_et);
  clone_tree->Branch("met_muon_phi",   &new_met_muon_phi);
  clone_tree->Branch("met_jet_et",   &new_met_jet_et);
  clone_tree->Branch("met_jet_phi",   &new_met_jet_phi);

  clone_tree->Branch("ht0", &new_ht0);				
  clone_tree->Branch("ht", &new_ht);				
  clone_tree->Branch("meff", &new_meff);				

  clone_tree->Branch("rt1", &new_rt1);				
  clone_tree->Branch("rt2", &new_rt2);				
  clone_tree->Branch("rt3", &new_rt3);				
  clone_tree->Branch("rt4", &new_rt4);				

  clone_tree->Branch("mt_gam", &new_mt_gam);				

  clone_tree->Branch("dphi_gamjet", &new_dphi_gamjet);
  clone_tree->Branch("dphi_jetmet", &new_dphi_jetmet);
  clone_tree->Branch("dphi_jet1met", &new_dphi_jet1met);
  clone_tree->Branch("dphi_gammet", &new_dphi_gammet);

  clone_tree->Branch("dphi_met_trackmet", &new_dphi_met_trackmet);

}

void MiniClone::Clear()
{
  new_ph_pt->clear();
  new_ph_eta->clear();
  new_ph_etas2->clear();
  new_ph_phi->clear();
  new_ph_etcone40->clear();
  new_ph_ptcone20->clear();
  new_ph_iso->clear();
  new_ph_trackiso->clear();
  new_ph_conv->clear();
  new_ph_w->clear();
  
  new_ph_truth_pt->clear();
  new_ph_truth_eta->clear();
  new_ph_truth_phi->clear();
  new_ph_truth_id->clear();
  new_ph_truth_type->clear();
  new_ph_truth_origin->clear();

  new_ph_noniso_pt->clear();
  new_ph_noniso_eta->clear();
  new_ph_noniso_etas2->clear();
  new_ph_noniso_phi->clear();
  new_ph_noniso_etcone40->clear();
  new_ph_noniso_ptcone20->clear();
  new_ph_noniso_iso->clear();
  new_ph_noniso_trackiso->clear();
  new_ph_noniso_conv->clear();
  new_ph_noniso_w->clear();

  new_jet_pt->clear();
  new_jet_eta->clear();
  new_jet_phi->clear();
  new_jet_e->clear();
  new_jet_isb->clear();
  new_jet_w->clear();
    
  new_el_pt->clear();
  new_el_eta->clear();
  new_el_etas2->clear();
  new_el_phi->clear();
  new_el_ch->clear();
  new_el_w->clear();

  new_el_medium_pt->clear();
  new_el_medium_eta->clear();
  new_el_medium_etas2->clear();
  new_el_medium_phi->clear();
  new_el_medium_ch->clear();
  
  new_mu_pt->clear();
  new_mu_eta->clear();
  new_mu_phi->clear();
  new_mu_ch->clear();
  new_mu_w->clear();

}

void MiniClone::CopyAllBlocks()
{
  CopyPhotonsBlock();
  CopyElectronsBlock();
  CopyMuonsBlock();
  CopyJetsBlock();
  
  CopyPhotonsNonIsoBlock();
  CopyElectronsMediumBlock();

  CopyEventBlock();
  CopyWeightBlock();
  CopyMetBlock();
  CopyOthersBlock();
}

void MiniClone::CopyPhotonsBlock()
{
  new_ph_n = ph_n;
  for (int i=0; i<ph_n; i++) {
    new_ph_pt->push_back((*ph_pt)[i]);
    new_ph_eta->push_back((*ph_eta)[i]);
    new_ph_etas2->push_back((*ph_etas2)[i]);
    new_ph_phi->push_back((*ph_phi)[i]);
    new_ph_etcone40->push_back((*ph_etcone40)[i]);
    new_ph_ptcone20->push_back((*ph_ptcone20)[i]);
    new_ph_iso->push_back((*ph_iso)[i]);
    new_ph_trackiso->push_back((*ph_trackiso)[i]);
    new_ph_conv->push_back((*ph_conv)[i]);
    new_ph_w->push_back((*ph_w)[i]);

    if (m_ismc) {
      new_ph_truth_pt->push_back((*ph_truth_pt)[i]);
      new_ph_truth_eta->push_back((*ph_truth_eta)[i]);
      new_ph_truth_phi->push_back((*ph_truth_phi)[i]);
      new_ph_truth_id->push_back((*ph_truth_id)[i]);
      new_ph_truth_type->push_back((*ph_truth_type)[i]);
      new_ph_truth_origin->push_back((*ph_truth_origin)[i]);
    }
    else {
      new_ph_truth_pt->push_back(0.);
      new_ph_truth_eta->push_back(0.);
      new_ph_truth_phi->push_back(0.);
      new_ph_truth_id->push_back(0);
      new_ph_truth_type->push_back(0);
      new_ph_truth_origin->push_back(0);
    }

  }

}

void MiniClone::CopyPhotonsNonIsoBlock()
{
  new_ph_noniso_n = ph_noniso_n;
  for (int i=0; i<ph_noniso_n; i++) {
    new_ph_noniso_pt->push_back((*ph_noniso_pt)[i]);
    new_ph_noniso_eta->push_back((*ph_noniso_eta)[i]);
    new_ph_noniso_etas2->push_back((*ph_noniso_etas2)[i]);
    new_ph_noniso_phi->push_back((*ph_noniso_phi)[i]);
    new_ph_noniso_etcone40->push_back((*ph_noniso_etcone40)[i]);
    new_ph_noniso_ptcone20->push_back((*ph_noniso_ptcone20)[i]);
    new_ph_noniso_iso->push_back((*ph_noniso_iso)[i]);
    new_ph_noniso_trackiso->push_back((*ph_noniso_trackiso)[i]);
    new_ph_noniso_conv->push_back((*ph_noniso_conv)[i]);
    new_ph_noniso_w->push_back((*ph_noniso_w)[i]);
  }

}

void MiniClone::CopyElectronsBlock()
{
  new_el_n = el_n;
  for (int i=0; i<el_n; i++) {
    new_el_pt->push_back((*el_pt)[i]);
    new_el_eta->push_back((*el_eta)[i]);
    new_el_etas2->push_back((*el_etas2)[i]);
    new_el_phi->push_back((*el_phi)[i]);
    new_el_ch->push_back((*el_ch)[i]);
    new_el_w->push_back((*el_w)[i]);
  }
}

void MiniClone::CopyElectronsMediumBlock()
{
  new_el_medium_n = el_medium_n;
  for (int i=0; i<el_medium_n; i++) {
    new_el_medium_pt->push_back((*el_medium_pt)[i]);
    new_el_medium_eta->push_back((*el_medium_eta)[i]);
    new_el_medium_etas2->push_back((*el_medium_etas2)[i]);
    new_el_medium_phi->push_back((*el_medium_phi)[i]);
    new_el_medium_ch->push_back((*el_medium_ch)[i]);
  }
}

void MiniClone::CopyMuonsBlock()
{
  new_mu_n = mu_n;
  for (int i=0; i<mu_n; i++) {
    new_mu_pt->push_back((*mu_pt)[i]);
    new_mu_eta->push_back((*mu_eta)[i]);
    new_mu_phi->push_back((*mu_phi)[i]);
    new_mu_ch->push_back((*mu_ch)[i]);
    new_mu_w->push_back((*mu_w)[i]);
  }

}

void MiniClone::CopyJetsBlock()
{
  new_jet_n = jet_n;
  new_bjet_n = bjet_n;

  for (int i=0; i<jet_n; i++) {
    new_jet_pt->push_back((*jet_pt)[i]);
    new_jet_eta->push_back((*jet_eta)[i]);
    new_jet_phi->push_back((*jet_phi)[i]);
    new_jet_e->push_back((*jet_e)[i]);
    new_jet_w->push_back((*jet_w)[i]);
    new_jet_isb->push_back((*jet_isb)[i]); 
  }
}

void MiniClone::CopyEventBlock()
{
  new_run = run;
  new_lb = lb;
  new_event = event;
  new_avgmu = avgmu;

  if (m_ismc) {
    new_fs = fs;
    new_mcveto = mcveto;
  }

  new_year = year;
  new_PRWHash = PRWHash;

  new_pass_g120 = pass_g120;
  new_pass_g140 = pass_g140;
  new_pass_g70_xe70 = pass_g70_xe70;
}

void MiniClone::CopyWeightBlock()
{
  new_weight_mc = weight_mc;
  new_weight_pu = weight_pu;
  new_weight_pu_down = weight_pu_down;
  new_weight_pu_up = weight_pu_up;
  new_weight_sf = weight_sf;

  if (m_efake_sample) {
    new_weight_feg = weight_feg;
    new_weight_feg_dn = weight_feg_dn;
    new_weight_feg_up = weight_feg_up;
  }
  if (m_efake_sample) {
    new_weight_fjg = weight_fjg;
    new_weight_fjg_dn = weight_fjg_dn;
    new_weight_fjg_up = weight_fjg_up;
  }
}

void MiniClone::CopyMetBlock()
{
  new_met_et  = met_et;
  new_met_phi = met_phi;

  if (!m_is_old_version) {
    new_met_sumet = met_sumet;
    new_met_sig   = met_sig;

    new_met_track_et  = met_track_et;
    new_met_track_phi = met_track_phi;

    new_met_soft_et  = met_soft_et;
    new_met_soft_phi = met_soft_phi;
    new_met_ele_et  = met_ele_et;
    new_met_ele_phi = met_ele_phi;
    new_met_gam_et  = met_gam_et;
    new_met_gam_phi = met_gam_phi;
    new_met_muon_et  = met_muon_et;
    new_met_muon_phi = met_muon_phi;
    new_met_jet_et  = met_jet_et;
    new_met_jet_phi = met_jet_phi;
  }
}

void MiniClone::CopyOthersBlock()
{
  new_ht0 = ht0;
  new_ht = ht;
  new_meff = meff;
  new_rt1	= rt1;
  new_rt2	= rt2;
  new_rt3	= rt3;
  new_rt4	= rt4;

  new_mt_gam = mt_gam;

  new_dphi_gamjet = dphi_gamjet;
  new_dphi_jetmet = dphi_jetmet;
  new_dphi_jet1met = dphi_jet1met;
  new_dphi_gammet = dphi_gammet;

  new_dphi_met_trackmet = dphi_met_trackmet;
}


void MiniClone::AddNewBranch(TString name, Float_t *address)
{
  clone_tree->Branch(name, address);
}

void MiniClone::AddNewBranch(TString name, Int_t *address)
{
  clone_tree->Branch(name, address);
}

void MiniClone::Fill()
{
  clone_tree->Fill();
}

void MiniClone::Save()
{
  // just for consistency. They are wrong
  TH1D *events = (TH1D*)orig_file->Get("events");
  TH1D *cutflow = (TH1D*)orig_file->Get("cutflow");
  TH1D *cutflow_w = (TH1D*)orig_file->Get("cutflow_w");
  TH1D *susy_sumw = (TH1D*)orig_file->Get("susy_sumw");
  
  clone_file->cd();
  if (events) events->Write(); 
  if (cutflow) cutflow->Write(); 
  if (cutflow_w) cutflow_w->Write(); 
  if (susy_sumw) susy_sumw->Write(); 

  clone_tree->Write();
  clone_file->Close();
}

// void MiniClone::Loop()
// {
//    if (fChain == 0) return;

//    Long64_t nentries = fChain->GetEntriesFast();

//    Long64_t nbytes = 0, nb = 0;
//    for (Long64_t jentry=0; jentry<nentries;jentry++) {
//       Long64_t ientry = LoadTree(jentry);
//       if (ientry < 0) break;
//       nb = fChain->GetEntry(jentry);   nbytes += nb;
//       // if (Cut(ientry) < 0) continue;
//    }
// }
