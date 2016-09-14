#ifndef mini_h
#define mini_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class MiniClone {

public:

  TTree *orig_tree;   //!
  TTree *clone_tree;   //

  TFile *orig_file; //!
  TFile *clone_file; //!

  // Original variables
  Int_t           run;
  Int_t           lb;
  ULong64_t       event;
  Float_t         avgmu;
  UInt_t          fs;
  UInt_t          pass_tst_cleaning;
  UInt_t          year;
  UInt_t          mcveto;
  UInt_t          pass_g120;
  UInt_t          pass_g140;
  Float_t         weight_mc;
  Float_t         weight_pu;
  Float_t         weight_pu_down;
  Float_t         weight_pu_up;
  Float_t         weight_sf;
  Float_t         weight_btag;

  Float_t weight_feg;
  Float_t weight_feg_dn;
  Float_t weight_feg_up;

  Float_t weight_fjg;
  Float_t weight_fjg_dn;
  Float_t weight_fjg_up;

   ULong64_t       PRWHash;
   Int_t           ph_n;
   vector<float>   *ph_pt;
   vector<float>   *ph_eta;
   vector<float>   *ph_etas2;
   vector<float>   *ph_phi;
   vector<float>   *ph_iso40;
   vector<float>   *ph_iso;
   vector<float>   *ph_w;
   vector<float>   *ph_truth_pt;
   vector<float>   *ph_truth_eta;
   vector<float>   *ph_truth_phi;
   vector<int>     *ph_truth_id;
   vector<int>     *ph_truth_type;
   vector<int>     *ph_truth_origin;
   Int_t           ph_noniso_n;
   vector<float>   *ph_noniso_pt;
   vector<float>   *ph_noniso_eta;
   vector<float>   *ph_noniso_phi;
   vector<float>   *ph_noniso_iso;
   vector<float>   *ph_noniso_w;
   Int_t           jet_n;
   Int_t           bjet_n;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<float>   *jet_pt;
   vector<float>   *jet_e;
   vector<bool>    *jet_isb;
   vector<float>   *jet_w;
   Int_t           el_n;
   vector<float>   *el_eta;
   vector<float>   *el_etas2;
   vector<float>   *el_phi;
   vector<float>   *el_pt;
   vector<int>     *el_ch;
   vector<float>   *el_w;
   Int_t           el_medium_n;
   vector<float>   *el_medium_eta;
   vector<float>   *el_medium_etas2;
   vector<float>   *el_medium_phi;
   vector<float>   *el_medium_pt;
   vector<int>     *el_medium_ch;
   Int_t           mu_n;
   vector<float>   *mu_eta;
   vector<float>   *mu_phi;
   vector<float>   *mu_pt;
   vector<int>     *mu_ch;
   vector<float>   *mu_w;
   Float_t         met_et;
   Float_t         met_phi;
   Float_t         met_sumet;
   Float_t         met_sig;
   Float_t         tst_et;
   Float_t         tst_phi;
   Float_t         ht;
   Float_t         meff;
   Float_t         rt1;
   Float_t         rt2;
   Float_t         rt3;
   Float_t         rt4;
   Float_t         dphi_gamjet;
   Float_t         dphi_jetmet;
   Float_t         dphi_gammet;

   // Clone variables
   Int_t           new_run;
   Int_t           new_lb;
   ULong64_t       new_event;
   Float_t         new_avgmu;
   UInt_t          new_fs;
   UInt_t          new_pass_tst_cleaning;
   UInt_t          new_year;
   UInt_t          new_mcveto;
   UInt_t          new_pass_g120;
   UInt_t          new_pass_g140;
   Float_t         new_weight_mc;
   Float_t         new_weight_pu;
   Float_t         new_weight_pu_down;
   Float_t         new_weight_pu_up;
   Float_t         new_weight_sf;
   Float_t         new_weight_btag;

  Float_t new_weight_feg;
  Float_t new_weight_feg_dn;
  Float_t new_weight_feg_up;

  Float_t new_weight_fjg;
  Float_t new_weight_fjg_dn;
  Float_t new_weight_fjg_up;

   ULong64_t       new_PRWHash;
   Int_t           new_ph_n;
   vector<float>   *new_ph_pt;
   vector<float>   *new_ph_eta;
   vector<float>   *new_ph_etas2;
   vector<float>   *new_ph_phi;
   vector<float>   *new_ph_iso40;
   vector<float>   *new_ph_iso;
   vector<float>   *new_ph_w;
   vector<float>   *new_ph_truth_pt;
   vector<float>   *new_ph_truth_eta;
   vector<float>   *new_ph_truth_phi;
   vector<int>     *new_ph_truth_id;
   vector<int>     *new_ph_truth_type;
   vector<int>     *new_ph_truth_origin;
   Int_t           new_ph_noniso_n;
   vector<float>   *new_ph_noniso_pt;
   vector<float>   *new_ph_noniso_eta;
   vector<float>   *new_ph_noniso_phi;
   vector<float>   *new_ph_noniso_iso;
   vector<float>   *new_ph_noniso_w;
   Int_t           new_jet_n;
   Int_t           new_bjet_n;
   vector<float>   *new_jet_eta;
   vector<float>   *new_jet_phi;
   vector<float>   *new_jet_pt;
   vector<float>   *new_jet_e;
   vector<bool>    *new_jet_isb;
   vector<float>   *new_jet_w;
   Int_t           new_el_n;
   vector<float>   *new_el_eta;
   vector<float>   *new_el_etas2;
   vector<float>   *new_el_phi;
   vector<float>   *new_el_pt;
   vector<int>     *new_el_ch;
   vector<float>   *new_el_w;
   Int_t           new_el_medium_n;
   vector<float>   *new_el_medium_eta;
   vector<float>   *new_el_medium_etas2;
   vector<float>   *new_el_medium_phi;
   vector<float>   *new_el_medium_pt;
   vector<int>     *new_el_medium_ch;
   Int_t           new_mu_n;
   vector<float>   *new_mu_eta;
   vector<float>   *new_mu_phi;
   vector<float>   *new_mu_pt;
   vector<int>     *new_mu_ch;
   vector<float>   *new_mu_w;
   Float_t         new_met_et;
   Float_t         new_met_phi;
   Float_t         new_met_sumet;
   Float_t         new_met_sig;
   Float_t         new_tst_et;
   Float_t         new_tst_phi;
   Float_t         new_ht;
   Float_t         new_meff;
   Float_t         new_rt1;
   Float_t         new_rt2;
   Float_t         new_rt3;
   Float_t         new_rt4;
   Float_t         new_dphi_gamjet;
   Float_t         new_dphi_jetmet;
   Float_t         new_dphi_gammet;

   // List of original branches
   TBranch        *b_run;   //!
   TBranch        *b_lb;   //!
   TBranch        *b_event;   //!
   TBranch        *b_avgmu;   //!
   TBranch        *b_final_state;   //!
   TBranch        *b_pass_tst_cleaning;   //!
   TBranch        *b_year;   //!
   TBranch        *b_mcveto;   //!
   TBranch        *b_pass_g120;   //!
   TBranch        *b_pass_g140;   //!
   TBranch        *b_weight_mc;   //!
   TBranch        *b_weight_pu;   //!
   TBranch        *b_weight_pu_down;   //!
   TBranch        *b_weight_pu_up;   //!
   TBranch        *b_weight_sf;   //!
   TBranch        *b_weight_btag;   //!
   TBranch        *b_weight_feg;   //!
   TBranch        *b_weight_feg_dn;   //!
   TBranch        *b_weight_feg_up;   //!
   TBranch        *b_weight_fjg;   //!
   TBranch        *b_weight_fjg_dn;   //!
   TBranch        *b_weight_fjg_up;   //!
   TBranch        *b_PRWHash;   //!
   TBranch        *b_ph_n;   //!
   TBranch        *b_ph_pt;   //!
   TBranch        *b_ph_eta;   //!
   TBranch        *b_ph_etas2;   //!
   TBranch        *b_ph_phi;   //!
   TBranch        *b_ph_iso40;   //!
   TBranch        *b_ph_iso;   //!
   TBranch        *b_ph_w;   //!
   TBranch        *b_ph_truth_pt;   //!
   TBranch        *b_ph_truth_eta;   //!
   TBranch        *b_ph_truth_phi;   //!
   TBranch        *b_ph_truth_id;   //!
   TBranch        *b_ph_truth_type;   //!
   TBranch        *b_ph_truth_origin;   //!
   TBranch        *b_ph_noniso_n;   //!
   TBranch        *b_ph_noniso_pt;   //!
   TBranch        *b_ph_noniso_eta;   //!
   TBranch        *b_ph_noniso_phi;   //!
   TBranch        *b_ph_noniso_iso;   //!
   TBranch        *b_ph_noniso_w;   //!
   TBranch        *b_jet_n;   //!
   TBranch        *b_bjet_n;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_e;   //!
   TBranch        *b_jet_isb;   //!
   TBranch        *b_jet_w;   //!
   TBranch        *b_el_n;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_etas2;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_ch;   //!
   TBranch        *b_el_w;   //!
   TBranch        *b_el_medium_n;   //!
   TBranch        *b_el_medium_eta;   //!
   TBranch        *b_el_medium_etas2;   //!
   TBranch        *b_el_medium_phi;   //!
   TBranch        *b_el_medium_pt;   //!
   TBranch        *b_el_medium_ch;   //!
   TBranch        *b_mu_n;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_pt;   //!
   TBranch        *b_mu_ch;   //!
   TBranch        *b_mu_w;   //!
   TBranch        *b_met_et;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_met_sig;   //!
   TBranch        *b_tst_et;   //!
   TBranch        *b_tst_phi;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_meff;   //!
   TBranch        *b_rt1;   //!
   TBranch        *b_rt2;   //!
   TBranch        *b_rt3;   //!
   TBranch        *b_rt4;   //!
   TBranch        *b_dphi_gamjet;   //!
   TBranch        *b_dphi_jetmet;   //!
   TBranch        *b_dphi_gammet;   //!


  bool m_ismc;
  bool m_efake_sample;
  bool m_jfake_sample;

  MiniClone(TString, TString);
  virtual ~MiniClone();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual void     InitOriginalTree();
  virtual void     CreateCloneTree();
  
  void AddNewBranch(TString, Float_t*);
  void AddNewBranch(TString, Int_t*);
  void Fill();
  void Save();
  Long64_t GetEntries() { return orig_tree->GetEntriesFast(); };

  void SetMC() { m_ismc = true; };
  void SetEfakeSample() { m_efake_sample = true; };
  void SetJfakeSample() { m_jfake_sample = true; };

  void Clear();
  void CopyPhotonsBlock();
  void CopyPhotonsNonIsoBlock();
  void CopyElectronsBlock();
  void CopyElectronsMediumBlock();
  void CopyMuonsBlock();
  void CopyJetsBlock();
  void CopyEventBlock();

};

#endif
