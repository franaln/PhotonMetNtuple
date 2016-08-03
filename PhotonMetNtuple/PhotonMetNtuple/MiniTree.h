#ifndef MiniTree_h
#define MiniTree_h

#include <vector>
#include <map>
#include <string>
#include "TTree.h"
#include "AsgTools/AsgMetadataTool.h"
#include "PATInterfaces/SystematicSet.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "SUSYTools/ISUSYObjDef_xAODTool.h"

struct AnalysisCollections {

  // Containers
  xAOD::PhotonContainer* photons;
  xAOD::ElectronContainer* electrons;
  xAOD::MuonContainer* muons;
  xAOD::JetContainer* jets;
  xAOD::MissingETContainer* met;

  // Aux containers too
  xAOD::ShallowAuxContainer* photons_aux;
  xAOD::ShallowAuxContainer* electrons_aux;
  xAOD::ShallowAuxContainer* muons_aux;
  xAOD::ShallowAuxContainer* jets_aux;
  xAOD::MissingETAuxContainer* met_aux;

  //Truth containers
  xAOD::TruthParticleContainer* truth_p;
  xAOD::JetContainer* truth_jets; 
  TLorentzVector truth_met;
};

class MiniTree : public asg::AsgMetadataTool {

 public:
  MiniTree(const std::string& name);
  ~MiniTree();

  StatusCode initialize();
  void clear();
  bool process(AnalysisCollections collections, std::string sysname);
  StatusCode FillTree();
  TString BookName(TString branch, TString sys_name);

  bool PassEtaCut(const xAOD::IParticle *part, Float_t maxeta=99.);
  bool PassEtaCut(const xAOD::Photon *part, Float_t maxeta=99.);
  bool PassEtaCut(const xAOD::Electron *part, Float_t maxeta=99.);

  void SetRunNumber(int arg)    { run_number = arg; };
  void SetLumiBlock(int arg)    { lumi_block = arg; };
  void SetEventNumber(ULong64_t arg) { event_number = arg; };
  void SetAvgMu(float arg)      { avg_mu = arg; };

  void SetMCFinalState(unsigned int fs) { final_state = fs; };
  
  void SetWeightMC(float w) { weight_mc = w; };
  void SetWeightBtag(std::string sysname, float w) { weight_btag_map[sysname] = w; };
  void SetPRWHash(ULong64_t hash) { PRWHash = hash; };

  void set_weight_pu(float w) { weight_pu = w; };
  void set_weight_pu_down(float w) { weight_pu_down = w; };
  void set_weight_pu_up  (float w) { weight_pu_up = w; };

  void SetPassTSTCleaning(int w) { pass_tst_cleaning = w; };
  void SetYear(unsigned int y) { year = y; };

  void set_pass_g120(unsigned int fs) { pass_g120 = fs; };
  void set_pass_g140(unsigned int fs) { pass_g140 = fs; };
  void set_mcveto(unsigned int fs) { mcveto = fs; };

  TTree* tree;
  
  std::map<const std::string, int> ph_n_map;
  std::map<const std::string, std::vector<float>*> ph_pt_map;
  std::map<const std::string, std::vector<float>*> ph_eta_map;
  std::map<const std::string, std::vector<float>*> ph_etas2_map;
  std::map<const std::string, std::vector<float>*> ph_phi_map;
  std::map<const std::string, std::vector<float>*> ph_iso40_map;
  std::map<const std::string, std::vector<float>*> ph_iso_map;
  std::map<const std::string, std::vector<float>*> ph_w_map;

  std::map<const std::string, int> ph_noniso_n_map;
  std::map<const std::string, std::vector<float>*> ph_noniso_pt_map;
  std::map<const std::string, std::vector<float>*> ph_noniso_eta_map;
  std::map<const std::string, std::vector<float>*> ph_noniso_phi_map;
  std::map<const std::string, std::vector<float>*> ph_noniso_iso_map;
  std::map<const std::string, std::vector<float>*> ph_noniso_w_map;

  std::map<const std::string, int> jet_n_map;
  std::map<const std::string, int> bjet_n_map;
  std::map<const std::string, std::vector<float>*> jet_pt_map;
  std::map<const std::string, std::vector<float>*> jet_eta_map;
  std::map<const std::string, std::vector<float>*> jet_phi_map;
  std::map<const std::string, std::vector<float>*> jet_e_map;
  std::map<const std::string, std::vector<bool>*>  jet_isb_map;
  std::map<const std::string, std::vector<float>*> jet_w_map;

  std::map<const std::string, int> mu_n_map;
  std::map<const std::string, std::vector<float>*> mu_pt_map;
  std::map<const std::string, std::vector<float>*> mu_eta_map;
  std::map<const std::string, std::vector<float>*> mu_phi_map;
  std::map<const std::string, std::vector<int>*>   mu_ch_map;
  std::map<const std::string, std::vector<float>*> mu_w_map;

  std::map<const std::string, int> el_n_map;
  std::map<const std::string, std::vector<float>*> el_pt_map;
  std::map<const std::string, std::vector<float>*> el_eta_map;
  std::map<const std::string, std::vector<float>*> el_etas2_map;
  std::map<const std::string, std::vector<float>*> el_phi_map;
  std::map<const std::string, std::vector<int>*>   el_ch_map;
  std::map<const std::string, std::vector<float>*> el_w_map;

  std::map<const std::string, int> el_medium_n_map;
  std::map<const std::string, std::vector<float>*> el_medium_pt_map;
  std::map<const std::string, std::vector<float>*> el_medium_eta_map;
  std::map<const std::string, std::vector<float>*> el_medium_etas2_map;
  std::map<const std::string, std::vector<float>*> el_medium_phi_map;
  std::map<const std::string, std::vector<int>*>   el_medium_ch_map;

  std::map<const std::string, float> met_phi_map;
  std::map<const std::string, float> met_et_map;

  std::map<const std::string, float> tst_phi_map;
  std::map<const std::string, float> tst_et_map;

  std::map<const std::string, float> ht_map;
  std::map<const std::string, float> meff_map;

  std::map<const std::string, float> rt1_map;
  std::map<const std::string, float> rt2_map;
  std::map<const std::string, float> rt3_map;
  std::map<const std::string, float> rt4_map;

  std::map<const std::string, float> dphi_jetmet_map;
  std::map<const std::string, float> dphi_gamjet_map;
  std::map<const std::string, float> dphi_gammet_map;

  std::map<const std::string, float> weight_sf_map;
  std::map<const std::string, float> weight_btag_map;
    

protected:
  TDirectory *m_outfile;    
  std::vector<ST::SystInfo> m_sysList;
  Bool_t m_ismc;

  std::vector<float> *ph_pt; 
  std::vector<float> *ph_eta;
  std::vector<float> *ph_etas2;
  std::vector<float> *ph_phi;
  std::vector<float> *ph_iso40;
  std::vector<float> *ph_iso;
  std::vector<float> *ph_w;

  std::vector<float> *ph_truth_pt; 
  std::vector<float> *ph_truth_eta;
  std::vector<float> *ph_truth_phi;
  std::vector<int>   *ph_truth_origin;
  std::vector<int>   *ph_truth_type;
  std::vector<int>   *ph_truth_id;

  std::vector<float> *ph_noniso_pt; 
  std::vector<float> *ph_noniso_eta;
  std::vector<float> *ph_noniso_phi;
  std::vector<float> *ph_noniso_iso;
  std::vector<float> *ph_noniso_w;

  std::vector<float> *jet_pt; 
  std::vector<float> *jet_eta;
  std::vector<float> *jet_phi;
  std::vector<float> *jet_e;
  std::vector<float> *jet_w;
  std::vector<bool>  *jet_isb;
  
  std::vector<float> *el_pt; 
  std::vector<float> *el_eta;
  std::vector<float> *el_etas2;
  std::vector<float> *el_phi;
  std::vector<int>   *el_ch;
  std::vector<float> *el_w;

  std::vector<float> *el_medium_pt; 
  std::vector<float> *el_medium_eta;
  std::vector<float> *el_medium_etas2;
  std::vector<float> *el_medium_phi;
  std::vector<int>   *el_medium_ch;
    
  std::vector<float> *mu_pt; 
  std::vector<float> *mu_eta;
  std::vector<float> *mu_phi;
  std::vector<int>   *mu_ch;
  std::vector<float> *mu_w;
  
  int run_number;
  int lumi_block;
  ULong64_t event_number;
  float avg_mu;

  unsigned int final_state;
  unsigned int pass_tst_cleaning;
  unsigned int year;

  unsigned int pass_g120;
  unsigned int pass_g140;

  unsigned int mcveto;

  float weight_mc;
  float weight_sf;
  float weight_btag;
  float weight_pu;
  float weight_pu_down;
  float weight_pu_up;

  ULong64_t PRWHash;

};

#endif
