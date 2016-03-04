#ifndef OutTree_h
#define OutTree_h

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

  Int_t event_number;

  Double_t weight_pu;
  Float_t weight_mc;
  Float_t weight_btag;
  Float_t cross_section;
  Float_t primvx_z;
  Int_t avg_mu;

  // Containers
  xAOD::PhotonContainer* photons;
  xAOD::ElectronContainer* electrons;
  xAOD::MuonContainer* muons;
  xAOD::JetContainer* jets;
  xAOD::MissingETContainer* met;
  xAOD::MissingETContainer* met_tst;

  // Aux containers too
  xAOD::ShallowAuxContainer* photons_aux;
  xAOD::ShallowAuxContainer* electrons_aux;
  xAOD::ShallowAuxContainer* muons_aux;
  xAOD::ShallowAuxContainer* jets_aux;
  xAOD::MissingETAuxContainer* met_aux;
  xAOD::MissingETAuxContainer* met_tst_aux;

  //Truth containers
  xAOD::TruthParticleContainer* truth_p;
  xAOD::JetContainer* truth_jets; 
  TLorentzVector truth_met;
};

class OutTree : public asg::AsgMetadataTool {

 public:
  OutTree(const std::string& name);
  ~OutTree();
  StatusCode initialize();
  bool process(AnalysisCollections collections, std::string sysname);
  StatusCode FillTree();

  TString BookName(TString branch, TString sys_name);
  
  TTree* tree;
  
  std::map<const std::string, int> ph_n_map;
  std::map<const std::string, std::vector<float>*> ph_pt_map;
  std::map<const std::string, std::vector<float>*> ph_eta_map;
  std::map<const std::string, std::vector<float>*> ph_phi_map;
  std::map<const std::string, std::vector<float>*> ph_w_map;

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
  std::map<const std::string, std::vector<float>*> el_phi_map;
  std::map<const std::string, std::vector<int>*>   el_ch_map;
  std::map<const std::string, std::vector<float>*> el_w_map;

  std::map<const std::string, float> met_phi_map;
  std::map<const std::string, float> met_et_map;

  std::map<const std::string, float> ht_map;
  std::map<const std::string, float> meff_map;
  std::map<const std::string, float> rt2_map;
  std::map<const std::string, float> rt4_map;

  std::map<const std::string, float> dphi_jetmet_map;
  std::map<const std::string, float> dphi_gamjet_map;
  std::map<const std::string, float> dphi_gammet_map;

  std::map<const std::string, float> mgj_map;
  std::map<const std::string, float> mgjj_map;
  std::map<const std::string, float> mgjjj_map;

protected:
  TDirectory *m_outfile;    
  std::vector<ST::SystInfo> m_sysList;

  std::vector<float> *ph_loose_pt; 
  std::vector<float> *ph_loose_eta;
  std::vector<float> *ph_loose_phi;
  std::vector<unsigned int> *ph_loose_isem;
  std::vector<float> *ph_loose_iso20;
  std::vector<float> *ph_loose_iso40;

  std::vector<float> *ph_pt; 
  std::vector<float> *ph_eta;
  std::vector<float> *ph_phi;
  std::vector<float> *ph_w;

  std::vector<float> *jet_pt; 
  std::vector<float> *jet_eta;
  std::vector<float> *jet_phi;
  std::vector<float> *jet_e;
  std::vector<float> *jet_w;
  std::vector<bool>  *jet_isb;
  
  std::vector<float> *el_pt; 
  std::vector<float> *el_eta;
  std::vector<float> *el_phi;
  std::vector<int>   *el_ch;
  std::vector<float> *el_w;
  
  std::vector<float> *mu_pt; 
  std::vector<float> *mu_eta;
  std::vector<float> *mu_phi;
  std::vector<int>   *mu_ch;
  std::vector<float> *mu_w;
  
  int ph_loose_n;

  int event_number;
  float weight_mc;
  float weight_pu;
  float weight_sf;

  int avg_mu;
};

#endif
