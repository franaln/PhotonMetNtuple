#ifndef BaselineTree_h
#define BaselineTree_h

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

struct AnalysisBaselineCollections {

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
  
};

class BaselineTree : public asg::AsgMetadataTool {

 public:
  BaselineTree(const std::string& name);
  ~BaselineTree();

  StatusCode initialize();
  void clear();
  bool process(AnalysisBaselineCollections collections);
  StatusCode FillTree();

  void SetEventNumber(int en) { event_number = en; };
  void SetRunNumber(int en) { run_number = en; };

  TTree* tree;
  
protected:
  TDirectory *m_outfile;    
  std::vector<ST::SystInfo> m_sysList;
  Bool_t m_ismc;

  int ph_n;
  int el_n;
  int mu_n;
  int jet_n;
  int bjet_n;

  std::vector<float> *ph_pt; 
  std::vector<float> *ph_eta;
  std::vector<float> *ph_etas2;
  std::vector<float> *ph_phi;
  std::vector<float> *ph_iso;
  std::vector<int> *ph_passOR;
  std::vector<int> *ph_signal;
  std::vector<int> *ph_isol;

  std::vector<float> *jet_pt; 
  std::vector<float> *jet_eta;
  std::vector<float> *jet_phi;
  std::vector<float> *jet_e;
  std::vector<float> *jet_w;
  std::vector<bool>  *jet_isb;
  std::vector<int>  *jet_passOR;
  std::vector<int>  *jet_signal;
  
  std::vector<float> *el_pt; 
  std::vector<float> *el_eta;
  std::vector<float> *el_etas2;
  std::vector<float> *el_phi;
  std::vector<int>   *el_ch;
  std::vector<int> *el_passOR;
  std::vector<int> *el_signal;
  std::vector<int> *el_isol;
  std::vector<float> *mu_pt; 
  std::vector<float> *mu_eta;
  std::vector<float> *mu_phi;
  std::vector<int>   *mu_ch;
  std::vector<int> *mu_passOR;
  std::vector<int> *mu_signal;
  std::vector<int> *mu_isol;
  
  int event_number;
  int run_number;

  float met_et;
  float met_phi;

};

#endif
