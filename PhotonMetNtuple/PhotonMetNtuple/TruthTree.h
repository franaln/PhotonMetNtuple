#ifndef TruthTree_h
#define TruthTree_h

#include <iomanip>
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "AsgTools/AsgMetadataTool.h"

#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>

class TruthTree : public asg::AsgMetadataTool {
  
 public :
  TruthTree(const std::string& name);
  ~TruthTree();
  
  StatusCode initialize();
  bool process();
  inline StatusCode Fill() { tree->Fill(); return StatusCode::SUCCESS; };
  

  TTree *tree;

  void CreateTree(TString); 

  void CreateBranches();
  void Clear();

  // multiplicities
  Int_t ph_n;
  Int_t el_n;
  Int_t mu_n;
  Int_t jet_n;
  Int_t bjet_n;
  
  // photon
  std::vector<float> *ph_pt;
  std::vector<float> *ph_iso;
  std::vector<float> *ph_eta;
  std::vector<float> *ph_phi;
  std::vector<int>   *ph_type;

  // electrons
  std::vector<float> *el_pt;
  std::vector<float> *el_eta;
  std::vector<float> *el_phi;

  // muons
  std::vector<float> *mu_pt;
  std::vector<float> *mu_eta;
  std::vector<float> *mu_phi;

  // jets
  std::vector<float> *jet_pt;
  std::vector<float> *jet_eta;
  std::vector<float> *jet_phi;

  // met
  Float_t met_et;
  Float_t met_phi;
  Float_t met_sumet;

  Float_t met_truth_et;
  Float_t met_truth_phi;

  // other
  Float_t meff;
  Float_t ht;
  Float_t rt2;
  Float_t rt4;

  // dphi
  Float_t dphi_gamjet;
  Float_t dphi_jetmet;
  Float_t dphi_gammet;

  std::vector<float> *weight_pdf1;
  std::vector<float> *weight_pdf2;
  std::vector<float> *weight_pdf3;

 protected:
  TDirectory *m_outfile;    
  bool m_pdfrw;

};
#endif
