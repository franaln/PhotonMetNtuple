/* single photon analysis
   MiniWriter.h */

#ifndef MiniWriter_h
#define MiniWriter_h


#include <iomanip>
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>

#include <TROOT.h>
#include <TChain.h>
#include <TString.h>
#include <TFile.h>
#include <TH1.h>

using namespace std;

class MiniWriter {
  
 public :
  MiniWriter();
  ~MiniWriter();
  
  void Create(TString, TString);
  void CreateFile(TString);
  void CreateTree(TString, Long64_t=-1);

  inline void Fill() { tree->Fill(); };
  void Save();
  inline void Close() { file->Close(); };
  void SaveHist(TH1*);
  void RemoveObjectIfExists(TString);

  void CreateBranches();
  void Clear();

  // multiplicities
  Int_t ph_n;
  Int_t el_n;
  Int_t mu_n;
  Int_t jet_n;
  Int_t bjet_n;
  
  // photon
  vector<float> *ph_pt;
  vector<float> *ph_iso;
  vector<float> *ph_eta;
  vector<float> *ph_phi;
  vector<int> *ph_type;

  // electrons
  vector<float> *el_pt;
  vector<float> *el_eta;
  vector<float> *el_phi;

  // muons
  vector<float> *mu_pt;
  vector<float> *mu_eta;
  vector<float> *mu_phi;

  // jets
  vector<float> *jet_pt;
  vector<float> *jet_eta;
  vector<float> *jet_phi;

  // met
  Float_t met_et;
  Float_t met_phi;
  Float_t met_sumet;

  Float_t met_truth_et;
  Float_t met_truth_phi;

  // other
  Float_t ht;
  Float_t rt2;
  Float_t rt4;

  // dphi
  Float_t dphi_gamjet;
  Float_t dphi_jetmet;
  Float_t dphi_gammet;

  // Long64_t event;
  // Int_t smeared;

  // // weights
  // Float_t weight;
  // Float_t weight_prwdn;
  // Float_t weight_prwup;
  // Float_t weight_feg;
  // vector<float> *btag_weight;

 protected:
  TFile *file;
  TTree *tree;

};
#endif
