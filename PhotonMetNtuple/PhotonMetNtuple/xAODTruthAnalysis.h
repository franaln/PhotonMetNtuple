#ifndef MyAnalysis_xAODTruthAnalysis_H
#define MyAnalysis_xAODTruthAnalysis_H

#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#include <TH1.h>
#include <TTree.h>

#include <PhotonMetNtuple/TruthTree.h>


namespace LHAPDF {
  class PDF;
}

class xAODTruthAnalysis : public EL::Algorithm
{
  
private:

public:
  
  bool do_pdfrw;
  bool do_lhe3;
  bool is_truth3;

  TruthTree *ntuple; //!

  TH1D *h_events; //!

  std::map<std::string, int> map_lhe3; //!
  TH1D *h_lhe3_sumw; //!

  xAOD::TEvent *m_event;  //!

  // this is a standard constructor
  xAODTruthAnalysis();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();


  //PDF reweighting
  std::string pdf1; //!
  std::string pdf2; //!
  std::string pdf3; //!
  std::vector<LHAPDF::PDF*> m_pdfs_1; //!
  std::vector<LHAPDF::PDF*> m_pdfs_2; //!
  std::vector<LHAPDF::PDF*> m_pdfs_3; //!
  // std::vector<double> *weights_pdf_1; //!
  // std::vector<double> *weights_pdf_2; //!
  // std::vector<double> *weights_pdf_3; //!
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(xAODTruthAnalysis, 1);
};

#endif
