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

class xAODTruthAnalysis : public EL::Algorithm
{
  
  #ifndef __CINT__
    TruthTree *mem_leaker; //!
  #endif // not __CINT__

private:

public:
  TruthTree *ntuple; //!

  TH1D *h_events; //!

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
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(xAODTruthAnalysis, 1);
};

#endif
