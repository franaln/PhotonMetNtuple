#ifndef xAODCountEwkProcesses_H
#define xAODCountEwkProcesses_H

#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#include "xAODTruth/TruthParticleContainer.h"

#include <TH1.h>

class xAODCountEwkProcesses : public EL::Algorithm
{

private:

  unsigned int GetFinalState(const int SUSY_Spart1_pdgId, const int SUSY_Spart2_pdgId);

  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:

  TH1D *h_hp; //!

  xAOD::TEvent *m_event;  //!

  // this is a standard constructor
  xAODCountEwkProcesses ();

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
  ClassDef(xAODCountEwkProcesses, 1);
};

#endif
