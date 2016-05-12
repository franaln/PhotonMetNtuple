#ifndef MyAnalysis_xAODAnalysis_H
#define MyAnalysis_xAODAnalysis_H

#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"

#include <TH1.h>
#include <TTree.h>
#include "PATInterfaces/SystematicCode.h"
#include "PATInterfaces/SystematicSet.h"
#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicVariation.h"
#include <SUSYTools/SUSYCrossSection.h>
//#include <PhotonMetNtuple/MiniTree.h>
#include <PhotonMetNtuple/MiniTree2.h>
#include "SUSYTools/ISUSYObjDef_xAODTool.h"

#include <PhotonMetNtuple/MCFilter.h>

// GRL
class GoodRunsListSelectionTool;
namespace CP{
  class PileupReweightingTool;
}
using namespace CP;
class JetCleaningTool;

namespace ST{
  class SUSYObjDef_xAOD;
}
using namespace ST;

class JetVertexTagger;

class xAODAnalysis : public EL::Algorithm
{
  
  SUSY::CrossSectionDB *my_XsecDB;  //!

#ifndef __CINT__
  GoodRunsListSelectionTool *m_grl; //!
  PileupReweightingTool *m_pileupReweightingTool; //!
  SUSYObjDef_xAOD *objTool; //!
  //  MiniTree *mem_leaker; //!
#endif // not __CINT__

private:

  int event_number; //!
  int run_number; //!

  double m_xsec; //!
  bool isDerived; //!

  double weight_mc; //!
  double weight_pu; //!

  Int_t avg_mu; //!
  Int_t tbcid; //!
  Int_t tLB; //!
  Double_t tEventWeight; //!
  Int_t tChannelNumber; //!
  Int_t tDetError; //!

  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  bool is_data;
  bool is_atlfast;
  bool is_susy;
  bool is_susy_ewk;
  bool do_syst;

  //  std::vector<CP::SystematicSet> sysList; //!

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  MiniTree2 *outtree; //!

  MCFilter *mc_filter; //!

  TH1D *h_events; //!
  TH1D *h_cutflow; //!
  TH1D *h_events_subproceses; //!
  TH1D *h_sumw_subproceses; //!

  xAOD::TEvent *m_event;  //!

  float  m_initialSumOfWeights; //!
  std::vector<ST::SystInfo> systInfoList; //!

  // this is a standard constructor
  xAODAnalysis ();

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
  ClassDef(xAODAnalysis, 1);
};

#endif
