#include <TSystem.h>

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <PhotonMetNtuple/xAODBaselineAnalysis.h>
#include "EventLoop/OutputStream.h"
#include <TTreeFormula.h>

// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODBTaggingEfficiency/BTaggingEfficiencyTool.h"
#include "xAODBTagging/BTagging.h"
#include "FourMomUtils/xAODP4Helpers.h"

// Tools
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "PileupReweighting/PileupReweightingTool.h"
#include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "SUSYTools/SUSYCrossSection.h"

// Amg include
#include "EventPrimitives/EventPrimitivesHelpers.h"
#include "xAODEgamma/EgammaxAODHelpers.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

static const char *APP_NAME = "";

// this is needed to distribute the algorithm to the workers
ClassImp(xAODBaselineAnalysis)


xAODBaselineAnalysis::xAODBaselineAnalysis() : 
  m_grl(0), susytools(0),
  config_file(""),
  is_data(false), 
  is_atlfast(false), 
  is_susy(false),
  is_susy_ewk(false),
  do_syst(false)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}

EL::StatusCode xAODBaselineAnalysis::setupJob(EL::Job &job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  job.useXAOD();
  
  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init("xAODBaselineAnalysis").ignore(); // call before opening first file
  
  // tell EventLoop about our output:
  EL::OutputStream out("output");
  job.outputAdd(out);
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODBaselineAnalysis::histInitialize()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  
  h_events = new TH1D("events", "Events", 6, 0.5, 6.5);
  h_events->GetXaxis()->SetBinLabel(1, "# events initial");
  h_events->GetXaxis()->SetBinLabel(2, "# events selected");
  h_events->GetXaxis()->SetBinLabel(3, "sumw initial");
  h_events->GetXaxis()->SetBinLabel(4, "sumw selected");
  h_events->GetXaxis()->SetBinLabel(5, "sumw2 initial");
  h_events->GetXaxis()->SetBinLabel(6, "sumw2 selected");

  h_cutflow = new TH1D("cutflow", "Cutflow", 8, 0.5, 8.5);
  h_cutflow->GetXaxis()->SetBinLabel(1, "All");
  h_cutflow->GetXaxis()->SetBinLabel(2, "GRL/MC filter");
  h_cutflow->GetXaxis()->SetBinLabel(3, "DQ");
  h_cutflow->GetXaxis()->SetBinLabel(4, "Trigger");
  h_cutflow->GetXaxis()->SetBinLabel(5, "Good Vertex");
  h_cutflow->GetXaxis()->SetBinLabel(6, "Cosmic Muon (not applied)");
  h_cutflow->GetXaxis()->SetBinLabel(7, "Bad Jet Cleaning");
  h_cutflow->GetXaxis()->SetBinLabel(8, "Skim (1 baseline photon) ");

  h_cutflow_w = new TH1D("cutflow_w", "Cutflow_w", 8, 0.5, 8.5);
  h_cutflow_w->GetXaxis()->SetBinLabel(1, "All");
  h_cutflow_w->GetXaxis()->SetBinLabel(2, "GRL/MC filter");
  h_cutflow_w->GetXaxis()->SetBinLabel(3, "DQ");
  h_cutflow_w->GetXaxis()->SetBinLabel(4, "Trigger");
  h_cutflow_w->GetXaxis()->SetBinLabel(5, "Good Vertex");
  h_cutflow_w->GetXaxis()->SetBinLabel(6, "Cosmic Muon (not applied)");
  h_cutflow_w->GetXaxis()->SetBinLabel(7, "Bad Jet Cleaning");
  h_cutflow_w->GetXaxis()->SetBinLabel(8, "Skim (1 baseline photon) ");

 TDirectory *out_dir = (TDirectory*) wk()->getOutputFile("output");
 h_events->SetDirectory(out_dir);
 h_cutflow->SetDirectory(out_dir);
 h_cutflow_w->SetDirectory(out_dir);
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODBaselineAnalysis::fileExecute()
{
  // Here you do everything that n<eeds to be done exactly once for every
  // single file, e.g. collect a st of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODBaselineAnalysis::changeInput(bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  // Check file's metadata:    
  m_event = wk()->xaodEvent(); 
  
  TTree *metadata = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
  if (!metadata) {
    Error(APP_NAME, "MetaData tree not found! Exiting.");
    return EL::StatusCode::FAILURE;
  }
  metadata->LoadTree(0);

  //check if file is from a DxAOD
  bool is_derivation = !metadata->GetBranch("StreamAOD");

  ULong64_t m_initial_events = 0;
  Double_t m_initial_sumw = 0.;
  Double_t m_initial_sumw2 = 0.;

  if (is_derivation) {
    //Read the CutBookkeeper container
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    if (!m_event->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()) {
      Error(APP_NAME, "Failed to retrieve CutBookkeepers from MetaData! Exiting.");
      return EL::StatusCode::FAILURE;
    }
    
    // Now, let's actually find the right one that contains all the needed info...
    const xAOD::CutBookkeeper* all_events_cbk = 0;
    //const xAOD::CutBookkeeper* dxaod_events_cbk = 0;
      
    int maxCycle = -1;
    for (const auto& cbk :  *completeCBC) {
      if (cbk->cycle() > maxCycle && cbk->name() == "AllExecutedEvents" && cbk->inputStream() == "StreamAOD") {
        all_events_cbk = cbk;
        maxCycle = cbk->cycle();
      }
    }
      
    m_initial_events = all_events_cbk->nAcceptedEvents();
    m_initial_sumw   = all_events_cbk->sumOfEventWeights();
    m_initial_sumw2  = all_events_cbk->sumOfEventWeightsSquared();
  }
  else {
    TTree* CollectionTree = dynamic_cast<TTree*>( wk()->inputFile()->Get("CollectionTree") );
    
    m_initial_events  = CollectionTree->GetEntries(); 
    m_initial_sumw    = CollectionTree->GetWeight() * CollectionTree->GetEntries();
    m_initial_sumw2   = (CollectionTree->GetWeight() * CollectionTree->GetWeight()) * CollectionTree->GetEntries();
  }
  
  std::cout << "Initial events = " << m_initial_events << ", Sumw = " << m_initial_sumw << std::endl;
  
  h_events->Fill(1, m_initial_events);
  h_events->Fill(3, m_initial_sumw);
  h_events->Fill(5, m_initial_sumw2);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODBaselineAnalysis::initialize()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  m_event = wk()->xaodEvent();

  // Data dir
  m_data_dir = gSystem->ExpandPathName("$ROOTCOREBIN/data/PhotonMetNtuple/");
  
  // Config file
  if (config_file.empty()) {
    Error(APP_NAME, "No configfile. Exiting... ");
    return EL::StatusCode::FAILURE;
  }
  config_file = m_data_dir + config_file;

  ReadConfiguration();
  DumpConfiguration();

  // ST Options
  ST::ISUSYObjDef_xAODTool::DataSource datasource = (is_data ? ST::ISUSYObjDef_xAODTool::Data : (is_atlfast ? ST::ISUSYObjDef_xAODTool::AtlfastII : ST::ISUSYObjDef_xAODTool::FullSim));
  
  susytools = new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");
  susytools->msg().setLevel(MSG::FATAL);

  CHECK(susytools->setProperty("DataSource", datasource)); 

  // ST Config file
  CHECK(susytools->setProperty("ConfigFile", m_st_config_file));

  // Pile Up Reweighting
  CHECK(susytools->setProperty("PRWConfigFiles", m_prw_mc_files));
  CHECK(susytools->setProperty("PRWLumiCalcFiles", m_prw_lumicalc_files));
 
  if (susytools->initialize() != StatusCode::SUCCESS) {
    Error(APP_NAME, "Cannot intialize SUSYObjDef_xAOD...");
    Error(APP_NAME, "Exiting... ");
    return EL::StatusCode::FAILURE;
  } 
  else {
    Info(APP_NAME, "SUSYObjDef_xAOD initialized... ");
  }
  
  // GRL
  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  CHECK(m_grl->setProperty("GoodRunsListVec", m_grl_files));
  CHECK(m_grl->setProperty("PassThrough", false) ); // if true (default) will ignore result of GRL and will just pass all events
  if (!m_grl->initialize().isSuccess()) { // check this isSuccess
    Error(APP_NAME, "Failed to properly initialize the GRL. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  TDirectory *out_dir = (TDirectory*) wk()->getOutputFile("output");
  
  outtree = new BaselineTree("BaselineTree");

  //Here you can select a subsect of the systematic uncertainties
  CHECK(outtree->setProperty("OutFile", out_dir));
  CHECK(outtree->setProperty("IsMC", !is_data));
  CHECK(outtree->initialize());

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODBaselineAnalysis::execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // clear tree
  outtree->clear();
  
  //----------------------------
  // Event information
  //---------------------------
 
  const xAOD::EventInfo *eventInfo = 0;
  if (!m_event->retrieve(eventInfo, "EventInfo").isSuccess()) {
    Error(APP_NAME, "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  if (eventInfo->eventNumber() != 1774431567)
    return EL::StatusCode::SUCCESS;
  
  std::cout << wk()->inputFile()->GetName() << std::endl;

  outtree->SetEventNumber(eventInfo->eventNumber());
  outtree->SetRunNumber(eventInfo->runNumber());

  // check if the event is data or MC
  bool is_mc = !is_data; 

  if (is_data && eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION)) {
    Error(APP_NAME, "Bad configuration. Is DATA or MC?. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // float mc_weight = 1.;
  // if (is_mc) {
  //   mc_weight = eventInfo->mcEventWeight();
  //   outtree->SetWeightMC(mc_weight);

  CHECK(susytools->ApplyPRWTool());  
  //   outtree->SetWeightPU(susytools->GetPileupWeight());
  //   outtree->SetPRWHash(susytools->GetPileupWeightHash());
  // }
  // else {
  //   outtree->SetWeightMC(1.);
  //   outtree->SetWeightPU(1.);
  // }
  
  //--------------------
  // CLEANING CUTS HERE
  //--------------------
  h_cutflow->Fill(1);

  if (is_data && !m_grl->passRunLB(*eventInfo)) {
    return EL::StatusCode::SUCCESS; // go to next event
  }
  if (is_mc && !mc_filter->accept_event(eventInfo->mcChannelNumber(), *m_event)) {
    return EL::StatusCode::SUCCESS; // go to next event
  }
  h_cutflow->Fill(2);

  if (is_data &&
      ((eventInfo->errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error) ||
       (eventInfo->errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error) ||
       (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18)))) {
    return EL::StatusCode::SUCCESS; // go to next event
  }
  h_cutflow->Fill(3);
  

  // Trigger
  bool passed = susytools->IsTrigPassed("HLT_g140_loose") || susytools->IsTrigPassed("HLT_g120_loose");
  //float prescale = susytools->GetTrigPrescale(trigitem);
  //Info(APP_NAME, "passing %s trigger? %d, prescale %f", trigitem.c_str(), (int)passed, prescale);
  if (!passed)
    return EL::StatusCode::SUCCESS;
  h_cutflow->Fill(4);


  //Get the event Vertex
  const xAOD::Vertex *prim_vx = 0;
  prim_vx = susytools->GetPrimVtx();
  if (!prim_vx) 
    return EL::StatusCode::SUCCESS;
  h_cutflow->Fill(5);

  //---------------------------------
  // RETRIEVE THE NOMINAL CONTAINERS
  //---------------------------------
  
  // Electrons
  xAOD::ElectronContainer* electrons_nominal(0);
  xAOD::ShallowAuxContainer* electrons_nominal_aux(0);
  CHECK(susytools->GetElectrons(electrons_nominal, electrons_nominal_aux));
  
  //  electrons_nominal->sort(ptsorter);

  // Photons
  xAOD::PhotonContainer* photons_nominal(0);
  xAOD::ShallowAuxContainer* photons_nominal_aux(0);
  CHECK(susytools->GetPhotons(photons_nominal, photons_nominal_aux));

  //  photons_nominal->sort(ptsorter);

  // Muons
  xAOD::MuonContainer* muons_nominal(0);
  xAOD::ShallowAuxContainer* muons_nominal_aux(0);
  CHECK(susytools->GetMuons(muons_nominal, muons_nominal_aux));

  //  muons_nominal->sort(ptsorter);
  
  // Jets
  xAOD::JetContainer* jets_nominal(0);
  xAOD::ShallowAuxContainer* jets_nominal_aux(0);
  CHECK(susytools->GetJets(jets_nominal, jets_nominal_aux));

  //  jets_nominal->sort(ptsorter);

  
  // MET (remember,you can pick either CST or TST)
  xAOD::MissingETContainer* met_nominal = new xAOD::MissingETContainer;
  xAOD::MissingETAuxContainer* met_nominal_aux = new xAOD::MissingETAuxContainer;
  met_nominal->setStore(met_nominal_aux);

  //-------------------------
  // NOMINAL TREE PROCESSING
  //-------------------------
  
  // Overlap removal
  CHECK(susytools->OverlapRemoval(electrons_nominal, muons_nominal, jets_nominal, photons_nominal));
  
  // MET
  CHECK(susytools->GetMET(*met_nominal, jets_nominal, electrons_nominal, muons_nominal, photons_nominal));
  

  h_cutflow->Fill(7);

  AnalysisBaselineCollections collections;
  collections.photons = photons_nominal;
  collections.electrons = electrons_nominal;
  collections.muons = muons_nominal;
  collections.jets = jets_nominal;
  collections.met = met_nominal;
  collections.photons_aux = photons_nominal_aux;
  collections.electrons_aux = electrons_nominal_aux;
  collections.muons_aux = muons_nominal_aux;   
  collections.jets_aux = jets_nominal_aux;
  collections.met_aux = met_nominal_aux;
  
  outtree->process(collections);
  
  h_cutflow->Fill(8);

  CHECK(outtree->FillTree());
  
  delete jets_nominal;
  delete jets_nominal_aux;
  delete muons_nominal;
  delete muons_nominal_aux;
  delete electrons_nominal;
  delete electrons_nominal_aux;
  delete photons_nominal;
  delete photons_nominal_aux;
  delete met_nominal;
  delete met_nominal_aux;
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODBaselineAnalysis::postExecute()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODBaselineAnalysis::finalize()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.
  
  // GRL
  if (m_grl) {
    delete m_grl;
    m_grl = 0;
  }
  
  if (susytools) {
    delete susytools;
    susytools = 0;
  }
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODBaselineAnalysis::histFinalize()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}

std::vector<std::string> xAODBaselineAnalysis::SplitString(TString line){
  
  std::vector<std::string> vtokens;
  TObjArray* tokens = TString(line).Tokenize(","); 
  
  if(tokens->GetEntriesFast()) {
    TIter iString(tokens);
    TObjString* os=0;
    while ((os=(TObjString*)iString())) {
      vtokens.push_back( os->GetString().Data() );
    }
  }
  delete tokens;
  return vtokens;
}

void xAODBaselineAnalysis::ReadConfiguration()
{
  TEnv env(config_file.c_str());

  m_st_config_file = m_data_dir + env.GetValue("ST.ConfigFile", "");
  
  TString ilumicalc_files = env.GetValue("PRW.LumiCalcFile", "");
  for (auto s : SplitString(ilumicalc_files))
    m_prw_lumicalc_files.push_back(m_data_dir + s);
        
  TString mc_files = env.GetValue("PRW.MCFile", "");
  for (auto s : SplitString(mc_files))
    m_prw_mc_files.push_back(m_data_dir + s);

  TString grl_files = env.GetValue("GRL.File", "");
  for (auto s : SplitString(grl_files))
    m_grl_files.push_back(m_data_dir + s);

}

void xAODBaselineAnalysis::DumpConfiguration()
{
  Info(APP_NAME, "-- DumpConfiguration");

  Info(APP_NAME, "Using configfile: %s", config_file.c_str());

  Info(APP_NAME, "ST.ConfigFile: %s", m_st_config_file.c_str());
  
  for (auto s : m_prw_lumicalc_files)
    Info(APP_NAME, "PRW.LumiCalcFile: %s", s.c_str());

  for (auto s : m_prw_mc_files)
    Info(APP_NAME, "PRW.MCFile: %s", s.c_str());

  for (auto s : m_grl_files)
    Info(APP_NAME, "GRL.File: %s", s.c_str());
}

