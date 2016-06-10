#include <TSystem.h>

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <PhotonMetNtuple/xAODAnalysis.h>
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

// GRL
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

//PU Reweighting
#include "PileupReweighting/PileupReweightingTool.h"

#include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"

// Amg include
#include "EventPrimitives/EventPrimitivesHelpers.h"
#include "xAODEgamma/EgammaxAODHelpers.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

const char *APP_NAME = "PhotonMetNtuple";
const char *APP_VERSION = "v31";


// this is needed to distribute the algorithm to the workers
ClassImp(xAODAnalysis)


bool ptsorter( const xAOD::IParticle* j1, const xAOD::IParticle* j2 ) {
  return ( j1->pt() > j2->pt() );
}

xAODAnalysis::xAODAnalysis() : 
  my_XsecDB(0), m_grl(0), susytools(0),
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

EL::StatusCode xAODAnalysis::setupJob(EL::Job &job)
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
  xAOD::Init("xAODAnalysis").ignore(); // call before opening first file
  
  // tell EventLoop about our output:
  EL::OutputStream out("output");
  job.outputAdd(out);
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis::histInitialize()
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

 TDirectory *out_dir = (TDirectory*) wk()->getOutputFile("output");
 h_events->SetDirectory(out_dir);
 h_cutflow->SetDirectory(out_dir);
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis::fileExecute()
{
  // Here you do everything that n<eeds to be done exactly once for every
  // single file, e.g. collect a st of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis::changeInput(bool firstFile)
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

  if (!is_data) {

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
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis::initialize()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  Info(APP_NAME, APP_VERSION);

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
  // std::string st_config_file = m_data_dir+"ST_PhotonMet.conf";
  CHECK(susytools->setProperty("ConfigFile", m_st_config_file));

  // Pile Up Reweighting
  //std::vector<std::string> prwFiles;
  //prwFiles.push_back(data_dir+"merged_prw_mc15c.root");
  //  std::vector<std::string> lumicalcFiles;
  //  lumicalcFiles.push_back(data_dir+"ilumicalc_histograms_None_276262-284154.root");
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
  
  // my_XsecDB = new SUSY::CrossSectionDB("susy_crosssections_13TeV.txt");
  
  // GRL
  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  // std::vector<std::string> vecStringGRL;
  //vecStringGRL.push_back(data_dir +  "grl.xml");
  CHECK(m_grl->setProperty("GoodRunsListVec", m_grl_files));
  CHECK(m_grl->setProperty("PassThrough", false) ); // if true (default) will ignore result of GRL and will just pass all events
  if (!m_grl->initialize().isSuccess()) { // check this isSuccess
    Error(APP_NAME, "Failed to properly initialize the GRL. Exiting." );
    return EL::StatusCode::FAILURE;
  }


  // Now we can look at systematics:    
  if (do_syst) {
    systInfoList = susytools->getSystInfoList();
  }
  else {
    ST::SystInfo infodef;
    infodef.affectsKinematics = false;
    infodef.affectsWeights = false;
    infodef.affectsType = ST::Unknown;
    systInfoList.push_back(infodef);
  } 

  TDirectory *out_dir = (TDirectory*) wk()->getOutputFile("output");
  
  outtree = new MiniTree2("MiniTree");

  //Here you can select a subsect of the systematic uncertainties
  CHECK(outtree->setProperty("SystematicList", systInfoList));
  CHECK(outtree->setProperty("OutFile", out_dir));
  CHECK(outtree->setProperty("IsMC", !is_data));
  CHECK(outtree->initialize());

  mc_filter = new MCFilter;

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis :: execute ()
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
  
  outtree->SetEventNumber(eventInfo->eventNumber());
  outtree->SetRunNumber(eventInfo->runNumber());
  outtree->SetAvgMu(eventInfo->averageInteractionsPerCrossing());

  // check if the event is data or MC
  bool is_mc = !is_data; 

  if (is_data && eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION)) {
    Error(APP_NAME, "Bad configuration. Is DATA or MC?. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  if (is_mc) {
    outtree->SetWeightMC(eventInfo->mcEventWeight());

    CHECK(susytools->ApplyPRWTool());  
	outtree->SetWeightPU(susytools->GetPileupWeight());
    outtree->SetPRWHash(susytools->GetPileupWeightHash());
  }
  else {
    outtree->SetWeightMC(1.);
    outtree->SetWeightPU(1.);
  }
  
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
  
  electrons_nominal->sort(ptsorter);

  // Photons
  xAOD::PhotonContainer* photons_nominal(0);
  xAOD::ShallowAuxContainer* photons_nominal_aux(0);
  CHECK(susytools->GetPhotons(photons_nominal, photons_nominal_aux));

  photons_nominal->sort(ptsorter);

  // Muons
  xAOD::MuonContainer* muons_nominal(0);
  xAOD::ShallowAuxContainer* muons_nominal_aux(0);
  CHECK(susytools->GetMuons(muons_nominal, muons_nominal_aux));

  muons_nominal->sort(ptsorter);
  
  // Jets
  xAOD::JetContainer* jets_nominal(0);
  xAOD::ShallowAuxContainer* jets_nominal_aux(0);
  CHECK(susytools->GetJets(jets_nominal, jets_nominal_aux));

  jets_nominal->sort(ptsorter);

  
  // MET (remember,you can pick either CST or TST)
  xAOD::MissingETContainer* met_nominal = new xAOD::MissingETContainer;
  xAOD::MissingETAuxContainer* met_nominal_aux = new xAOD::MissingETAuxContainer;
  met_nominal->setStore(met_nominal_aux);

  //--------------------
  // LOOP ON SYSTEMATICS
  //--------------------
  int ret = 0;
  for (const auto& sysInfo : systInfoList) {
      
    if (sysInfo.affectsKinematics || sysInfo.affectsWeights) {
      const CP::SystematicSet& sys = sysInfo.systset;
      if (susytools->applySystematicVariation(sys) != CP::SystematicCode::Ok) {
        Error(APP_NAME, "Cannot configure SUSYTools for systematic var. %s", (sys.name()).c_str() );
      }
      else {
        //Info(APP_NAME, "SUSYTools configured for systematic var. %s", (sys.name()).c_str() );

        // Generic pointers for either nominal or systematics copy
        xAOD::ElectronContainer* electrons(electrons_nominal);
        xAOD::PhotonContainer* photons(photons_nominal);
        xAOD::MuonContainer* muons(muons_nominal);
        xAOD::JetContainer* jets(jets_nominal);
        xAOD::MissingETContainer* met(met_nominal);
        
        // Aux containers too
        xAOD::ShallowAuxContainer* electrons_aux(electrons_nominal_aux);
        xAOD::ShallowAuxContainer* photons_aux(photons_nominal_aux);
        xAOD::ShallowAuxContainer* muons_aux(muons_nominal_aux);
        xAOD::ShallowAuxContainer* jets_aux(jets_nominal_aux);
        xAOD::MissingETAuxContainer* met_aux(met_nominal_aux);
        
        // If necessary (kinematics affected), make a shallow copy with the variation applied
        bool syst_affectsElectrons = ST::testAffectsObject(xAOD::Type::Electron, sysInfo.affectsType);
        bool syst_affectsMuons = ST::testAffectsObject(xAOD::Type::Muon, sysInfo.affectsType);
        bool syst_affectsTaus = ST::testAffectsObject(xAOD::Type::Tau, sysInfo.affectsType);
        bool syst_affectsPhotons = ST::testAffectsObject(xAOD::Type::Photon, sysInfo.affectsType);
        bool syst_affectsJets = ST::testAffectsObject(xAOD::Type::Jet, sysInfo.affectsType);
        bool syst_affectsBTag = ST::testAffectsObject(xAOD::Type::BTag, sysInfo.affectsType);
        
        if (syst_affectsTaus) 
          continue;
        
        if (syst_affectsElectrons) {
          xAOD::ElectronContainer* electrons_syst(0);
          xAOD::ShallowAuxContainer* electrons_syst_aux(0);

          CHECK(susytools->GetElectrons(electrons_syst, electrons_syst_aux));
          electrons_syst->sort(ptsorter);

          electrons = electrons_syst;
          electrons_aux = electrons_syst_aux;
        }
		
        if (syst_affectsMuons) {
          xAOD::MuonContainer* muons_syst(0);
          xAOD::ShallowAuxContainer* muons_syst_aux(0);

          CHECK(susytools->GetMuons(muons_syst, muons_syst_aux));
          muons_syst->sort(ptsorter);

          muons = muons_syst;
          muons_aux = muons_syst_aux;
        }
		
        if (syst_affectsPhotons) {
          xAOD::PhotonContainer* photons_syst(0);
          xAOD::ShallowAuxContainer* photons_syst_aux(0);

          CHECK(susytools->GetPhotons(photons_syst, photons_syst_aux));
          photons_syst->sort(ptsorter);

          photons = photons_syst;
          photons_aux = photons_syst_aux;
        }
		
        if (syst_affectsJets || syst_affectsBTag) {
          xAOD::JetContainer* jets_syst(0);
          xAOD::ShallowAuxContainer* jets_syst_aux(0);

          CHECK(susytools->GetJets(jets_syst, jets_syst_aux));
          jets_syst->sort(ptsorter);

          jets = jets_syst;
          jets_aux = jets_syst_aux;
        }
        
        //-----------------
        // OVERLAP REMOVAL
        //-----------------
        CHECK(susytools->OverlapRemoval(electrons, muons, jets, photons));
        
        xAOD::MissingETContainer*    met_syst = new xAOD::MissingETContainer;
        xAOD::MissingETAuxContainer* met_syst_aux = new xAOD::MissingETAuxContainer;
        met_syst->setStore(met_syst_aux);

        // MET
        CHECK(susytools->GetMET(*met_syst, jets, electrons, muons, photons));
        
        met = met_syst;
        met_aux = met_syst_aux;
        
        for (const auto& el : *electrons) {
          if (is_mc && el->auxdata<char>("baseline") == 1 &&
              el->auxdata<char>("passOR") == 1 &&
              el->auxdata<char>("signal") == 1)
            susytools->GetSignalElecSF(*el);
        }

        for (const auto& ph : *photons) {
          if (is_mc && ph->auxdata<char>("baseline") == 1 &&
              ph->auxdata<char>("passOR") == 1 &&
              ph->auxdata<char>("signal") == 1)
            susytools->GetSignalPhotonSF(*ph, true, false);
        }
        
        bool skip = false;
        for (const auto& mu : *muons) {
          if (is_mc && mu->auxdata<char>("baseline") == 1 &&
              mu->auxdata<char>("passOR") == 1 && 
              mu->auxdata<char>("signal") == 1)
            susytools->GetSignalMuonSF(*mu);
        }
        
        for (const auto& jet : *jets) {
          if (jet->auxdata<char>("baseline") == 1 &&
              jet->auxdata<char>("passOR") == 1 &&
              jet->auxdata<char>("bad") == 1) {
            skip = true;
            break;
          }
        }
        
        if (is_mc) {
          outtree->SetWeightBtag(sys.name().c_str(), susytools->BtagSF(jets));
        }

        if (!skip) {
          AnalysisCollections collections;

          collections.photons = photons;
          collections.electrons = electrons;
          collections.muons = muons;
          collections.jets = jets;
          collections.met = met;
          collections.photons_aux = photons_aux;
          collections.electrons_aux = electrons_aux;
          collections.muons_aux = muons_aux;   
          collections.jets_aux = jets_aux;
          collections.met_aux = met_aux;
          
          ret += outtree->process(collections, (sys.name()).c_str());
        }
	  
        if (syst_affectsElectrons) {
          delete electrons;
          delete electrons_aux;
        }
        if (syst_affectsMuons) {
          delete muons;
          delete muons_aux;
        }
        if (syst_affectsPhotons) {
          delete photons;
          delete photons_aux;
        }
        if (syst_affectsJets || syst_affectsBTag) {
          delete jets;
          delete jets_aux;
        } 
        delete met;
        delete met_aux;
	  
      }
    }
  } //end loop over systematics affecting kinematics or weights


  //-------------------------
  // NOMINAL TREE PROCESSING
  //-------------------------
  
  // Overlap removal
  CHECK(susytools->OverlapRemoval(electrons_nominal, muons_nominal, jets_nominal, photons_nominal));
  
  // MET
  CHECK(susytools->GetMET(*met_nominal, jets_nominal, electrons_nominal, muons_nominal, photons_nominal));
  
  // FIX: TST cleaning, until bug fix
  if (susytools->passTSTCleaning(*met_nominal))
    outtree->SetPassTSTCleaning(1);
  else
    outtree->SetPassTSTCleaning(0);

  // electrons
  for (const auto& el : *electrons_nominal) {
    if (is_mc && el->auxdata<char>("baseline") == 1 &&
        el->auxdata<char>("passOR") == 1 &&
        el->auxdata<char>("signal") == 1)
      susytools->GetSignalElecSF(*el);
  }
  
  // photons
  for (const auto& ph : *photons_nominal) {
    if (is_mc && ph->auxdata<char>("baseline") == 1 &&
        ph->auxdata<char>("passOR") == 1 && 
        ph->auxdata<char>("signal") == 1)
      susytools->GetSignalPhotonSF(*ph, true, false);
  }
  
  bool skip = false;
  // muons
  for (const auto& mu : *muons_nominal) {
    if (is_mc && mu->auxdata<char>("baseline") == 1 &&
        mu->auxdata<char>("passOR") == 1 && 
        mu->auxdata<char>("signal") == 1)
      susytools->GetSignalMuonSF(*mu);
  }
  if (!skip)
    h_cutflow->Fill(6);

  // Badj jet veto
  for (const auto& jet : *jets_nominal) {  
    if (jet->auxdata<char>("baseline") == 1 &&
        jet->auxdata<char>("passOR") == 1 &&
        jet->auxdata<char>("bad") == 1) {
      skip = true;
      break;
    }
  }
  
  if (is_mc)
    outtree->SetWeightBtag("Nominal", susytools->BtagSF(jets_nominal));


  // MC final state
  if (is_mc && is_susy_ewk) {
    
    const xAOD::TruthParticleContainer* truth_particles = 0;
    if(!m_event->retrieve(truth_particles, "TruthParticles").isSuccess()) {
      Error(APP_NAME, "Failed to retrieve truth particles collection. Exiting." );
      return EL::StatusCode::FAILURE;
    }
  
    int pdg1, pdg2;
    susytools->FindSusyHardProc(truth_particles, pdg1, pdg2);

    outtree->SetMCFinalState(SUSY::finalState(pdg1, pdg2));
  }


  if (!skip) {
    h_cutflow->Fill(7);

    AnalysisCollections collections;
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
    
    ret += outtree->process(collections, "Nominal");
  }
  
  if (ret > 0) {
    h_cutflow->Fill(8);
    CHECK(outtree->FillTree());
  }

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

EL::StatusCode xAODAnalysis::postExecute()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis::finalize()
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
  
  if (my_XsecDB) {
    delete my_XsecDB;
    my_XsecDB = 0;
  }
  
  if (susytools) {
    delete susytools;
    susytools = 0;
  }
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis::histFinalize()
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

std::vector<std::string> xAODAnalysis::SplitString(TString line){
  
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

void xAODAnalysis::ReadConfiguration()
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

void xAODAnalysis::DumpConfiguration()
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
