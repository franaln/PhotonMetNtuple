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

#include "PhotonMetNtuple/OutTree.h"

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

// this is needed to distribute the algorithm to the workers
ClassImp(xAODAnalysis)

xAODAnalysis::xAODAnalysis() : 
my_XsecDB(0), m_grl(0), objTool(0)
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
  h_cutflow->GetXaxis()->SetBinLabel(2, "GRL");
  h_cutflow->GetXaxis()->SetBinLabel(3, "DQ");
  h_cutflow->GetXaxis()->SetBinLabel(4, "Trigger");
  h_cutflow->GetXaxis()->SetBinLabel(5, "Good Vertex");
  h_cutflow->GetXaxis()->SetBinLabel(6, "Cosmic Muon (not applied)");
  h_cutflow->GetXaxis()->SetBinLabel(7, "Bad Jet");
  h_cutflow->GetXaxis()->SetBinLabel(8, "Skim");
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis::fileExecute()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis::changeInput(bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  
  const char *APP_NAME = "changeInput()";
  isDerived = false;

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
    const xAOD::CutBookkeeper* dxaod_events_cbk = 0;

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

  std::cout << m_initial_events << " " << m_initial_sumw << " " << m_initial_sumw2 << std::endl;
  h_events->Fill(1, m_initial_events);
  h_events->Fill(3, m_initial_sumw);
  h_events->Fill(5, m_initial_sumw2);

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

  const char *APP_NAME = "xAODAnalaysis";

  const char *APP_VERSION = "PhotonMetNtuple v16";
 
  Info(APP_NAME, APP_VERSION);

  m_event = wk()->xaodEvent();

  // ST Options
  // ST::SettingDataSource datasource = isData ? ST::Data : (isAtlfast ? ST::AtlfastII : ST::FullSim);
  ST::ISUSYObjDef_xAODTool::DataSource datasource = (isData ? ST::ISUSYObjDef_xAODTool::Data : (isAtlfast ? ST::ISUSYObjDef_xAODTool::AtlfastII : ST::ISUSYObjDef_xAODTool::FullSim));

  objTool = new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");
  objTool->msg().setLevel(MSG::FATAL);

  CHECK(objTool->setProperty("DataSource", datasource)); 

  // Data dir
  std::string data_dir = gSystem->ExpandPathName("$ROOTCOREBIN/data/PhotonMetNtuple/");

  // Config file
  CHECK(objTool->setProperty("ConfigFile", data_dir+"SUSYTools_Default.conf"));


  // Pile Up Reweighting
  // std::vector<std::string> prwFiles;
  // prwFiles.push_back(data_dir+"signal_prw.root"); 
  // prwFiles.push_back(data_dir+"signal_ewk_prw.root"); 
  // prwFiles.push_back(data_dir+"merged_prw.root"); 
  // prwFiles.push_back(data_dir+"ttgam_prw.root"); 
  // CHECK( objTool->setProperty("PRWConfigFiles", prwFiles) );

  // std::vector<std::string> lumicalcFiles;
  // lumicalcFiles.push_back(data_dir+"ilumicalc_histograms_None_276262-282712.root");
  // CHECK(objTool->setProperty("PRWLumiCalcFiles", lumicalcFiles));
 
  
  if (objTool->initialize() != StatusCode::SUCCESS) {
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
  std::vector<std::string> vecStringGRL;
  vecStringGRL.push_back(data_dir +  "grl.xml");
  CHECK(m_grl->setProperty("GoodRunsListVec", vecStringGRL));
  CHECK(m_grl->setProperty("PassThrough", false) ); // if true (default) will ignore result of GRL and will just pass all events
  if (!m_grl->initialize().isSuccess()) { // check this isSuccess
    Error(APP_NAME, "Failed to properly initialize the GRL. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  
  
  // Now we can look at systematics:    
  doSyst = false;
  if (!doSyst) {
    ST::SystInfo infodef;
    infodef.affectsKinematics = false;
    infodef.affectsWeights = false;
    infodef.affectsType = ST::Unknown;
    systInfoList.push_back(infodef);
  } 
  else {
    systInfoList = objTool->getSystInfoList();
  }

  TDirectory *out_dir = (TDirectory*) wk()->getOutputFile("output");
  
  outtree = new OutTree("outtree");

  //Here you can select a subsect of the systematic uncertainties
  CHECK(outtree->setProperty("SystematicList", systInfoList));
  CHECK(outtree->setProperty("OutFile", out_dir));
  CHECK(outtree->setProperty("IsMC", !isData));
  CHECK(outtree->initialize());

  h_events->SetDirectory(out_dir);
  h_cutflow->SetDirectory(out_dir);
    
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODAnalysis :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  const char *APP_NAME = "xAODAnalysis::execute()";
  
  //----------------------------
  // Event information
  //---------------------------
  
  const xAOD::EventInfo *eventInfo = 0;
  if (!m_event->retrieve(eventInfo, "EventInfo").isSuccess()) {
    Error(APP_NAME, "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  event_number = eventInfo->eventNumber();
  run_number   = eventInfo->runNumber();
  avg_mu       = eventInfo->averageInteractionsPerCrossing();

  // check if the event is data or MC
  bool is_mc = false;
  if (eventInfo->eventType(xAOD::EventInfo::IS_SIMULATION)) {
    is_mc = true;

    // CHECK(objTool->ApplyPRWTool());  

	// m_xsec = my_XsecDB->xsectTimesEff(eventInfo->mcChannelNumber());
    weight_mc = (eventInfo->mcEventWeights()).at(0);
    weight_pu = 1.;//objTool->GetPileupWeight();
  }
  else {
    weight_mc = 1.;
    weight_pu = 1.;
	//m_xsec = 1.;
  }
  
  //--------------------
  // CLEANING CUTS HERE
  //--------------------
  h_cutflow->Fill(1);

  if (!is_mc && !m_grl->passRunLB(*eventInfo)) {
    return EL::StatusCode::SUCCESS; // go to next event
  }
  h_cutflow->Fill(2);

  if (!is_mc &&
      ((eventInfo->errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error ) ||
       (eventInfo->errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error ) ||
       (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18)))) {
    return EL::StatusCode::SUCCESS; // go to next event
  }
  h_cutflow->Fill(3);
  

  // Trigger
  bool passed = objTool->IsTrigPassed("HLT_g140_loose") || objTool->IsTrigPassed("HLT_g120_loose");
  //float prescale = objTool->GetTrigPrescale(trigitem);
  //Info(APP_NAME, "passing %s trigger? %d, prescale %f", trigitem.c_str(), (int)passed, prescale);
  if (!is_mc && !passed)
    return EL::StatusCode::SUCCESS;
  h_cutflow->Fill(4);


  //Get the event Vertex
  const xAOD::Vertex *prim_vx = 0;
  prim_vx = objTool->GetPrimVtx();
  if (!prim_vx) 
    return EL::StatusCode::SUCCESS;
  h_cutflow->Fill(5);
    
  //---------------------------------
  // RETRIEVE THE NOMINAL CONTAINERS
  //---------------------------------
  
  // Electrons
  xAOD::ElectronContainer* electrons_nominal(0);
  xAOD::ShallowAuxContainer* electrons_nominal_aux(0);
  CHECK(objTool->GetElectrons(electrons_nominal, electrons_nominal_aux));
  
  // Photons
  xAOD::PhotonContainer* photons_nominal(0);
  xAOD::ShallowAuxContainer* photons_nominal_aux(0);
  CHECK(objTool->GetPhotons(photons_nominal, photons_nominal_aux));

  // Muons
  xAOD::MuonContainer* muons_nominal(0);
  xAOD::ShallowAuxContainer* muons_nominal_aux(0);
  CHECK(objTool->GetMuons(muons_nominal, muons_nominal_aux));
  
  // Jets
  xAOD::JetContainer* jets_nominal(0);
  xAOD::ShallowAuxContainer* jets_nominal_aux(0);
  CHECK(objTool->GetJets(jets_nominal, jets_nominal_aux));
  
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
      if (objTool->applySystematicVariation(sys) != CP::SystematicCode::Ok) {
        Error(APP_NAME, "Cannot configure SUSYTools for systematic var. %s", (sys.name()).c_str() );
      }
      else {
        Info(APP_NAME, "SUSYTools configured for systematic var. %s", (sys.name()).c_str() );
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
          CHECK(objTool->GetElectrons(electrons_syst,electrons_syst_aux));
          electrons = electrons_syst;
          electrons_aux = electrons_syst_aux;
        }
		
        if (syst_affectsMuons) {
          xAOD::MuonContainer* muons_syst(0);
          xAOD::ShallowAuxContainer* muons_syst_aux(0);
          CHECK(objTool->GetMuons(muons_syst,muons_syst_aux));
          muons = muons_syst;
          muons_aux = muons_syst_aux;
        }
		
        // if (syst_affectsTaus) {
        //   xAOD::TauJetContainer* taus_syst(0);
        //   xAOD::ShallowAuxContainer* taus_syst_aux(0);
        //   CHECK(objTool->GetTaus(taus_syst,taus_syst_aux));
        //   taus = taus_syst;
        //   taus_aux = taus_syst_aux;
        // }
        
        if (syst_affectsPhotons) {
          xAOD::PhotonContainer* photons_syst(0);
          xAOD::ShallowAuxContainer* photons_syst_aux(0);
          CHECK(objTool->GetPhotons(photons_syst,photons_syst_aux));
          photons = photons_syst;
          photons_aux = photons_syst_aux;
        }
		
        if (syst_affectsJets || syst_affectsBTag) {
          xAOD::JetContainer* jets_syst(0);
          xAOD::ShallowAuxContainer* jets_syst_aux(0);
          CHECK(objTool->GetJets(jets_syst,jets_syst_aux));
          jets = jets_syst;
          jets_aux = jets_syst_aux;
        }
        
        //-----------------
        // OVERLAP REMOVAL
        //-----------------
        float weight_sf = 1.;
        CHECK(objTool->OverlapRemoval(electrons, muons, jets, photons));
        
        xAOD::MissingETContainer*    met_syst = new xAOD::MissingETContainer;
        xAOD::MissingETAuxContainer* met_syst_aux = new xAOD::MissingETAuxContainer;
        met_syst->setStore(met_syst_aux);


        // CHECK(objTool->GetMET(*met_syst, jets, electrons, muons));
        // CHECK(objTool->GetMET(*met_syst, jets_nominal, electrons_nominal, muons_nominal, photons_nominal));
        CHECK(objTool->GetMET(*met_syst, jets, electrons, muons, photons));

        met = met_syst;
        met_aux = met_syst_aux;
        
        for (const auto& el : *electrons) {
          if (!isData && el->auxdata<char>("baseline") == 1 &&
              el->auxdata<char>("passOR") == 1 &&
              el->auxdata<char>("signal") == 1)
            weight_sf *= objTool->GetSignalElecSF(*el);
        }
        
        bool skip = false;
        for (const auto& mu : *muons) {
          if (!isData && mu->auxdata<char>("baseline") == 1 &&
              mu->auxdata<char>("passOR") == 1 && 
              mu->auxdata<char>("signal") == 1)
            objTool->GetSignalMuonSF(*mu);
          weight_sf *= objTool->GetSignalMuonSF(*mu);
        }
        
        for (const auto& jet : *jets) {
          if ((int)jet->auxdata<char>("bad")) {
            skip = true;
            break;
          }
        }
        
        float weight_btag = objTool->BtagSF(jets);
        
        if (!skip) {
          AnalysisCollections collections;
          collections.event_number = event_number;
          collections.avg_mu = avg_mu;

          collections.weight_pu = weight_pu;
          collections.weight_mc = weight_mc;

          //collections.SFweight = weight_sf;	
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
        if (syst_affectsJets || syst_affectsBTag ) {
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
  //float weight_sf = 1.;
  //float weight_btag = 1.;
  

  CHECK(objTool->OverlapRemoval(electrons_nominal, muons_nominal, jets_nominal, photons_nominal));

  // MET
  CHECK(objTool->GetMET(*met_nominal, jets_nominal, electrons_nominal, muons_nominal, photons_nominal));
  
  // electrons
  for (const auto& el : *electrons_nominal) {
    if (!isData && el->auxdata<char>("baseline") == 1 &&
        el->auxdata<char>("passOR") == 1 &&
        el->auxdata<char>("signal") == 1)
        objTool->GetSignalElecSF(*el);
  }
  
  // // photons
  // for (const auto& ph : *photons_nominal) {
  //   //if (!isData && ph->auxdata<char>("baseline") == 1 &&
  //   //       ph->auxdata<char>("passOR") == 1 && 
  //   //       ph->auxdata<char>("signal") == 1)
  //   //     objTool->GetSignalPhotonSF(*ph);
  // }
  
  bool skip = false;
  // muons
  for (const auto& mu : *muons_nominal) {
    if (!isData && mu->auxdata<char>("baseline") == 1 &&
        mu->auxdata<char>("passOR") == 1 && 
        mu->auxdata<char>("signal") == 1)
      objTool->GetSignalMuonSF(*mu);
  }
  if (!skip)
    h_cutflow->Fill(6);

  // jets
  for (const auto& jet : *jets_nominal) {  
    if ((int)jet->auxdata<char>("bad")) {
      skip = true;
      break;
    }
  }
  
  if (!skip) {
    h_cutflow->Fill(7);

    AnalysisCollections collections;
    collections.event_number = event_number;
    collections.avg_mu = avg_mu;

    collections.weight_pu = weight_pu;
    collections.weight_mc = weight_mc;

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
  
  if (objTool) {
    delete objTool;
    objTool = 0;
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


// static SG::AuxElement::Decorator<char> dec_baseline("baseline");
// static SG::AuxElement::Decorator<char> dec_signal("signal");
// static SG::AuxElement::Decorator<char> dec_isol("isol");
// static SG::AuxElement::Decorator<char> dec_passOR("passOR");                                                 
// static SG::AuxElement::Decorator<char> dec_passSignalID("passSignalID");
// static SG::AuxElement::Decorator<char> dec_bad("bad");
// static SG::AuxElement::Decorator<char> dec_bjet("bjet");
// static SG::AuxElement::Decorator<char> dec_bjet_loose("bjet_loose");
// static SG::AuxElement::Decorator<char> dec_cosmic("cosmic");
// static SG::AuxElement::Decorator<char> dec_passedHighPtCuts("passedHighPtCuts");


// EL::StatusCode xAODAnalysis::OverlapRemoval(const xAOD::ElectronContainer *electrons, const xAOD::MuonContainer *muons, const xAOD::JetContainer *jets, const xAOD::PhotonContainer *photons, const bool useSignalObjects, const bool useIsolObjects, const bool doBjetOR, const bool doBoostedMuonOR, const double dRejet, const double dRjetmu, const double dRjete, double dRemu, double dRee, double dRphjet, double dReph, double dRmuph) {

//   for (const auto& jet : *jets) {
//     dec_passOR( *jet ) = dec_baseline( *jet );
//   }
//   for (const auto& mu : *muons) {
//     bool mu_sel;
//     if (useSignalObjects) mu_sel = dec_signal(*mu);
//     else mu_sel = dec_baseline(*mu);
//     if (useIsolObjects) mu_sel &= dec_isol(*mu);
//     dec_passOR( *mu ) = mu_sel;
//   }

//   // remove jets overlapping with (baseline/signal) electrons
//   for (const auto& el : *electrons) {
//     bool el_sel;
//     if (useSignalObjects) el_sel = dec_signal(*el);
//     else el_sel = dec_baseline(*el);
//     if (useIsolObjects) el_sel &= dec_isol(*el);
//     dec_passOR( *el ) = el_sel;

//     if ( !el_sel ) continue;

//     for (const auto& jet : *jets) {
//       if ( !dec_passOR(*jet) ) continue;
//       if ( doBjetOR && dec_bjet_loose(*jet) ) continue;

//       if ( xAOD::P4Helpers::deltaR2(el, jet) < dRejet * dRejet ) {
//         ATH_MSG_VERBOSE( " Rejecting jet at (eta,phi)=(" << jet->eta() << "," << jet->phi() << ") "
//                          << " due to electron at (eta,phi)=(" << el->eta() << "," << el->phi() << ")" );
//         dec_passOR( *jet ) = false;
//       }
//     }
//   } // END loop over electrons

//   // remove jets overlapping with (baseline/signal) photons
//   for (const auto& ph : *photons) {
//     bool ph_sel;
//     if (useSignalObjects) ph_sel = dec_signal(*ph);
//     else ph_sel = dec_baseline(*ph);
//     if (useIsolObjects) ph_sel &= dec_isol(*ph);
//       dec_passOR( *ph ) = ph_sel;

//     if ( !ph_sel ) continue;

//     for (const auto& jet : *jets) {
//       if ( !dec_passOR(*jet) ) continue;

//       if ( xAOD::P4Helpers::deltaR2(ph, jet) < dRphjet * dRphjet ) {
//         ATH_MSG_VERBOSE( " Rejecting jet at (eta,phi)=(" << jet->eta() << "," << jet->phi() << ") "
//                          << " due to photon at (eta,phi)=(" << ph->eta() << "," << ph->phi() << ")" );
//         dec_passOR( *jet ) = false;
//       }
//     }
//   } // END loop over photons

//   // Remove electrons and muons overlapping with jets
//   for (const auto& el : *electrons) {
//     if ( !dec_passOR(*el) ) continue;

//     for (const auto& jet : *jets) {
//       if ( !dec_passOR( *jet ) ) continue;

//       if ( xAOD::P4Helpers::deltaR2(el, jet) < dRjete * dRjete ) {
//         ATH_MSG_VERBOSE( " Rejecting electron at (eta,phi)=(" << el->eta() << "," << el->phi() << ") "
//                          << " due to jet at (eta,phi)=(" << jet->eta() << "," << jet->phi() << ")" );
//         dec_passOR( *el ) = false;
//       }
//     }              
//   }

//   for (const auto& mu : *muons) {
//     if ( !dec_passOR(*mu) ) continue;

//     for (const auto& jet : *jets) {
//       if ( !dec_passOR( *jet ) ) continue;

//       std::vector<int> nTrkVec;
//       jet->getAttribute(xAOD::JetAttribute::NumTrkPt500, nTrkVec);
//       int jet_nTrk = (objTool->GetPrimVtx() == 0 || nTrkVec.size() == 0) ? 0 : nTrkVec[objTool->GetPrimVtx()->index()];

//       bool muInJet = xAOD::P4Helpers::deltaR2(*mu, *jet) < dRjetmu * dRjetmu;
//       if(doBoostedMuonOR) muInJet = xAOD::P4Helpers::deltaR2(*mu, *jet) < 0.04+10e3/mu->pt();
//       if ( muInJet ) {
//         if (jet_nTrk < 3) {
//           ATH_MSG_VERBOSE( " Rejecting jet at (pT,eta,phi)=(" << jet->pt() << "," << jet->eta() << "," << jet->phi() << ") with only nTrk=" << jet_nTrk
//                            << " due to muon at (pT,eta,phi)=(" << mu->pt() << "," << mu->eta() << "," << mu->phi() << ")" );
//           dec_passOR( *jet ) = false;
//         } else {
//           ATH_MSG_VERBOSE( " Rejecting muon at (pT,eta,phi)=(" << mu->pt() << "," << mu->eta() << "," << mu->phi() << ")"
//                            << " due to jet at (pT,eta,phi)=(" << jet->pt() << "," << jet->eta() << "," << jet->phi() << ") with nTrk=" << jet_nTrk );
//           dec_passOR( *mu ) = false;
//         }
//       }
//     }        
//   }

//   // Remove electrons and muons overlapping with each other
//   for (const auto& el : *electrons) {
//     if ( !dec_passOR(*el) ) continue;

//     for (const auto& mu : *muons) {
//       if ( !dec_passOR( *mu ) ) continue;

//       if ( xAOD::P4Helpers::deltaR2(el, mu) < dRemu * dRemu ) {
//         ATH_MSG_VERBOSE( " Rejecting both electron at (eta,phi)=(" << el->eta() << "," << el->phi() << ") "
//                          << " and muon at (eta,phi)=(" << mu->eta() << "," << mu->phi() << ")" );
//         dec_passOR( *el ) = false;
//         dec_passOR( *mu ) = false;
//       }
//     }
//   }

//   // Remove electrons overlapping with each other
//   for (auto el_itr = electrons->begin(); el_itr != electrons->end(); ++el_itr) {
//     const auto& el = *el_itr;
//     if ( !dec_passOR(*el) ) continue;

//     for (auto el2_itr = el_itr + 1; el2_itr != electrons->end(); ++el2_itr) {
//       if (el2_itr == el_itr) continue;
//       const auto& el2 = *el2_itr;
//       if ( !dec_passOR( *el2 ) ) continue;

//       if ( xAOD::P4Helpers::deltaR2(el, el2) < dRee * dRee ) {
//         if (el->pt() < el2->pt()) {
//           ATH_MSG_VERBOSE( " Rejecting electron at (eta,phi)=(" << el->eta() << "," << el->phi() << ") "
//                            << " and muon at (eta,phi)=(" << el2->eta() << "," << el2->phi() << ")" );
//           dec_passOR( *el ) = false;
//         } else {
//           ATH_MSG_VERBOSE( " Rejecting electron at (eta,phi)=(" << el2->eta() << "," << el2->phi() << ") "
//                            << " and muon at (eta,phi)=(" << el->eta() << "," << el->phi() << ")" );
//           dec_passOR( *el2 ) = false;
//         }
//       }
//     }
//   }

//   // Remove photons overlapping with electrons
//   for (const auto& ph : *photons) {
//     if ( !dec_passOR( *ph ) )
//       continue;

//     for (const auto& el : *electrons) {
//       if ( !dec_passOR(*el) ) continue;

//       if ( xAOD::P4Helpers::deltaR2(el, ph) < dReph * dReph ) {
//         ATH_MSG_VERBOSE( " Rejecting photon at (eta,phi)=(" << ph->eta() << "," << ph->phi() << ") "
//                          << " due to electron at (eta,phi)=(" << el->eta() << "," << el->phi() << ")" );
//         dec_passOR( *ph ) = 0;
//       }
//     }
//   }

//   // Remove photons overlapping with muons
//   for (const auto& ph : *photons) {
//     if ( !dec_passOR( *ph ) )
//       continue;

//     for (const auto& mu : *muons) {
//       if ( !dec_passOR(*mu) ) continue;

//       if ( xAOD::P4Helpers::deltaR2(mu, ph) < dRmuph * dRmuph ) {
//         ATH_MSG_VERBOSE( " Rejecting photon at (eta,phi)=(" << ph->eta() << "," << ph->phi() << ") "
//                          << " due to muon at (eta,phi)=(" << mu->eta() << "," << mu->phi() << ")" );
//         dec_passOR( *ph ) = 0;
//       }
//     }
//   }

//   // Count number of objects after overlap removal
//   int Nel = 0;
//   for (const auto& el : *electrons) {
//     if (dec_passOR( *el )) Nel++;
//   }

//   int Nmu = 0;
//   for (const auto& mu : *muons) {
//     if (dec_passOR( *mu )) Nmu++;
//   }

//   int Njet = 0;
//   for (const auto& jet : *jets) {
//     if (dec_passOR( *jet )) Njet++;
//   }

//   int Nph = 0;
//   for (const auto& ph : *photons) {
//     if (dec_passOR( *ph )) Nph++;
//   }

//   ATH_MSG_VERBOSE( " After overlap removal: Nel=" << Nel << ", Nmu=" << Nmu << ", Njet=" << Njet << ", Nph=" << Nph);
//   return StatusCode::SUCCESS;
// }
