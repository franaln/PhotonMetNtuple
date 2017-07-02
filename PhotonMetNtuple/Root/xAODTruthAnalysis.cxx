#include <TSystem.h>
#include <iostream>

// ROOT include(s):
#include <TChain.h>
#include <TError.h>
#include <TH1D.h>
#include <TFile.h>

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <PhotonMetNtuple/xAODTruthAnalysis.h>
#include "EventLoop/OutputStream.h"
#include <TTreeFormula.h>

#include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"

// EDM includes:
#include "xAODMissingET/MissingETContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

#include "EventPrimitives/EventPrimitivesHelpers.h"

#include "PhotonMetNtuple/Common.h"
#include "PhotonMetNtuple/Utils.h"
#include "PhotonMetNtuple/TruthUtils.h"
#include "PhotonMetNtuple/TruthTree.h"

// pdf reweighting
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Reweighting.h"


ClassImp(xAODTruthAnalysis)

struct pt_sort {
  bool operator()(const TruthParticle& v1, const TruthParticle& v2) const {
    return v1.Pt() > v2.Pt();
  }
};

xAODTruthAnalysis::xAODTruthAnalysis() :
  do_pdfrw(false),
  do_lhe3(false),
  is_truth3(false)
{
}

EL::StatusCode xAODTruthAnalysis::setupJob(EL::Job &job)
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
  xAOD::Init("xAODTruthAnalysis").ignore(); // call before opening first file
  
  // tell EventLoop about our output:
  EL::OutputStream out("output");
  job.outputAdd(out);
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODTruthAnalysis::histInitialize()
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

  TDirectory *out_dir = (TDirectory*) wk()->getOutputFile("output");
  h_events->SetDirectory(out_dir);

  if (do_pdfrw) {
    pdf1 = "CT10";
    pdf2 = "NNPDF30_lo_as_0130"; //# 
    pdf3 = "MMHT2014lo68cl"; //cteq66"; // MMHT2014

    Info(APP_NAME, "Loading PDF set %s", pdf1.c_str());
    const LHAPDF::PDFSet pdf_set_1(pdf1);
    m_pdfs_1 = pdf_set_1.mkPDFs();
 
    Info(APP_NAME, "Loading PDF set %s", pdf2.c_str());
    const LHAPDF::PDFSet pdf_set_2(pdf2);
    m_pdfs_2 = pdf_set_2.mkPDFs();
  
    Info(APP_NAME, "Loading PDF set %s", pdf3.c_str());
    const LHAPDF::PDFSet pdf_set_3(pdf3);
    m_pdfs_3 = pdf_set_3.mkPDFs();
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODTruthAnalysis::fileExecute()
{
  // Here you do everything that n<eeds to be done exactly once for every
  // single file, e.g. collect a st of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODTruthAnalysis::changeInput(bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  
  m_event = wk()->xaodEvent(); 

  TTree *metadata = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
  if (!metadata) {
    Error(APP_NAME, "MetaData tree not found! Exiting.");
    return EL::StatusCode::FAILURE;
  }
  metadata->LoadTree(0);

  //check if file is from a DxAOD
  bool is_derivation = !metadata->GetBranch("StreamAOD");

  std::cout << "is derivation" << is_derivation << std::endl;

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
    
    int maxCycle = -1;
    for (const auto& cbk :  *completeCBC) {
      
      // std::cout << cbk->nameIdentifier() << " : " << cbk->name() << " : desc = " << cbk->description()
      //           << " : inputStream = " << cbk->inputStream()  << " : outputStreams = " << (cbk->outputStreams().size() ? cbk->outputStreams()[0] : "")
      //           << " : cycle = " << cbk->cycle() << " :  allEvents = " << cbk->nAcceptedEvents()
      //           << std::endl;
      
      if (cbk->cycle() > maxCycle && cbk->name() == "AllExecutedEvents" && (cbk->inputStream() == "StreamAOD" || cbk->inputStream() == "StreamDAOD_TRUTH3" || cbk->inputStream() == "StreamDAOD_TRUTH1")) {
        all_events_cbk = cbk;
        maxCycle = cbk->cycle();
        
        if (cbk->inputStream() == "StreamDAOD_TRUTH3")
          is_truth3 = true;
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

EL::StatusCode xAODTruthAnalysis::initialize()
{
  // create "mini" ntuple
  TDirectory *out_dir = (TDirectory*) wk()->getOutputFile("output");
  
  ntuple = new TruthTree("truthtree");
  CHECK(ntuple->setProperty("OutFile", out_dir));
  CHECK(ntuple->setProperty("SavePDF", do_pdfrw));
  CHECK(ntuple->setProperty("SaveLHE3", do_lhe3));
  CHECK(ntuple->initialize());

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODTruthAnalysis::execute ()
{
  
  // const xAOD::TruthParticleContainer* truth_particles = 0;
  // if(!m_event->retrieve(truth_particles, "TruthParticles").isSuccess()) {
  //   Error(APP_NAME, "Failed to retrieve truth particles collection. Exiting." );
  //   return EL::StatusCode::FAILURE;
  // }

  const xAOD::JetContainer *truth_jets = 0;
  if(!m_event->retrieve(truth_jets, "AntiKt4TruthJets").isSuccess()) {
    Error(APP_NAME, "Failed to retrieve truth jets collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::TruthParticleContainer *truth_photons = 0;
  if (is_truth3) {
    if (!m_event->retrieve(truth_photons, "Truth3Photons").isSuccess()) {
      Error(APP_NAME, "Failed to retrieve truth photons collection. Exiting." );
      return EL::StatusCode::FAILURE;
    }
  }
  else {
    if (!m_event->retrieve(truth_photons, "TruthPhotons").isSuccess()) {
      Error(APP_NAME, "Failed to retrieve truth photons collection. Exiting." );
      return EL::StatusCode::FAILURE;
    }
  }

  const xAOD::TruthParticleContainer *truth_electrons = 0;
  if(!m_event->retrieve(truth_electrons, "TruthElectrons").isSuccess()) {
    Error(APP_NAME, "Failed to retrieve truth electrons collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::TruthParticleContainer *truth_muons = 0;
  if(!m_event->retrieve(truth_muons, "TruthMuons").isSuccess()) {
    Error(APP_NAME, "Failed to retrieve truth photons collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::MissingETContainer *truth_met = 0;
  if(!m_event->retrieve(truth_met, "MET_Truth").isSuccess()) {
    Error(APP_NAME, "Failed to retrieve truth met collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::MissingET* met_truth = 0;
  if( truth_met ) met_truth = (*truth_met)["NonInt"];

  ntuple->Clear();

  std::vector<TruthParticle> b_photons;
  std::vector<TruthParticle> b_electrons;
  std::vector<TruthParticle> b_muons;
  std::vector<TruthParticle> b_jets;
  
  std::vector<TruthParticle> photons;
  std::vector<TruthParticle> electrons;
  std::vector<TruthParticle> muons;
  std::vector<TruthParticle> jets;

  b_photons.clear();
  b_electrons.clear();
  b_muons.clear();
  b_jets.clear();

  photons.clear();
  electrons.clear();
  muons.clear();
  jets.clear();
  
  //-- Baseline selection
  // Photons
  xAOD::TruthParticleContainer::const_iterator ph_itr(truth_photons->begin());
  xAOD::TruthParticleContainer::const_iterator ph_end(truth_photons->end());
  for ( ; ph_itr != ph_end; ++ph_itr) {
      
    if ((*ph_itr)->pt() < 25000. || fabs((*ph_itr)->eta()) > 2.37) 
      continue;

    b_photons.push_back(TruthParticle(*ph_itr));
  }

  // Electrons
  xAOD::TruthParticleContainer::const_iterator el_itr(truth_electrons->begin());
  xAOD::TruthParticleContainer::const_iterator el_end(truth_electrons->end());
  for ( ; el_itr != el_end; ++el_itr) {
      
    if ((*el_itr)->pt() < 10000. || fabs((*el_itr)->eta()) > 2.47) 
        continue;

      b_electrons.push_back(TruthParticle(*el_itr));
  }

  // Muons
  xAOD::TruthParticleContainer::const_iterator mu_itr(truth_muons->begin());
  xAOD::TruthParticleContainer::const_iterator mu_end(truth_muons->end());
  for ( ; mu_itr != mu_end; ++mu_itr) {
      
    if ((*mu_itr)->pt() < 10000. || fabs((*mu_itr)->eta()) > 2.7) 
        continue;

      b_muons.push_back(TruthParticle(*mu_itr));
  }

  // Jets
  xAOD::JetContainer::const_iterator jet_itr(truth_jets->begin());
  xAOD::JetContainer::const_iterator jet_end(truth_jets->end());
  for ( ; jet_itr != jet_end; ++jet_itr) {

    if ((*jet_itr)->pt() < 20000. || fabs((*jet_itr)->eta()) > 2.8)
      continue;
        
    int flavor = 0;
    if ((*jet_itr)->isAvailable<int>("PartonTruthLabelID"))
      flavor = abs((*jet_itr)->auxdata<int>("PartonTruthLabelID"));

    bool isb = (flavor == 5); //(*jet_itr)->auxdata<int>("GhostBHadronsFinalCount") >= 1;

    b_jets.push_back(TruthParticle(*jet_itr, isb));
  }

  std::sort(b_photons.begin(), b_photons.end(), pt_sort());
  std::sort(b_electrons.begin(), b_electrons.end(), pt_sort());
  std::sort(b_muons.begin(), b_muons.end(), pt_sort());
  std::sort(b_jets.begin(), b_jets.end(), pt_sort());
  
  // itearatos
  std::vector<TruthParticle>::iterator ph_it;
  std::vector<TruthParticle>::iterator el_it;
  std::vector<TruthParticle>::iterator mu_it;
  std::vector<TruthParticle>::iterator jet_it;
  
  //-- Overlap Removal
  Double_t ph_el_deltaR_cut    = 0.01;
  Double_t eg_jet_deltaR_cut   = 0.2;
  Double_t jet_egmu_deltaR_cut = 0.4;
  
  // Remove photons overlapping with electrons
  for (ph_it=b_photons.begin(); ph_it!=b_photons.end(); ++ph_it) {
    if (TruthUtils::OverlapsOthers((*ph_it), b_electrons, ph_el_deltaR_cut) ) (*ph_it).good=false;
  }
  TruthUtils::CleanBads(b_photons);
  
  /// Remove jets overlapping with electrons and photons
  for (jet_it=b_jets.begin(); jet_it!=b_jets.end(); ++jet_it) {
    if (TruthUtils::OverlapsOthers((*jet_it), b_electrons, eg_jet_deltaR_cut)) (*jet_it).good=false;
  }
  TruthUtils::CleanBads(b_jets);
  for (jet_it=b_jets.begin(); jet_it!=b_jets.end(); ++jet_it) {
    if (TruthUtils::OverlapsOthers((*jet_it), b_photons,   eg_jet_deltaR_cut)) (*jet_it).good=false;  
  }
  TruthUtils::CleanBads(b_jets);
  
  /// Remove electrons and photons comming from jets
  for (el_it=b_electrons.begin(); el_it!=b_electrons.end(); ++el_it) {
      if (TruthUtils::OverlapsOthers((*el_it), b_jets, jet_egmu_deltaR_cut)) (*el_it).good = false;
  }
  TruthUtils::CleanBads(b_electrons);
  for (ph_it=b_photons.begin(); ph_it!=b_photons.end(); ++ph_it) {
    if (TruthUtils::OverlapsOthers((*ph_it), b_jets, jet_egmu_deltaR_cut)) (*ph_it).good = false;
  }
  TruthUtils::CleanBads(b_photons);
  for(mu_it=b_muons.begin(); mu_it!=b_muons.end(); ++mu_it) {
    if (TruthUtils::OverlapsOthers((*mu_it), b_jets, jet_egmu_deltaR_cut)) (*mu_it).good=false;
  }
  TruthUtils::CleanBads(b_muons);


  // Signal objects
  ntuple->ph_n = 0;
  ntuple->ph_pt->clear();
  ntuple->ph_eta->clear();
  ntuple->ph_phi->clear();
  ntuple->ph_iso->clear();
  for (ph_it=b_photons.begin(); ph_it!=b_photons.end(); ++ph_it) {

    Float_t ph_eta = fabs((*ph_it).Eta());

    if ((*ph_it).Pt() < 75. ||  (ph_eta>1.37 && ph_eta<1.52))
      continue;

    photons.push_back(*ph_it);
      
    ntuple->ph_n++;
      
    ntuple->ph_pt->push_back((*ph_it).Pt());
    ntuple->ph_eta->push_back((*ph_it).Eta());
    ntuple->ph_phi->push_back((*ph_it).Phi());
    ntuple->ph_iso->push_back(0.);
  }

  ntuple->el_n = 0;
  ntuple->el_pt->clear();
  ntuple->el_eta->clear();
  ntuple->el_phi->clear();
  for (el_it=b_electrons.begin(); el_it!=b_electrons.end(); ++el_it) {

    Float_t el_eta = fabs((*el_it).Eta());

    if ((*el_it).Pt() < 25. || (el_eta>1.37 && el_eta<1.52))
        continue;

    electrons.push_back(*el_it);
    
    ntuple->el_n++;
    ntuple->el_pt->push_back((*el_it).Pt());
    ntuple->el_eta->push_back((*el_it).Eta());
    ntuple->el_phi->push_back((*el_it).Phi());
  }

  ntuple->mu_n = 0;
  ntuple->mu_pt->clear();
  ntuple->mu_eta->clear();
  ntuple->mu_phi->clear();
  for (mu_it=b_muons.begin(); mu_it!=b_muons.end(); ++mu_it) {

    if ((*mu_it).Pt() < 25.)
        continue;

    muons.push_back(*mu_it);
    
    ntuple->mu_n++;
    ntuple->mu_pt->push_back((*mu_it).Pt());
    ntuple->mu_eta->push_back((*mu_it).Eta());
    ntuple->mu_phi->push_back((*mu_it).Phi());
  }

  ntuple->jet_n = 0;
  ntuple->bjet_n = 0;
  ntuple->jet_pt->clear();
  ntuple->jet_eta->clear();
  ntuple->jet_phi->clear();
  for (jet_it=b_jets.begin(); jet_it!=b_jets.end(); ++jet_it) {

    if ((*jet_it).Pt() < 30. || fabs((*jet_it).Eta()) > 2.8)
        continue;

    jets.push_back(*jet_it);

    ntuple->jet_n++;
    ntuple->jet_pt->push_back((*jet_it).Pt());
    ntuple->jet_eta->push_back((*jet_it).Eta());
    ntuple->jet_phi->push_back((*jet_it).Phi());

    if ((*jet_it).isbjet)
      ntuple->bjet_n++;
  }

  //-- Construct event variables
  
  /// Met
  Double_t ex = 0.0;
  Double_t ey = 0.0;
  
  // Met truth
  ex = met_truth->mpx() * 0.001;
  ey = met_truth->mpy() * 0.001;
  
  Double_t met_et = TMath::Hypot(ex, ey);
  Double_t met_phi = TMath::ATan2(ey, ex);

  ntuple->met_et = met_et;
  ntuple->met_phi = met_phi;

  /// Ht
  Double_t sum_jet_pt = 0.;
  for (jet_it=jets.begin(); jet_it!=jets.end(); ++jet_it)
    sum_jet_pt += (*jet_it).Pt();
  
  Double_t ht = sum_jet_pt;
  if (photons.size() > 0)
    ht += photons[0].Pt();
  
  ntuple->ht = ht;
  ntuple->meff = ht + ntuple->met_et;
  
  /// Rt
  if (jets.size() >= 2)
    ntuple->rt2 = (jets[0].Pt() + jets[1].Pt())/sum_jet_pt;
  
  if (jets.size() >= 4) 
    ntuple->rt4 = (jets[0].Pt() + jets[1].Pt() + jets[2].Pt() + jets[3].Pt())/sum_jet_pt;
  
  /// dphi's
  // Delta phi between met and closest jet
  Double_t dphi1 = 999.;
  Double_t dphi2 = 999.;
  
  if (jets.size() > 0) dphi1 = get_dphi(jets[0].Phi(), met_phi);
  if (jets.size() > 1) dphi2 = get_dphi(jets[1].Phi(), met_phi);
  ntuple->dphi_jetmet = TMath::Min(dphi1, dphi2);
    
  // dphi beteen leading photon and leading jet
  if (photons.size() > 0 && jets.size() > 0) ntuple->dphi_gamjet = get_dphi(photons[0].Phi(), jets[0].Phi());
  
  // dphi beteen leading photon and MET
  if (photons.size() > 0) ntuple->dphi_gammet = get_dphi(photons[0].Phi(), met_phi);


  // pdf
  if (do_pdfrw) { 
    const xAOD::TruthEventContainer* truthEvent = 0;
    if( m_event->retrieve( truthEvent, "TruthEvents" ).isFailure()) {
      Error("Execute()", "TruthEvent not found! pdf calculation will be wrong");
    } 
    
    xAOD::TruthEvent::PdfInfo pdata = truthEvent->at(0)->pdfInfo();
    
    for (size_t iPdf=0; iPdf<m_pdfs_1.size(); iPdf++) {
      ntuple->weight_pdf1->push_back( LHAPDF::weightxxQ(pdata.pdgId1, pdata.pdgId2, pdata.x1, pdata.x2, pdata.Q, m_pdfs_1[0], m_pdfs_1[iPdf]) );
    }
    
    for (size_t iPdf=0; iPdf<m_pdfs_2.size(); iPdf++) {
      ntuple->weight_pdf2->push_back( LHAPDF::weightxxQ(pdata.pdgId1, pdata.pdgId2, pdata.x1, pdata.x2, pdata.Q, m_pdfs_2[0], m_pdfs_2[iPdf]) );
    }
    
    for (size_t iPdf=0; iPdf<m_pdfs_3.size(); iPdf++) {
      ntuple->weight_pdf3->push_back( LHAPDF::weightxxQ(pdata.pdgId1, pdata.pdgId2, pdata.x1, pdata.x2, pdata.Q, m_pdfs_3[0], m_pdfs_3[iPdf]) );
    }
  }

  //Get LHE3 weights
  const xAOD::TruthEventContainer* truth_evt_cont;
  if ( !m_event->retrieve(truth_evt_cont, "TruthEvents").isSuccess() ) {
    Error(APP_NAME, "Could not retrieve truth event container with key TruthEvents");
  }
  const xAOD::TruthEvent *truth_event = (*truth_evt_cont)[0];
  const std::vector<float> weights  = truth_event->weights();
  for (auto w : weights)
    ntuple->weight_lhe3->push_back(w);
 

  // fill ntuple only there is at least one photon
  // if (photons.size() > 0)
  ntuple->Fill();

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODTruthAnalysis::postExecute()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODTruthAnalysis::finalize()
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode xAODTruthAnalysis::histFinalize()
{
  return EL::StatusCode::SUCCESS;
}


