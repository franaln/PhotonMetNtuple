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

xAODTruthAnalysis::xAODTruthAnalysis() :
do_pdfrw(false)
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
    pdf2 = "NNPDF30_lo_as_0130";
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

  TTree* CollectionTree = dynamic_cast<TTree*>( wk()->inputFile()->Get("CollectionTree") );
  

  ULong64_t m_initial_events = 0;
  Double_t m_initial_sumw = 0.;
  Double_t m_initial_sumw2 = 0.;
  
  m_initial_events  = CollectionTree->GetEntries(); 
  m_initial_sumw    = CollectionTree->GetWeight() * CollectionTree->GetEntries();
  m_initial_sumw2   = (CollectionTree->GetWeight() * CollectionTree->GetWeight()) * CollectionTree->GetEntries();
  
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
  CHECK(ntuple->initialize());

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode xAODTruthAnalysis::execute ()
{
  
  const xAOD::TruthParticleContainer* truth_particles = 0;
  if(!m_event->retrieve(truth_particles, "TruthParticles").isSuccess()) {
    Error(APP_NAME, "Failed to retrieve truth particles collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  const xAOD::JetContainer *truth_jets = 0;
  if(!m_event->retrieve(truth_jets, "AntiKt4TruthJets").isSuccess()) {
    Error(APP_NAME, "Failed to retrieve truth jets collection. Exiting." );
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
  
  //-- Baseline selection
  for (auto truth = truth_particles->begin(); truth!=truth_particles->end(); ++truth) {

    if (!is_stable(*truth))
      continue;
    
    Int_t pid = (*truth)->pdgId();
    Int_t abspid = abs((*truth)->pdgId());
    
    Int_t st = (*truth)->status();
    
    Double_t pt = (*truth)->pt() * 0.001;
    Double_t eta = (*truth)->eta();
    Double_t phi = (*truth)->phi();
    
    Double_t abseta = fabs(eta);
    
    // Photons
    if (pid == 22 && st == 1) {
      
      if (pt < 25. || abseta > 2.37) 
        continue;
      
      b_photons.push_back(TruthParticle(*truth));
      
      if (pt < 75.)
        continue;
      
      photons.push_back(TruthParticle(*truth));
      
      ntuple->ph_n++;
      
      ntuple->ph_pt->push_back(pt);
      ntuple->ph_eta->push_back(eta);
      ntuple->ph_phi->push_back(phi);
      ntuple->ph_iso->push_back(0.);
      ntuple->ph_type->push_back(get_photon_type(*truth));
    }
    
    // Electrons
    else if (abspid == 11 && st == 1) {
      
      if (pt < 10 || abseta > 2.47)
        continue;
      
      b_electrons.push_back(TruthParticle(*truth));
      
      if (pt < 25)
        continue;
      
      electrons.push_back(TruthParticle(*truth));
      
      ntuple->el_n++;
      ntuple->el_pt->push_back(pt);
      ntuple->el_eta->push_back(eta);
      ntuple->el_phi->push_back(phi);
      
    }
    
    // Muons
    else if (abspid == 13 && st == 1) {
      
      if (pt < 10 || abseta > 2.4)
        continue;
      
      b_muons.push_back(TruthParticle(*truth));        
      
      if (pt < 25)
        continue;
      
      muons.push_back(TruthParticle(*truth));        
      
      ntuple->mu_n++;
      ntuple->mu_pt->push_back(pt);
      ntuple->mu_eta->push_back(eta);
      ntuple->mu_phi->push_back(phi);
      
    }
    
  } // end truth particle container loop 


  // Jets
  xAOD::JetContainer::const_iterator jet_itr(truth_jets->begin());
  xAOD::JetContainer::const_iterator jet_end(truth_jets->end());

  for ( ; jet_itr != jet_end; ++jet_itr) {
    
    Double_t pt = (*jet_itr)->pt() * 0.001;
    Double_t eta = (*jet_itr)->eta();
    Double_t phi = (*jet_itr)->phi();
    Double_t abseta = fabs(eta);
    
    if (pt < 20 || abseta > 2.8)
      continue;
    
    b_jets.push_back(TruthParticle(*jet_itr));
    
    if (pt < 30 || abseta > 2.5)
      continue;
        
    jets.push_back(TruthParticle(*jet_itr));
    
    ntuple->jet_n++;
    ntuple->jet_pt->push_back(pt);
    ntuple->jet_eta->push_back(eta);
    ntuple->jet_phi->push_back(phi);
    
  }
      
  // we dont want events without photons
  // if (photons.size() == 0)
  //   continue;
  
  // itearatos
  std::vector<TruthParticle>::iterator ph_it;
  std::vector<TruthParticle>::iterator el_it;
  std::vector<TruthParticle>::iterator mu_it;
  std::vector<TruthParticle>::iterator jet_it;
  
 
  //-- Overlap Removal
  Double_t ph_el_deltaR_cut    = 0.01;
  Double_t ph_jet_deltaR_cut   = 0.2;
  Double_t el_jet_deltaR_cut   = 0.2;
  Double_t jet_egmu_deltaR_cut = 0.4;
  
  // Remove photons overlapping with electrons
  for (ph_it=b_photons.begin(); ph_it!=b_photons.end(); ++ph_it) {
    if (TruthUtils::OverlapsOthers((*ph_it), b_electrons, ph_el_deltaR_cut) ) (*ph_it).good=false;
  }
  TruthUtils::CleanBads(b_photons);
  
  /// Remove jets overlapping with electrons and photons
  for (jet_it=b_jets.begin(); jet_it!=b_jets.end(); ++jet_it) {
    if (TruthUtils::OverlapsOthers((*jet_it), b_electrons, el_jet_deltaR_cut)) (*jet_it).good=false;
    if (TruthUtils::OverlapsOthers((*jet_it), b_photons, ph_jet_deltaR_cut)) (*jet_it).good=false;  
  }
  TruthUtils::CleanBads(b_jets);
  
  /// Remove electrons and photons comming from jets
  for (el_it=b_electrons.begin(); el_it!=b_electrons.end(); ++el_it) {
      if (TruthUtils::OverlapsOthers((*el_it), b_jets, jet_egmu_deltaR_cut)) (*el_it).good = false;
  }
  for (ph_it=b_photons.begin(); ph_it!=b_photons.end(); ++ph_it) {
    if (TruthUtils::OverlapsOthers((*ph_it), b_jets, jet_egmu_deltaR_cut)) (*ph_it).good = false;
  }
  TruthUtils::CleanBads(b_electrons);
  TruthUtils::CleanBads(b_photons);
  
  /// Remove muons overlapping with remaining jets
  for(mu_it=b_muons.begin(); mu_it!=b_muons.end(); ++mu_it) {
    if (TruthUtils::OverlapsOthers((*mu_it), b_jets, jet_egmu_deltaR_cut)) (*mu_it).good=false;
  }
  TruthUtils::CleanBads(b_muons);
  

  //-- Construct event variables
  
  /// Met
  Double_t ex = 0.0;
  Double_t ey = 0.0;
  
  for (ph_it=b_photons.begin(); ph_it!=b_photons.end(); ++ph_it) {
    ex -= (*ph_it).Et() * TMath::Cos((*ph_it).Phi());
    ey -= (*ph_it).Et() * TMath::Sin((*ph_it).Phi());
  }
  for (el_it=b_electrons.begin(); el_it!=b_electrons.end(); ++el_it) {
    ex -= (*el_it).Et() * TMath::Cos((*el_it).Phi());
    ey -= (*el_it).Et() * TMath::Sin((*el_it).Phi());
  }
  for (mu_it=b_muons.begin(); mu_it!=b_muons.end(); ++mu_it) {
      ex -= (*mu_it).Et() * TMath::Cos((*mu_it).Phi());
      ey -= (*mu_it).Et() * TMath::Sin((*mu_it).Phi());
  }
  for (jet_it=b_jets.begin(); jet_it!=b_jets.end(); ++jet_it) {
    ex -= (*jet_it).Et() * TMath::Cos((*jet_it).Phi());
    ey -= (*jet_it).Et() * TMath::Sin((*jet_it).Phi());
  }
  
  Double_t met_et  = TMath::Hypot(ex, ey);
  Double_t met_phi = TMath::ATan2(ey, ex);
  
  ntuple->met_truth_et = met_et;
  ntuple->met_truth_phi = met_phi;
  
  // Met truth
  ex = met_truth->mpx() * 0.001;
  ey = met_truth->mpy() * 0.001;
  
  ntuple->met_et = TMath::Hypot(ex, ey);
  ntuple->met_phi = TMath::ATan2(ey, ex);
  
  /// Ht
  Double_t sum_jet_pt = 0.;
  for (jet_it=jets.begin(); jet_it!=jets.end(); ++jet_it)
    sum_jet_pt += (*jet_it).Pt();
  
  Double_t ht = sum_jet_pt;
  if (photons.size() > 0)
    ht += photons[0].Pt();
  
  ntuple->ht = ht;
  ntuple->meff = ht + ntuple->met_truth_et;
  
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
      ntuple->weight_pdf2->push_back( LHAPDF::weightxxQ(pdata.pdgId1,pdata.pdgId2,pdata.x1,pdata.x2,pdata.Q,m_pdfs_2[0],m_pdfs_2[iPdf]) );
    }
    
    for (size_t iPdf=0; iPdf<m_pdfs_3.size(); iPdf++) {
      ntuple->weight_pdf3->push_back( LHAPDF::weightxxQ(pdata.pdgId1,pdata.pdgId2,pdata.x1,pdata.x2,pdata.Q,m_pdfs_3[0],m_pdfs_3[iPdf]) );
    }
  }

  // fill ntuple
  if (photons.size()>0)
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


