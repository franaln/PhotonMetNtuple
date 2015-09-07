
#include <iostream>

// ROOT include(s):
#include <TChain.h>
#include <TError.h>
#include <TH1D.h>
#include <TFile.h>

// xAOD include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODRootAccess/tools/Message.h"

// EDM include(s):
#include "xAODMissingET/MissingETContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"

#include "PhotonMetNtuple/Utils.h"
#include "PhotonMetNtuple/TruthUtils.h"
#include "PhotonMetNtuple/MiniWriter.h"

static const char* APP_NAME = "xAODTruthAnalysis";

void show_usage() 
{
  std::cout << APP_NAME << " -i INPUT_FILE -o OUTPUT_FILE" << std::endl;
}

int main(int argc, char *argv[]) 
{

  if (argc < 2) {
    show_usage();
    return 0;
  }

  TString input = "";
  TString output = "";
  for (int i=1; i<argc; i++) {

    if (strcmp("-h", argv[i]) == 0 || strcmp("--help", argv[i]) == 0) {
      show_usage();
      return 0;
    }
    else if (i == argc-1) {
      show_usage();
      return 1;
    }

    // with argument
    else if (strcmp(argv[i], "-i") == 0) {
      input = argv[++i]; //split_cs_string(argv[++i]);
    }
    else if (strcmp(argv[i], "-o") == 0) {
      output = argv[++i];
    }
     else {
       show_usage();
       return 1;
     }
  }
  
  if (input.IsNull() || output.IsNull()) {
    show_usage();
    return 1;
  }
    
  // Initialise the environment for xAOD reading:
  RETURN_CHECK(APP_NAME, xAOD::Init(APP_NAME));

  // Create a TChain with the input file(s):
  TChain chain("CollectionTree");
  chain.Add(input);
  
  // Create a TEvent object with this TChain as input:
  xAOD::TEvent event(xAOD::TEvent::kBranchAccess);
  RETURN_CHECK(APP_NAME, event.readFrom(&chain));
  
  // create "mini" ntuple
  MiniWriter *ntuple = new MiniWriter();
  ntuple->Create(output, "mini");

  // Loop over the input file(s):
  const Long64_t entries = event.getEntries();
  
  for (Long64_t entry = 0; entry < entries; ++entry ) {
    
    // Load the event:
    if (event.getEntry(entry) < 0) {
      Error(APP_NAME, XAOD_MESSAGE("Failed to load entry %i"), static_cast<int>(entry));
      return 1;
    }
    
    const xAOD::TruthParticleContainer* truth_particles = 0;
    RETURN_CHECK(APP_NAME, event.retrieve(truth_particles, "TruthParticles"));

    const xAOD::JetContainer *truth_jets = 0;
    RETURN_CHECK(APP_NAME, event.retrieve(truth_jets, "AntiKt4TruthJets"));    

    const xAOD::MissingETContainer *truth_met = 0;
    RETURN_CHECK(APP_NAME, event.retrieve(truth_met, "MET_Truth"));    

    if (entry%1000 == 0) 
      Info(APP_NAME, "Event Number: %i", static_cast<int>(entry));
    
    ntuple->Clear();

    vector<TruthParticle> b_photons;
    vector<TruthParticle> b_electrons;
    vector<TruthParticle> b_muons;
    vector<TruthParticle> b_jets;

    vector<TruthParticle> photons;
    vector<TruthParticle> electrons;
    vector<TruthParticle> muons;
    vector<TruthParticle> jets;
    
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
      
      // // Jets
      // else if (abspid <= 6) {

      //   if (pt < 20 || abseta > 2.8)
      //     continue;
        
      //   b_jets.push_back(TruthParticle(*truth));

      //   if (pt < 40)
      //     continue;
        
      //   jets.push_back(TruthParticle(*truth));

      //   ntuple->jet_n++;
      //   ntuple->jet_pt->push_back(pt);
      //   ntuple->jet_eta->push_back(eta);
      //   ntuple->jet_phi->push_back(phi);

      // }

    } // end truth particle container loop 


    // loop over truth jets
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

      // if (pt < 40)
      //   continue;
        
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
    vector<TruthParticle>::iterator ph_it;
    vector<TruthParticle>::iterator el_it;
    vector<TruthParticle>::iterator mu_it;
    vector<TruthParticle>::iterator jet_it;

   
    //-- Overlap Removal
    

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
    
    Double_t met_et = TMath::Hypot(ex, ey);
    Double_t met_phi = TMath::ATan2(ey, ex);

    ntuple->met_et = met_et;
    ntuple->met_phi = met_phi;

    // Met truth
    ex = truth_met->at(0)->mpx() * 0.001;
    ey = truth_met->at(0)->mpy() * 0.001;

    ntuple->met_truth_et = TMath::Hypot(ex, ey);
    ntuple->met_truth_phi = TMath::ATan2(ey, ex);

    /// Ht
    Double_t sum_jet_pt = 0.;
    for (jet_it=jets.begin(); jet_it!=jets.end(); ++jet_it)
      sum_jet_pt += (*jet_it).Pt();

    Double_t ht = sum_jet_pt;
    if (photons.size() > 0)
      ht += photons[0].Pt();

    ntuple->ht = ht;

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
    

    // fill ntuple
    ntuple->Fill();

  } // end events loop


  ntuple->Save();
  ntuple->Close();
  
  return 0;
}

