#include <TROOT.h>
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>

#include <vector>
#include <iostream>

#include <PhotonMetNtuple/Utils.h>
#include <PhotonMetNtuple/MiniClone.h>

void loop(TString fake_rate_path, TString input_path, TString output_path)
{
  std::cout << "create_efake_mini: " << input_path << " -> " << output_path << std::endl;
  
  MiniClone *mini = new MiniClone(input_path, output_path);

  mini->Init();

  Int_t total_events = mini->GetEntries();
  if (total_events == 0) {
    std::cout << "no events" << std::endl;
    return;
  }

  float weight_feg;
  float weight_feg_dn;
  float weight_feg_up;

  mini->AddNewBranch("weight_feg", &weight_feg);
  mini->AddNewBranch("weight_feg_up", &weight_feg_up);
  mini->AddNewBranch("weight_feg_dn", &weight_feg_dn);

  // Load ff histograms (4 eta bins)
  TFile *file_ff = new TFile(fake_rate_path);

  TH2F *efake_ff    = (TH2F*)file_ff->Get("efake_ff");
  TH2F *efake_ff_up = (TH2F*)file_ff->Get("efake_ff_up");
  TH2F *efake_ff_dn = (TH2F*)file_ff->Get("efake_ff_dn");
  

  // Loop over all events
  int msg_interval = int(total_events/10);
  for (Long64_t jentry=0; jentry<total_events; jentry++) {

    mini->GetEntry(jentry);

    if (total_events>10 && jentry%msg_interval == 0) 
      std::cout << "Processing event " << jentry << " of " << total_events << std::endl;
    

    // skip event with one signal photon
    if (mini->ph_n > 0 && (*mini->ph_pt)[0] > 145.)
      continue;

    // skip events without electrons (?)
    if (mini->el_n == 0)
      continue;

    // clear
    mini->Clear();

    // Photon/Electron blocks
    // Interchange el_medium <-> photon

    // skip event if medium electron not in acceptance region
    float eleta = fabs((*mini->el_etas2)[0]);
    float elpt  = (*mini->el_pt)[0];
    float elphi = (*mini->el_phi)[0];
    
    if (elpt < 145. || eleta > 2.37)
      continue;
    
    if (eleta > 1.37 && eleta < 1.52)
      continue;
    
    mini->new_ph_n = 1;
    
    mini->new_ph_pt->push_back(elpt);
    mini->new_ph_eta->push_back((*mini->el_eta)[0]);
    mini->new_ph_etas2->push_back((*mini->el_etas2)[0]);
    mini->new_ph_phi->push_back(elphi);
    mini->new_ph_iso->push_back(0.);
    mini->new_ph_trackiso->push_back(0.);
    mini->new_ph_etcone40->push_back(0.);
    mini->new_ph_ptcone20->push_back(0.);
    mini->new_ph_w->push_back(1.);
    
    // // medium electrons/nominal elecrons matching
    // int match_idx = -1;
    // for (int i=0; i<mini->el_n; i++) {
    //   if (((*mini->el_pt)[i] - elpt) < 1 &&  ((*mini->el_phi)[i] - elphi) < 0.01 && (fabs((*mini->el_eta)[i]) - eleta) < 0.01) {
    //     match_idx = i;
    //     break;
    //   }
    // }
    
    int el_n = 0;
    for (int i=1; i<mini->el_n; i++) {
      el_n += 1;
      mini->new_el_pt->push_back((*mini->el_pt)[i]);
      mini->new_el_eta->push_back((*mini->el_eta)[i]);
      mini->new_el_phi->push_back((*mini->el_phi)[i]);
      mini->new_el_ch->push_back((*mini->el_ch)[i]);
      mini->new_el_w->push_back((*mini->el_w)[i]);
    }
    mini->new_el_n = el_n;

    mini->CopyMuonsBlock();
    mini->CopyJetsBlock();
    
    mini->CopyEventBlock();
    mini->CopyWeightBlock();
    mini->CopyMetBlock();
    mini->CopyOthersBlock();
    
    // Replace some variables with electron instead of photon
    mini->new_dphi_gammet = get_dphi(elphi, mini->met_phi);
    
    if (mini->jet_n > 0) mini->new_dphi_gamjet = get_dphi(elphi, (*mini->jet_phi)[0]);
    
    mini->new_ht = mini->ht + elpt;
    mini->new_meff = mini->new_ht + mini->met_et;

    Float_t mt2_gam = 2 * mini->met_et * elpt * (1 - TMath::Cos(mini->dphi_gammet));
    mini->mt_gam = TMath::Sqrt(mt2_gam);

    // f(e->gam) weight
    unsigned int pt_bin = efake_ff->GetYaxis()->FindBin(elpt); 
    unsigned int eta_bin = efake_ff->GetXaxis()->FindBin(eleta); 
    
    weight_feg    = efake_ff   ->GetBinContent(eta_bin, pt_bin); 
    weight_feg_dn = efake_ff_dn->GetBinContent(eta_bin, pt_bin);
    weight_feg_up = efake_ff_up->GetBinContent(eta_bin, pt_bin);
  
    mini->Fill();
  }

  mini->Save();
}

int main(int argc, char *argv[]) 
{
  if (argc < 3) {
    std::cout << "usage: create_efake_mini <fake_rate_file> <input_file> <output_file>" << std::endl;
    return 1;
  }

  TString fake_rate_file = argv[1];
  TString input_file = argv[2];
  TString output_file = argv[3];

  loop(fake_rate_file, input_file, output_file);

  return 0;
}
