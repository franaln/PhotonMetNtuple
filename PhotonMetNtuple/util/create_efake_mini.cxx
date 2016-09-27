#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>

#include <vector>
#include <iostream>

#include <PhotonMetNtuple/MiniClone.h>

float get_dphi(float phi1, float phi2)
{
  float  phi = fabs(phi1 - phi2);
  if(phi <= TMath::Pi())  return phi;
  else                    return (2 * TMath::Pi() - phi);
}

float feg_factor[]    = {0.012, 0.016, 0.028, 0.045};
float feg_factor_up[] = {0.016, 0.023, 0.037, 0.055};
float feg_factor_dn[] = {0.008, 0.009, 0.019, 0.035};

unsigned int get_eta_bin(double eta)
{
  float abseta = fabs(eta);

  if (abseta <= 0.6)
    return 0;
  else if (abseta > 0.6 && abseta <= 1.37)
    return 1;
  else if (abseta > 1.52 && abseta <= 1.82)
    return 2;
  else if (abseta > 1.82 && abseta <= 2.37)
    return 3;

  return 99;
}

void loop(TString input_path, TString output_path)
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
  

  // Loop over all events
  int msg_interval = int(total_events/10);
  for (Long64_t jentry=0; jentry<total_events; jentry++) {

    mini->GetEntry(jentry);

    if (total_events>10 && jentry%msg_interval == 0) 
      std::cout << "Processing event " << jentry << " of " << total_events << std::endl;
    
    // clear
    mini->Clear();


    // Photon/Electron blocks
    // Interchange el_medium <-> photon

    // skip event with one signal photon
    if (mini->ph_n > 0 && (*mini->ph_pt)[0] > 145.)
      continue;

    // skip event without medium electrons (?)
    if (mini->el_medium_n == 0) 
      continue;

    // skip event if medium electron not in acceptance region
    float eleta = fabs((*mini->el_medium_etas2)[0]);
    float elpt = (*mini->el_medium_pt)[0];
    float elphi = (*mini->el_medium_phi)[0];
    
    if (elpt < 145. || eleta > 2.37)
      continue;
    
    if (eleta > 1.37 && eleta < 1.52)
      continue;
    
    mini->new_ph_n = 1;
    
    mini->new_ph_pt->push_back(elpt);
    mini->new_ph_eta->push_back((*mini->el_medium_eta)[0]);
    mini->new_ph_etas2->push_back((*mini->el_medium_etas2)[0]);
    mini->new_ph_phi->push_back(elphi);
    mini->new_ph_iso->push_back(0.);
    mini->new_ph_iso40->push_back(0.);
    mini->new_ph_w->push_back(1.);
    
    // medium electrons/nominal elecrons matching
    int match_idx = -1;
    for (int i=0; i<mini->el_n; i++) {
      if (((*mini->el_pt)[i] - elpt) < 5 &&  ((*mini->el_phi)[i] - elphi) < 0.1) {
        match_idx = i;
        break;
      }
    }
    
    int el_n = 0;
    for (int i=0; i<mini->el_n; i++) {
      if (i == match_idx)
        continue;
      
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
    
    Float_t dphi1 = 99.;
    Float_t dphi2 = 99.;
    if (mini->jet_n > 0) dphi1 = get_dphi(elphi, (*mini->jet_phi)[0]);
    if (mini->jet_n > 1) dphi2 = get_dphi(elphi, (*mini->jet_phi)[1]);

    mini->new_dphi_jetmet = TMath::Min(dphi1, dphi2);
    if (mini->new_dphi_jetmet > 4.):
      mini->new_dphi_jetmet = -99.;

    mini->new_ht = mini->ht + elpt;
    mini->new_meff = mini->new_ht + mini->met_et;

    // f(e->gam) weight
    unsigned int bin = get_eta_bin(eleta);
    if (bin < 4) {
      weight_feg    = feg_factor[bin];
      weight_feg_dn = feg_factor_dn[bin];
      weight_feg_up = feg_factor_up[bin];
    }
  
    mini->Fill();
  }

  mini->Save();
}

int main(int argc, char *argv[]) 
{
  if (argc < 3) {
    std::cout << "usage: create_efake_mini <input_file> <output_file>" << std::endl;
    return 1;
  }

  TString input_file = argv[1];
  TString output_file = argv[2];

  loop(input_file, output_file);

  return 0;
}
