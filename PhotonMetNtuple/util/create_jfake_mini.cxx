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

#include <PhotonMetNtuple/Utils.h>
#include <PhotonMetNtuple/MiniClone.h>

float fjg_factor[2*3] = {
  0.147, 0.133, 0.124, // eta < 1.37:  145<pt<175, 175<pt<225, pt>225
  0.100, 0.096, 0.104  // eta > 1.52:  145<pt<175, 175<pt<225, pt>225
};

float fjg_syst_dn[2*3] = {
  0.035, 0.033, 0.032, // eta < 1.37:  145<pt<175, 175<pt<225, pt>225
  0.033, 0.032, 0.032  // eta > 1.52:  145<pt<175, 175<pt<225, pt>225
};

float fjg_syst_up[2*3] = {
  0.037, 0.036, 0.036, // eta < 1.37:  145<pt<175, 175<pt<225, pt>225
  0.032, 0.031, 0.031  // eta > 1.52:  145<pt<175, 175<pt<225, pt>225
};


unsigned int get_eta_bin(float eta)
{
  float abseta = fabs(eta);

  if (abseta < 1.37)
    return 0;
  else if (abseta > 1.52 && abseta <= 2.37)
    return 1;

  return 99;
}

unsigned int get_pt_bin(float pt)
{
  if (pt > 145. && pt <= 175.)
    return 0;
  else if (pt > 175. && pt <= 225.)
    return 1;
  else if (pt > 225.)
    return 2;

  return 99;
}

  
void loop(TString input_path, TString output_path)
{

  std::cout << "create_jfake_mini: " << input_path << " -> " << output_path << std::endl;
  
  MiniClone *mini = new MiniClone(input_path, output_path);

  mini->Init();

  Int_t total_events = mini->GetEntries();
  if (total_events == 0) {
    std::cout << "no events" << std::endl;
    return;
  }

  float weight_fjg;
  float weight_fjg_dn;
  float weight_fjg_up;

  mini->AddNewBranch("weight_fjg", &weight_fjg);
  mini->AddNewBranch("weight_fjg_dn", &weight_fjg_dn);
  mini->AddNewBranch("weight_fjg_up", &weight_fjg_up);

  // Loop over all events
  int msg_interval = int(total_events/10);
  for (Long64_t jentry=0; jentry<total_events; jentry++) {

    mini->GetEntry(jentry);

    if (total_events>10 && jentry%msg_interval == 0) 
      std::cout << "Processing event " << jentry << " of " << total_events << std::endl;
    
    // skip event with one signal iso photon 
    if (mini->ph_n > 0 && (*mini->ph_pt)[0] > 145.)
      continue;

    // skip event without noniso photons
    if (mini->ph_noniso_n == 0) 
      continue;


    // clear
    mini->Clear();

    // Photon/Electron blocks
    // Interchange ph_noniso <-> photon

    // skip event if non iso photon not in acceptance region
    float phpt  = (*mini->ph_noniso_pt)[0];
    float phphi = (*mini->ph_noniso_phi)[0];
    float pheta = fabs((*mini->ph_noniso_etas2)[0]); // etas2 only for mini >= v41
    
    if (phpt < 145. || pheta > 2.37)
      continue;
    
    if (pheta > 1.37 && pheta < 1.52) // crack region
      continue;

    // only keep events with noniso photons in B region: 5.45 < iso < 29.45
    if ((*mini->ph_noniso_iso)[0] < 5.45 || (*mini->ph_noniso_iso)[0] > 29.45)
      continue;

    mini->new_ph_n = 1;
      
    mini->new_ph_pt->push_back(phpt);
    mini->new_ph_eta->push_back((*mini->ph_noniso_eta)[0]);
    mini->new_ph_etas2->push_back((*mini->ph_noniso_etas2)[0]);
    mini->new_ph_phi->push_back(phphi);
    mini->new_ph_iso->push_back(0.);
    mini->new_ph_trackiso->push_back(0.);
    mini->new_ph_etcone40->push_back(0.);
    mini->new_ph_ptcone20->push_back(0.);
    mini->new_ph_w->push_back(1.);
    
    mini->CopyElectronsBlock();
    mini->CopyMuonsBlock();
    mini->CopyJetsBlock();
    
    mini->CopyEventBlock();
    mini->CopyWeightBlock();
    mini->CopyMetBlock();
    mini->CopyOthersBlock();

    // std::cout << "---" << std::endl;
    // std::cout << "dphi_jetmet mini: " << mini->dphi_jetmet << " copy: " << mini->new_dphi_jetmet << std::endl;
    // std::cout << "dphi_gammet mini: " << mini->dphi_gammet << " copy: " << mini->new_dphi_gammet << std::endl;
    // std::cout << "dphi_gamjet mini: " << mini->dphi_gamjet << " copy: " << mini->new_dphi_gamjet << std::endl;
    // std::cout << "---" << std::endl;
   
    // Replace some variables with noniso photon instead of photon
    mini->new_dphi_gammet = get_dphi(phphi, mini->met_phi);
    
    if (mini->jet_n > 0) 
      mini->new_dphi_gamjet = get_dphi(phphi, (*mini->jet_phi)[0]);

    // std::cout << "dphi_jetmet mini: " << mini->dphi_jetmet << " copy: " << mini->new_dphi_jetmet << std::endl;
    // std::cout << "dphi_gammet mini: " << mini->dphi_gammet << " copy: " << mini->new_dphi_gammet << std::endl;
    // std::cout << "dphi_gamjet mini: " << mini->dphi_gamjet << " copy: " << mini->new_dphi_gamjet << std::endl;
    // std::cout << "---" << std::endl;
   
    mini->new_ht = mini->ht + phpt;
    mini->new_meff = mini->new_ht + mini->met_et;
    
    // f(jet->gam) weight
    unsigned int pt_bin  = get_pt_bin(phpt);
    unsigned int eta_bin = get_eta_bin(pheta);

    if (pt_bin < 3 && eta_bin < 2) {

      unsigned int bin = eta_bin*3+pt_bin;
      
      weight_fjg    = fjg_factor[bin];
      weight_fjg_dn = fjg_factor[bin] - fjg_syst_dn[bin];
      weight_fjg_up = fjg_factor[bin] + fjg_syst_up[bin];
    }
    else
      std::cout << "pt = " << (*mini->ph_noniso_pt)[0] << "(" << pt_bin << "), eta = " << (*mini->ph_noniso_eta)[0] <<"(" << eta_bin << ")" << std::endl;

    mini->Fill();
  }

  mini->Save();
}

int main(int argc, char *argv[]) 
{
  if (argc < 3) {
    std::cout << "usage: create_jfake_mini <input_file> <output_file>" << std::endl;
    return 1;
  }

  TString input_file = argv[1];
  TString output_file = argv[2];

  loop(input_file, output_file);

  return 0;
}
