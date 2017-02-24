#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
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

  // Load ff histogram (x: pt, y: eta)
  TFile *file_ff = new TFile(fake_rate_path);
  TH2F  *jfake_ff = (TH2F*)file_ff->Get("jfake_ff");

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

    // only keep events with noniso photons in B region: 10.45 < iso < 29.45
    if ((*mini->ph_noniso_iso)[0] < 10.45 || (*mini->ph_noniso_iso)[0] > 29.45)
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

    // Replace some variables with noniso photon instead of photon
    mini->new_dphi_gammet = get_dphi(phphi, mini->met_phi);
    
    if (mini->jet_n > 0) 
      mini->new_dphi_gamjet = get_dphi(phphi, (*mini->jet_phi)[0]);

    mini->new_ht = mini->ht + phpt;
    mini->new_meff = mini->new_ht + mini->met_et;
    
    // f(jet->gam) weight
    unsigned int pt_bin   = jfake_ff->GetXaxis()->FindBin(phpt);  
    unsigned int eta_bin  = jfake_ff->GetYaxis()->FindBin(fabs(pheta)); 

    if (pt_bin > 3)
      pt_bin = 3;

    weight_fjg = jfake_ff->GetBinContent(pt_bin, eta_bin);
    //weight_fjg_dn = fjg_factor[bin] - fjg_syst_dn[bin];
    //weight_fjg_up = fjg_factor[bin] + fjg_syst_up[bin];

    mini->Fill();
  }

  mini->Save();
}

int main(int argc, char *argv[]) 
{
  if (argc < 4) {
    std::cout << "usage: create_jfake_mini <fake_rate_file> <input_file> <output_file>" << std::endl;
    return 1;
  }

  TString fake_rate_file = argv[1];
  TString input_file = argv[2];
  TString output_file = argv[3];

  loop(fake_rate_file, input_file, output_file);

  return 0;
}
