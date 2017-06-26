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
  TH2F  *jfake_ff_dn = (TH2F*)file_ff->Get("jfake_ff_dn");
  TH2F  *jfake_ff_up = (TH2F*)file_ff->Get("jfake_ff_up");

  // Loop over all events
  int msg_interval = int(total_events/10);
  for (Long64_t jentry=0; jentry<total_events; jentry++) {

    mini->GetEntry(jentry);

    if (total_events>10 && jentry%msg_interval == 0) 
      std::cout << "Processing event " << jentry << " of " << total_events << std::endl;
    
    // skip event without non-tight/iso photons
    if (mini->ph_n == 0 || (*mini->ph_pt)[0] < 75.)
      continue;

    // clear
    mini->Clear();

    float phpt  = (*mini->ph_pt)[0];
    float pheta = fabs((*mini->ph_etas2)[0]); // etas2 only for mini >= v41

    mini->CopyAllBlocks();
    
    // f(jet->gam) weight
    unsigned int pt_bin   = jfake_ff->GetXaxis()->FindBin(phpt);  
    unsigned int eta_bin  = jfake_ff->GetYaxis()->FindBin(fabs(pheta)); 

    if (pt_bin > 6)
      pt_bin = 6;

    weight_fjg    = jfake_ff->GetBinContent(pt_bin, eta_bin);
    weight_fjg_dn = jfake_ff_dn->GetBinContent(pt_bin, eta_bin);
    weight_fjg_up = jfake_ff_up->GetBinContent(pt_bin, eta_bin);

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
