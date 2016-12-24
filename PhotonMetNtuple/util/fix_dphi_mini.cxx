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


void loop(TString input_path, TString output_path)
{

  std::cout << "fix_mini: " << input_path << " -> " << output_path << std::endl;
  
  MiniClone *mini = new MiniClone(input_path, output_path);

  if (input_path.Contains("mc15_13TeV"))
    mini->SetMC();

  mini->Init();

  Int_t total_events = mini->GetEntries();
  if (total_events == 0) {
    std::cout << "no events" << std::endl;
    return;
  }

  // Loop over all events
  int msg_interval = int(total_events/10);
  for (Long64_t jentry=0; jentry<total_events; jentry++) {

    mini->GetEntry(jentry);

    if (total_events>10 && jentry%msg_interval == 0) 
      std::cout << "Processing event " << jentry << " of " << total_events << std::endl;
    
    mini->Clear();


    // copy blocks
    mini->CopyPhotonsBlock();
    mini->CopyPhotonsNonIsoBlock();

    mini->CopyElectronsBlock();
    mini->CopyElectronsMediumBlock();

    mini->CopyMuonsBlock();
    mini->CopyJetsBlock();
    
    mini->CopyEventBlock();
    mini->CopyWeightBlock();
    mini->CopyMetBlock();
    mini->CopyOthersBlock();
   
    // Fix dphis
    mini->new_dphi_gammet = fabs(mini->dphi_gammet);
    mini->new_dphi_gamjet = fabs(mini->dphi_gamjet);

    mini->new_dphi_met_trackmet = fabs(mini->dphi_met_trackmet);

    // dphi (jet, met)
    Float_t dphi1 = 4.;
    Float_t dphi2 = 4.;

    if (mini->jet_n > 0) dphi1 = get_dphi((*mini->jet_phi)[0], mini->met_phi);
    if (mini->jet_n > 1) dphi2 = get_dphi((*mini->jet_phi)[1], mini->met_phi);

    if (mini->jet_n > 0) {
      mini->new_dphi_jet1met = dphi1;
      mini->new_dphi_jetmet  = TMath::Min(dphi1, dphi2);
    }

    mini->Fill();
  }

  mini->Save();
}

int main(int argc, char *argv[]) 
{
  if (argc < 3) {
    std::cout << "usage: fix_mini <input_file> <output_file>" << std::endl;
    return 1;
  }

  loop(argv[1], argv[2]);

  return 0;
}
