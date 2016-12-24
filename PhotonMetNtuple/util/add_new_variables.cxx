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
  std::cout << "add_new_variables: " << input_path << " -> " << output_path << std::endl;

  MiniClone *mini = new MiniClone(input_path, output_path);

  if (input_path.Contains("mc15"))
    mini->SetMC();
  else if (input_path.Contains("efake"))
    mini->SetEfakeSample();
  else if (input_path.Contains("jfake"))
    mini->SetJfakeSample();
  
  mini->Init();

  Int_t total_events = mini->GetEntries();
  if (total_events == 0) {
    std::cout << "no events" << std::endl;
    return;
  }

  // New variables
  float ht0;
  float ht2;
  float meff0;
  float meff2;
  float stgam;

  float dphi_gammet2;
  float dphi_gamjet1;
  float dphi_gamjet2;
  float dphi_gamjet3;
  float dphi_gamjet4;
  float dphi_jetmet1;
  float dphi_jetmet2;
  float dphi_jetmet3;
  float dphi_jetmet4;

  float mt_gam;
  float mt_lep;

  float deta_gamjet;
  float deta_gamjet2;
  float dr_gamjet;
  float dr_gamjet2;
  float eta2_gamjet;
  float dpt_gamjet;

  float met_sqrtht;

  // New branches
  mini->AddNewBranch("ht0", &ht0);				
  mini->AddNewBranch("ht2", &ht2);				
  mini->AddNewBranch("meff0", &meff0);				
  mini->AddNewBranch("meff2", &meff2);				
  mini->AddNewBranch("stgam", &stgam);				

  mini->AddNewBranch("dphi_gammet2", &dphi_gammet2);
  mini->AddNewBranch("dphi_jetmet1", &dphi_jetmet1);
  mini->AddNewBranch("dphi_jetmet2", &dphi_jetmet2);
  mini->AddNewBranch("dphi_jetmet3", &dphi_jetmet3);  
  mini->AddNewBranch("dphi_jetmet4", &dphi_jetmet4);
  mini->AddNewBranch("dphi_gamjet1", &dphi_gamjet1);
  mini->AddNewBranch("dphi_gamjet2", &dphi_gamjet2);
  mini->AddNewBranch("dphi_gamjet3", &dphi_gamjet3);  
  mini->AddNewBranch("dphi_gamjet4", &dphi_gamjet4);

  mini->AddNewBranch("mt_gam", &mt_gam);
  mini->AddNewBranch("mt_lep", &mt_lep);

  mini->AddNewBranch("met_sqrtht", &met_sqrtht);


  mini->AddNewBranch("deta_gamjet", &deta_gamjet);
  mini->AddNewBranch("deta_gamjet2", &deta_gamjet2);
  mini->AddNewBranch("dr_gamjet", &dr_gamjet);
  mini->AddNewBranch("dr_gamjet2", &dr_gamjet2);
  mini->AddNewBranch("eta2_gamjet", &eta2_gamjet);
  mini->AddNewBranch("dpt_gamjet", &dpt_gamjet);
  
  // Loop over all events
  int msg_interval = int(total_events/10);
  for (Long64_t jentry=0; jentry<total_events; jentry++) {

    mini->GetEntry(jentry);

    if (total_events>10 && jentry%msg_interval == 0) 
      std::cout << "Processing event " << jentry << " of " << total_events << std::endl;


    // skip events without photons
    if (mini->ph_n < 1)
      continue;

    // clear
    mini->Clear();

    mini->CopyPhotonsBlock();
    mini->CopyElectronsBlock();
    mini->CopyMuonsBlock();
    mini->CopyJetsBlock();

    // mini->CopyPhotonsNonIsoBlock();
    // mini->CopyElectronsMediumBlock();
    
    mini->CopyEventBlock();
    mini->CopyWeightBlock();
    mini->CopyMetBlock();
    mini->CopyOthersBlock();

    // New variables
    float sum_jet_pt = 0.;
    for (auto jetpt : (*mini->jet_pt))
      sum_jet_pt += jetpt;

    //// Ht definition
    ht0 = sum_jet_pt;
    ht2 = sum_jet_pt;
    if (mini->ph_n > 0) ht2 += mini->ph_pt->at(0);
    if (mini->ph_n > 1) ht2 += mini->ph_pt->at(0);

    meff0 = ht0 + mini->met_et;
    meff2 = ht2 + mini->met_et;

    //// Stgam
    stgam = mini->met_et;
    if (mini->ph_n > 0) stgam += mini->ph_pt->at(0);

    // MET / sqrt(HT)
    met_sqrtht = mini->met_et / TMath::Sqrt(mini->ht);

    //// Dphi definition

    //// fix common definitions
    if (mini->ph_n == 0 || mini->jet_n == 0) mini->new_dphi_gamjet = -99.;
    if (mini->jet_n == 0) mini->new_dphi_jetmet = -99.;
    if (mini->ph_n == 0) mini->new_dphi_gammet = -99.;

    dphi_gamjet1 = 4.;
    dphi_gamjet2 = 4.;
    dphi_gamjet3 = 4.;
    dphi_gamjet4 = 4.;
    
    dphi_jetmet1 = 4.;
    dphi_jetmet2 = 4.;
    dphi_jetmet3 = 4.;
    dphi_jetmet4 = 4.;

    //// photon/met    
    dphi_gammet2 = mini->dphi_gammet;
    if (mini->ph_n > 1)
      dphi_gammet2 = TMath::Min(mini->dphi_gammet, get_dphi((*mini->ph_phi)[1], mini->met_phi));

    //// photon/jet
    Float_t dphi1 = 4.;
    Float_t dphi2 = 4.;
    Float_t dphi3 = 4.;
    Float_t dphi4 = 4.;
    if (mini->jet_n > 0) dphi1 = get_dphi((*mini->ph_phi)[0], (*mini->jet_phi)[0]);
    if (mini->jet_n > 1) dphi2 = get_dphi((*mini->ph_phi)[0], (*mini->jet_phi)[1]);
    if (mini->jet_n > 2) dphi3 = get_dphi((*mini->ph_phi)[0], (*mini->jet_phi)[2]);
    if (mini->jet_n > 3) dphi4 = get_dphi((*mini->ph_phi)[0], (*mini->jet_phi)[3]);
    
    dphi_gamjet1 = dphi1;
    dphi_gamjet2 = TMath::Min(dphi_gamjet1, dphi2);
    dphi_gamjet3 = TMath::Min(dphi_gamjet2, dphi3);
    dphi_gamjet4 = TMath::Min(dphi_gamjet3, dphi4);

    //// jet/met
    dphi1 = 4.;
    dphi2 = 4.;
    dphi3 = 4.;
    dphi4 = 4.;
    if (mini->jet_n > 0) dphi1 = get_dphi(mini->met_phi, (*mini->jet_phi)[0]);
    if (mini->jet_n > 1) dphi2 = get_dphi(mini->met_phi, (*mini->jet_phi)[1]);
    if (mini->jet_n > 2) dphi3 = get_dphi(mini->met_phi, (*mini->jet_phi)[2]);
    if (mini->jet_n > 3) dphi4 = get_dphi(mini->met_phi, (*mini->jet_phi)[3]);
    
    dphi_jetmet1 = dphi1;
    dphi_jetmet2 = TMath::Min(dphi_jetmet1, dphi2);
    dphi_jetmet3 = TMath::Min(dphi_jetmet2, dphi3);
    dphi_jetmet4 = TMath::Min(dphi_jetmet3, dphi4);

    //// Mt 
    Float_t mt2_gam = 2 * mini->met_et * mini->ph_pt->at(0) * (1 - TMath::Cos(mini->dphi_gammet));
    mt_gam = TMath::Sqrt(mt2_gam);

    Float_t mt2_lep = -99;
    if (mini->el_n > 0 && mini->mu_n > 0) {
      if ((*mini->el_pt)[0] > (*mini->mu_pt)[0])
        mt2_lep = 2 * mini->met_et * mini->el_pt->at(0) * (1 - TMath::Cos(get_dphi(mini->el_phi->at(0), mini->met_et)));
      else
        mt2_lep = 2 * mini->met_et * mini->mu_pt->at(0) * (1 - TMath::Cos(get_dphi(mini->mu_phi->at(0), mini->met_et)));
    }
    else if (mini->el_n > 0)
        mt2_lep = 2 * mini->met_et * mini->el_pt->at(0) * (1 - TMath::Cos(get_dphi(mini->el_phi->at(0), mini->met_et)));
    else if (mini->mu_n > 0)
        mt2_lep = 2 * mini->met_et * mini->mu_pt->at(0) * (1 - TMath::Cos(get_dphi(mini->mu_phi->at(0), mini->met_et)));

    mt_lep = TMath::Sqrt(mt2_lep);

    //// Separation photon-jet
    deta_gamjet = -99.;
    dr_gamjet   = -99.;
    if (mini->jet_n > 0) {
      deta_gamjet = get_deta(mini->ph_etas2->at(0), mini->jet_eta->at(0));
      dr_gamjet   = get_dr(mini->ph_etas2->at(0), mini->ph_phi->at(0), mini->jet_eta->at(0), mini->jet_phi->at(0));

      eta2_gamjet = (*mini->ph_etas2)[0] * (*mini->jet_eta)[0];

      dpt_gamjet = TMath::Abs((*mini->ph_pt)[0] - (*mini->jet_pt)[0]);
    }

    deta_gamjet2 = deta_gamjet;
    dr_gamjet2 = dr_gamjet;
    if (mini->jet_n > 1) {
      deta_gamjet2 = TMath::Min(deta_gamjet, get_deta(mini->ph_eta->at(0), mini->jet_eta->at(0)));
      dr_gamjet2   = TMath::Min(dr_gamjet, get_dr(mini->ph_eta->at(0), mini->ph_phi->at(0), mini->jet_eta->at(0), mini->jet_phi->at(0)));
    }


    mini->Fill();
  }


  mini->Save();
}

int main(int argc, char *argv[]) 
{
  if (argc < 3) {
    std::cout << "usage: add_new_variables <input_file> <output_file>" << std::endl;
    return 1;
  }

  TString input_file = argv[1];
  TString output_file = argv[2];

  loop(input_file, output_file);

  return 0;
}
