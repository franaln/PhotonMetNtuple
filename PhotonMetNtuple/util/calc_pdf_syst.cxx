#include <vector>
#include <iostream>
#include <iomanip>

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TH1F.h>

#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Reweighting.h"


void fill_histograms(std::vector<TString> paths, 
                     std::vector<double> xs,
                     double lumi,
                     std::map<TString, TH1F*> &h_pdf1_sel, 
                     std::map<TString, TH1F*> &h_pdf2_sel, 
                     std::map<TString, TH1F*> &h_pdf3_sel) 
{
  
  for(auto const &e : h_pdf1_sel)
    e.second->Reset();
  
  for(auto const &e : h_pdf2_sel)
    e.second->Reset();
  
  for(auto const &e : h_pdf3_sel)
    e.second->Reset();

  int counter = 0;
  for (auto path : paths) {

    TFile *f = TFile::Open(path);
    TH1F *events =  (TH1F*)f->Get("events");

    Double_t sumw = events->GetBinContent(3);
    Double_t weight = (xs[counter] * lumi) / sumw;

    f->Close();


    TChain *tree = new TChain("mini");
    tree->Add(path);
    
    Int_t total_events = tree->GetEntries();
    
    if (total_events == 0)
      return;
    
    TTreeReader reader(tree); 
    
    TTreeReaderValue<int> ph_n(reader, "ph_n");
    TTreeReaderValue<int> el_n(reader, "el_n");
    TTreeReaderValue<int> mu_n(reader, "mu_n");
    TTreeReaderValue<int> jet_n(reader, "jet_n");
    TTreeReaderValue<int> bjet_n(reader, "bjet_n");

    TTreeReaderValue<std::vector<float> > ph_pt(reader, "ph_pt");
    
    TTreeReaderValue<float> met_et(reader, "met_et");
    TTreeReaderValue<float> meff(reader, "meff");
    TTreeReaderValue<float> rt4(reader, "rt4");
    TTreeReaderValue<float> dphi_jetmet(reader, "dphi_jetmet");
    TTreeReaderValue<float> dphi_gamjet(reader, "dphi_gamjet");
    TTreeReaderValue<float> dphi_gammet(reader, "dphi_gammet");
    
    TTreeReaderValue< std::vector<float> > weight_pdf1(reader, "weight_pdf1");
    TTreeReaderValue< std::vector<float> > weight_pdf2(reader, "weight_pdf2");
    TTreeReaderValue< std::vector<float> > weight_pdf3(reader, "weight_pdf3");
    

    // Loop over all entries of the TTree or TChain.
    while (reader.Next()) {
      
      unsigned int n_pdf1 = (*weight_pdf1).size();
      unsigned int n_pdf2 = (*weight_pdf2).size();
      unsigned int n_pdf3 = (*weight_pdf3).size();

      // presel
      if ((*ph_n)>0 && (*el_n+*mu_n)==0 && (*ph_pt)[0]>145. && (*met_et)>200. && (*jet_n)>1) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["presel"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["presel"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["presel"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }

      // SRL
      if ((*ph_n)>0 && (*el_n+*mu_n)==0 && (*ph_pt)[0]>145. && (*met_et)>200. && (*jet_n)>4 && (*dphi_jetmet)>0.4 && (*dphi_gammet)>0.4 && (*meff)>2000. && (*rt4)<0.9) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["SRL"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["SRL"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["SRL"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }
      
      // SRH
      if ((*ph_n)>0 && (*el_n+*mu_n)==0 && (*ph_pt)[0]>400. && (*met_et)>400. && (*jet_n)>2 && (*dphi_jetmet)>0.4 && (*dphi_gammet)>0.4 && (*meff)>2000.) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["SRH"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["SRH"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["SRH"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }
      
      // CRQ
      if ((*ph_n)>0 && (*el_n+*mu_n)==0 && (*ph_pt)[0]>145. && (*jet_n)>2 && (*met_et)>100. && (*dphi_jetmet)<0.4 && (*dphi_gammet)>0.4 && (*meff)>2000.) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["CRQ"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["CRQ"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["CRQ"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }
      
      // CRW
      if ((*ph_n)>=1 && (*el_n+*mu_n)>=1 && (*ph_pt)[0]>145. && (*jet_n)>=1 && (*met_et)>100. && (*met_et)<200. && (*dphi_jetmet)>0.4 && (*meff)>500. && (*bjet_n)==0) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["CRW"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["CRW"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["CRW"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }

      // CRT
      if ((*ph_n)>=1 && (*el_n+*mu_n)>=1 && (*ph_pt)[0]>145. && (*jet_n)>=2 && (*met_et)>50. && (*met_et)<200. && (*dphi_jetmet)>0.4 && (*meff)>500. && (*bjet_n)>=2) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["CRT"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["CRT"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["CRT"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }

      // VRL1
      if ((*ph_n)>=1 && (*el_n+*mu_n)>=1 && (*ph_pt)[0]>145. && (*jet_n)>=2 && (*met_et)<200. && (*dphi_jetmet)>0.4 && (*meff)>1000.) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["VRL1"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["VRL1"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["VRL1"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }

      // VRL2
      if ((*ph_n)>=1 && (*el_n+*mu_n)>=1 && (*ph_pt)[0]>145. && (*jet_n)>=2 && (*met_et)<200. && (*dphi_jetmet)>0.4 && (*meff)>1500.) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["VRL2"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["VRL2"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["VRL2"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }

      // VRL3
      if ((*ph_n)>=1 && (*el_n+*mu_n)>=1 && (*ph_pt)[0]>145. && (*jet_n)>=2 && (*met_et)>200. && (*dphi_jetmet)>0.4 && (*meff)>1000. && (*meff)<2000.) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["VRL3"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["VRL3"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["VRL3"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }

      // VRL4
      if ((*ph_n)>=1 && (*el_n+*mu_n)>=1 && (*ph_pt)[0]>145. && (*jet_n)>=2 && (*met_et)>200. && (*dphi_jetmet)<0.4 && (*meff)>1500.) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["VRL4"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["VRL4"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["VRL4"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }

      // VRZ
      if ((*ph_n)>=1 && ((*el_n)==2 || (*mu_n==2)) && (*ph_pt)[0]>145. && (*jet_n)>=1 && (*met_et)<200. && (*meff)>1000. && (*bjet_n)==0) {
        for (size_t i=0; i< n_pdf1; i++)
          h_pdf1_sel["VRZ"]->Fill(i+0.5, (*weight_pdf1)[i] * weight);
        
        for (size_t i=0; i< n_pdf2; i++)
          h_pdf2_sel["VRZ"]->Fill(i+0.5, (*weight_pdf2)[i] * weight);
        
        for (size_t i=0; i< n_pdf3; i++)
          h_pdf3_sel["VRZ"]->Fill(i+0.5, (*weight_pdf3)[i] * weight);
      }

    }

    counter++;
  }

  return;
}


void calc_syst(TString sample, 
               std::vector<TString> regions, 
               std::map<TString, TH1F*> &h_pdf1_sel, 
               std::map<TString, TH1F*> &h_pdf2_sel, 
               std::map<TString, TH1F*> &h_pdf3_sel
               )
{
  
  std::string pdf1 = "CT10";
  std::string pdf2 = "NNPDF30_lo_as_0130";
  std::string pdf3 = "MMHT2014lo68cl";
  
  const LHAPDF::PDFSet pdf_set_1(pdf1);
  const LHAPDF::PDFSet pdf_set_2(pdf2);
  const LHAPDF::PDFSet pdf_set_3(pdf3);


  if (sample == "ttgamma")
    std::cout << "Nominal: PDF2\n" << std::endl;
  else if (sample == "wgamma" || sample == "zgamma" || sample == "znunugamma")
    std::cout << "Nominal: PDF1\n" << std::endl;
  
  for (auto reg : regions) {
    
    std::vector<double> pdf_yields_1;
    std::vector<double> pdf_yields_2;
    std::vector<double> pdf_yields_3;
    
    for (int i=1; i<=h_pdf1_sel[reg]->GetNbinsX(); i++)
      pdf_yields_1.push_back(h_pdf1_sel[reg]->GetBinContent(i));
    for (int i=1; i<=h_pdf2_sel[reg]->GetNbinsX(); i++)
      pdf_yields_2.push_back(h_pdf2_sel[reg]->GetBinContent(i));
    for (int i=1; i<=h_pdf3_sel[reg]->GetNbinsX(); i++)
      pdf_yields_3.push_back(h_pdf3_sel[reg]->GetBinContent(i));
    
    const LHAPDF::PDFUncertainty yield_uncertainty_1 = pdf_set_1.uncertainty(pdf_yields_1, 68); // all at 68% CL
    const LHAPDF::PDFUncertainty yield_uncertainty_2 = pdf_set_2.uncertainty(pdf_yields_2, 68);
    const LHAPDF::PDFUncertainty yield_uncertainty_3 = pdf_set_3.uncertainty(pdf_yields_3, 68);

    // CT10
    double nom_1  = yield_uncertainty_1.central;
    double up_1   = yield_uncertainty_1.errplus;
    double dn_1   = yield_uncertainty_1.errminus;
    
    // NNPDF
    double nom_2 = yield_uncertainty_2.central;
    double up_2  = yield_uncertainty_2.errplus;
    double dn_2  = yield_uncertainty_2.errminus;

    // MMHT2014
    double nom_3 = yield_uncertainty_3.central;
    double up_3  = yield_uncertainty_3.errplus;
    double dn_3  = yield_uncertainty_3.errminus;
    
    double error_combined = 0.5 * (std::max(up_1, std::max(up_2, up_3)) + std::max(dn_1, std::max(dn_2, dn_3)));

    double nominal = 0.;
    if (sample == "ttgamma")
      nominal = nom_2; // NNPDF
    else if (sample == "wgamma" || sample == "zgamma" || sample == "znunugamma")
      nominal = nom_1; // CT10


    std::cout << "Region: " << reg << std::endl;
    std::cout << "Stat error sqrt(nom): " << sqrt(nominal) << std::endl;
    std::cout << "PDF1     : " << nom_1 << " -" << dn_1 << " (" << 100*dn_1/nominal << "%) " << " +" << up_1  << " (" << 100*up_1/nominal << "%)" << std::endl;
    std::cout << "PDF2     : " << nom_2 << " -" << dn_2 << " (" << 100*dn_2/nominal << "%) " << " +" << up_2  << " (" << 100*up_2/nominal << "%)" << std::endl;
    std::cout << "PDF3     : " << nom_3 << " -" << dn_3 << " (" << 100*dn_3/nominal << "%) " << " +" << up_3  << " (" << 100*up_3/nominal << "%)" << std::endl;
    std::cout << "Combined : " << error_combined << " (" << 100*error_combined/nominal << "%)" <<  std::endl;
  }
}


int main()
{
  std::vector<TString> regions = {"presel", "CRQ", "CRW", "CRT", "VRL1", "VRL2", "VRL3", "VRL4", "VRZ", "SRL", "SRH"};
  
  double lumi = 36500.;
  
  unsigned int n_pdf1 = 53;
  unsigned int n_pdf2 = 101;
  unsigned int n_pdf3 = 51;

  std::map<TString, TH1F*> h_pdf1_sel;
  std::map<TString, TH1F*> h_pdf2_sel;
  std::map<TString, TH1F*> h_pdf3_sel;
  
  for (auto reg : regions) {
    h_pdf1_sel[reg] = new TH1F("h_pdf1_"+reg, "h_pdf1_"+reg, n_pdf1, 0, n_pdf1);
    h_pdf2_sel[reg] = new TH1F("h_pdf2_"+reg, "h_pdf2_"+reg, n_pdf2, 0, n_pdf2);
    h_pdf3_sel[reg] = new TH1F("h_pdf3_"+reg, "h_pdf3_"+reg, n_pdf3, 0, n_pdf3);
  }
  
  // ttgam
  std::cout << "\n PDF systematics for tt + gamma\n" << std::endl;
  
  std::vector<TString> paths = {
    "/raid/falonso/mini2/v51/mc15_13TeV.410084.MadGraphPythia8EvtGen_A14NNPDF23LO_ttgamma80_noallhad.truth.v51_output.root"
  };
  
  std::vector<Double_t> xs = { 0.493805781 };

  fill_histograms(paths, xs, lumi, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);
  calc_syst("ttgamma", regions, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);


  // ttgam NLO
  std::cout << "\n PDF systematics for tt + gamma NLO\n" << std::endl;
  
  paths = {
    "/raid/falonso/mini2/v51/mc15_13TeV.407320.aMcAtNloPythia8EvtGen_MEN30NLO_A14N23LO_tta140.truth.v51_output.root"
  };
  
  xs = { 0.21505 };

  fill_histograms(paths, xs, lumi, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);
  calc_syst("ttgamma", regions, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);


  //  wgamma
  std::cout << "\n PDF systematics for W + gamma\n" << std::endl;
  
  paths = {
    "/raid/falonso/mini2/v51/mc15_13TeV.301890.Sherpa_CT10_enugammaPt35_70.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301891.Sherpa_CT10_enugammaPt70_140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301892.Sherpa_CT10_enugammaPt140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301893.Sherpa_CT10_munugammaPt35_70.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301894.Sherpa_CT10_munugammaPt70_140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301895.Sherpa_CT10_munugammaPt140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301896.Sherpa_CT10_taunugammaPt35_70.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301897.Sherpa_CT10_taunugammaPt70_140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301898.Sherpa_CT10_taunugammaPt140.truth.v51_output.root"
  };

  xs = {
    15.348,
    1.5282,
    0.24155,
    15.272,
    1.5235,
    0.24183,
    15.297,
    1.529,
    0.2426
  };

  fill_histograms(paths, xs, lumi, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);
  calc_syst("wgamma", regions, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);


  // Z+gamma
  std::cout << "\n PDF systematics for Z+gamma\n" << std::endl;
  
  paths = {
    "/raid/falonso/mini2/v51/mc15_13TeV.301899.Sherpa_CT10_eegammaPt35_70.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301900.Sherpa_CT10_eegammaPt70_140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301901.Sherpa_CT10_eegammaPt140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301902.Sherpa_CT10_mumugammaPt35_70.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301903.Sherpa_CT10_mumugammaPt70_140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301904.Sherpa_CT10_mumugammaPt140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301905.Sherpa_CT10_tautaugammaPt35_70.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301906.Sherpa_CT10_tautaugammaPt70_140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301907.Sherpa_CT10_tautaugammaPt140.truth.v51_output.root"
  };

  xs = {
    5.242,
    0.38455,
    0.047209,
    5.2455,
    0.38548,
    0.04724,
    5.249,
    0.38482,
    0.047025
  };

  fill_histograms(paths, xs, lumi, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);
  calc_syst("zgamma", regions, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);

  //  Z(nunu)+gamma
  std::cout << "\n PDF systematics for Z(nunu) + gamma\n" << std::endl;
  
  paths = {
    "/raid/falonso/mini2/v51/mc15_13TeV.301908.Sherpa_CT10_nunugammaPt35_70.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301909.Sherpa_CT10_nunugammaPt70_140.truth.v51_output.root",
    "/raid/falonso/mini2/v51/mc15_13TeV.301910.Sherpa_CT10_nunugammaPt140.truth.v51_output.root"
  };

  xs = {
    4.0365,
    0.97151,
    0.17115
  };

  fill_histograms(paths, xs, lumi, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);
  calc_syst("znunugamma", regions, h_pdf1_sel, h_pdf2_sel, h_pdf3_sel);
}
