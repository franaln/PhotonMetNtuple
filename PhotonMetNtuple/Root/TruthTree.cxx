#define TruthTree_cxx
 
#include "PhotonMetNtuple/TruthTree.h"


TruthTree::TruthTree(const std::string& name): asg::AsgMetadataTool( name )
{
  declareProperty("OutFile", m_outfile); 
  declareProperty("SavePDF", m_pdfrw); 

  tree = 0;
}

TruthTree::~TruthTree()
{
  delete ph_pt;
  delete ph_iso;
  delete ph_eta;
  delete ph_phi;
  delete ph_type;
  delete el_pt;
  delete el_eta;
  delete el_phi;
  delete mu_pt;
  delete mu_eta;
  delete mu_phi;
  delete jet_pt;
  delete jet_eta;
  delete jet_phi;
}

StatusCode TruthTree::initialize()
{
  CreateTree("mini");
  CreateBranches();

  return StatusCode::SUCCESS;
}

void TruthTree::CreateTree(TString tree_name)
{
  tree = new TTree(tree_name, tree_name);
  tree->SetDirectory(m_outfile);
}

void TruthTree::CreateBranches()
{
  tree->Branch("ph_n", &ph_n);
  tree->Branch("el_n", &el_n);
  tree->Branch("mu_n", &mu_n);
  tree->Branch("jet_n", &jet_n);	        
  tree->Branch("bjet_n", &bjet_n);	        

  ph_pt = new std::vector<float>;
  ph_iso = new std::vector<float>;
  ph_eta = new std::vector<float>;
  ph_phi = new std::vector<float>;
  ph_type = new std::vector<int>;
  tree->Branch("ph_pt", &ph_pt);
  tree->Branch("ph_iso", &ph_iso);
  tree->Branch("ph_eta", &ph_eta);
  tree->Branch("ph_phi", &ph_phi);
  tree->Branch("ph_type", &ph_type);

  el_pt  = new std::vector<float>;
  el_eta = new std::vector<float>;
  el_phi = new std::vector<float>;
  tree->Branch("el_pt", &el_pt);
  tree->Branch("el_eta", &el_eta);
  tree->Branch("el_phi", &el_phi);
  
  mu_pt  = new std::vector<float>;
  mu_eta = new std::vector<float>;
  mu_phi = new std::vector<float>;
  tree->Branch("mu_pt", &mu_pt);
  tree->Branch("mu_eta", &mu_eta);
  tree->Branch("mu_phi", &mu_phi);

  jet_pt = new std::vector<float>;
  jet_eta = new std::vector<float>;
  jet_phi = new std::vector<float>;
  tree->Branch("jet_pt", &jet_pt);	      
  tree->Branch("jet_eta", &jet_eta);	      
  tree->Branch("jet_phi", &jet_phi);	      

  tree->Branch("met_et", &met_et);
  tree->Branch("met_phi", &met_phi);
  tree->Branch("met_sumet", &met_sumet);

  tree->Branch("met_truth_et", &met_truth_et);
  tree->Branch("met_truth_phi", &met_truth_phi);

  tree->Branch("meff", &meff);				
  tree->Branch("ht", &ht);				
  tree->Branch("rt2", &rt2);				
  tree->Branch("rt4", &rt4);				

  tree->Branch("dphi_gamjet", &dphi_gamjet);
  tree->Branch("dphi_jetmet", &dphi_jetmet);
  tree->Branch("dphi_gammet", &dphi_gammet);

  // tree->Branch("event", &event);
  // tree->Branch("smeared", &smeared);
  // tree->Branch("weight", &weight);

  // btag_weight = new vector<float>;
  // tree->Branch("btag_weight", &btag_weight);

  // tree->Branch("weight_prwdn", &weight_prwdn);
  // tree->Branch("weight_prwup", &weight_prwup);

  // tree->Branch("weight_feg", &weight_feg);

  if (m_pdfrw) {
    weight_pdf1 = new std::vector<float>;
    weight_pdf2 = new std::vector<float>;
    weight_pdf3 = new std::vector<float>;
    tree->Branch("weight_pdf1", &weight_pdf1);
    tree->Branch("weight_pdf2", &weight_pdf2);
    tree->Branch("weight_pdf3", &weight_pdf3);
  }

  Clear();
}

void TruthTree::Clear()
{
  ph_n = 0;
  el_n = 0;
  mu_n = 0;
  jet_n = 0;
  bjet_n = 0;

  ph_pt->clear();
  ph_iso->clear();
  ph_eta->clear();
  ph_phi->clear();
  ph_type->clear();

  el_pt->clear();
  el_eta->clear();
  el_phi->clear();

  mu_pt->clear();
  mu_eta->clear();
  mu_phi->clear();

  jet_pt->clear();
  jet_eta->clear();
  jet_phi->clear();

  met_et  = 0.;
  met_phi = 0.;

  met_truth_et  = 0.;
  met_truth_phi = 0.;

  meff = 0;
  ht  = 0.;
  rt2 = -1.;
  rt4 = -1.;

  dphi_jetmet = -1.;
  dphi_gamjet = -1.;
  dphi_gammet = -1.;

  if (m_pdfrw) {
    weight_pdf1->clear();
    weight_pdf2->clear();
    weight_pdf3->clear();
  }

}

