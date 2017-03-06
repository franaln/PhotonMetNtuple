#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"

#include "FourMomUtils/xAODP4Helpers.h"

#include "PhotonMetNtuple/MCFilter.h"

Bool_t MCFilter::accept_event(uint32_t did, xAOD::TEvent& event) 
{
  // // W/Z + jets
  // if ((did >= 361300 && did <= 361371) || (did >= 361374 && did <= 361467)) {
  //   return accept_vjets_event(event);
  // }

  // ttbar
  if (did == 410000) {
    return accept_ttbar_event(event);
  }

  // ttgamma
  else if ((did >= 410082 && did <= 410084) || did == 407320) {
    return accept_ttgamma_event(event);
  }

  // // singletop
  // else if (did == 410011 || did == 410012 || did == 410025 || did == 410026 || did == 410013 || did == 410014) {
  //   return accept_singletop_event(event);
  // }

  // else if (did >= 361042 && did <= 361061) {
  //   return accept_photonjet_event(event, did);
  // }

  return true;
}

Bool_t MCFilter::accept_vjets_event(xAOD::TEvent& event)
{
  /*
     Wjet/Wgamma
     veto events with at least one truth photon with pt>35 GeV y DR(photon,lepton)>0.1
     Zjet/Zgamma
     veto events with at least one truth photon with pt>35 GeV y DR(photon,lepton)>0.1
  */

  const xAOD::TruthParticleContainer* particles = 0;
  RETURN_CHECK("MCFilter", event.retrieve(particles, "TruthParticles"));


  Bool_t accept = true;
  for (auto ip = particles->begin(); ip!=particles->end(); ++ip) {

    Int_t pid_i = (*ip)->pdgId();
    Int_t st_i = (*ip)->status();
    
    Double_t pt_i = (*ip)->pt() * 0.001;
    
    if (pid_i == 22 && st_i == 1 && pt_i > 35) {
      
      for (auto jp = particles->begin(); jp!=particles->end(); ++jp) {
        
        Int_t abspid_j = abs((*jp)->pdgId());
        Int_t st_j = abs((*jp)->status());
        
        if ((abspid_j == 11 || abspid_j == 13 || abspid_j == 15 || abspid_j < 5)  && st_j == 3) {
          
          if (xAOD::P4Helpers::deltaR2(*ip, *jp) > 0.1 * 0.1) {
            accept = false;
            break;
          }
        }
      }
    }
  }
  
  return accept;
}

// ttbar/ttgamma
const xAOD::TruthParticle* MCFilter::get_last_truth_particle(const xAOD::TruthParticle* particle){ 
  auto pdgId = particle->pdgId();
  
  if (particle->nChildren() >= 1 && particle->child(0) &&  particle->child(0)->pdgId() == pdgId) //allow for t->t gluon
    return get_last_truth_particle(particle->child(0));
  
  return particle;
}

Bool_t MCFilter::has_me_photon(xAOD::TEvent& event)
{
  const xAOD::TruthParticleContainer* particles = 0;
  RETURN_CHECK("MCFilter", event.retrieve(particles, "TruthParticles"));

  const xAOD::IParticle *par = nullptr;
  int me_photon = 0;
  int pdgId = 22; // look at photons
  for (const auto& particle: *particles) {
    if (fabs(particle->pdgId()) == pdgId  && (particle->nParents()==0  || fabs(particle->parent(0)->pdgId()) != pdgId)) {// this particle is a photon
      par = get_last_truth_particle(particle);
      int motherPdgId = -1;
      if (particle->nParents() > 0) motherPdgId = particle->parent(0)->pdgId();
      if (abs(motherPdgId)<100 && particle->barcode() <2e5 && par->p4().Pt()>80e3)
        me_photon += 1;
    }
  }
  
  return (me_photon>0);
}

Bool_t MCFilter::accept_ttbar_event(xAOD::TEvent& event)
{
  // veto events with at least one ME photon with pt>80 GeV
  return !has_me_photon(event);
}

Bool_t MCFilter::accept_ttgamma_event(xAOD::TEvent& event)
{
  // accept events with at least one ME photon with pt>80 GeV
  return has_me_photon(event);
}

// single top
Bool_t MCFilter::accept_singletop_event(xAOD::TEvent& event)
{
  // accept events with at least one ME photon with pt>80 GeV

  const xAOD::TruthParticleContainer* particles = 0;
  RETURN_CHECK("MCFilter", event.retrieve(particles, "TruthParticles"));

  Bool_t accept = false;

  for (auto ip = particles->begin(); ip!=particles->end(); ++ip) {
    
    if ((*ip)->pdgId() == 22 && (*ip)->status() == 1 && (*ip)->pt() > 80000.) {

      // check if the the photon is ME
      Int_t mother_pdg = get_mother(*ip)->absPdgId();
      if (mother_pdg < 100 && mother_pdg != 11)  {
        accept = true;
        break;
      }
    }
  }
  
  return accept;
}


Bool_t MCFilter::accept_photonjet_event(xAOD::TEvent& event, uint32_t did)
{
  const xAOD::TruthParticleContainer* particles = 0;
  RETURN_CHECK("MCFilter", event.retrieve(particles, "TruthParticles"));

  Bool_t accept = false;

  Float_t lead_pt = 0.;
  for (auto ip = particles->begin(); ip!=particles->end(); ++ip) {
    if ((*ip)->pdgId() == 22 && (*ip)->status() == 1 && (*ip)->pt() > 80000. && (*ip)->pt() > lead_pt) {
      lead_pt = (*ip)->pt();
    }
  }

  if (did == 361042 || did == 361043 || did == 361044) {
    if (lead_pt >= 70000. && lead_pt < 140000.)
      accept = true;
  }
  else if (did == 361045 || did == 361046 || did == 361047) {
    if(lead_pt >= 140000. && lead_pt < 280000.)
      accept = true;
  }
  // else if (ph280sh_id == reader->mc_channel_number) { 
  //   if(truth_pt >= 350 * GeV && truth_pt < 600 * GeV)
  //     accept = true;
  // }
  // else if (ph500sh_id == reader->mc_channel_number) { 
  //   if(truth_pt >= 600 * GeV && truth_pt < 950 * GeV)
  //     accept = true;
  // }
  // else if (ph800sh_id == reader->mc_channel_number) { 
  //   if(truth_pt >= 950 * GeV && truth_pt < 1150 * GeV)
  //     accept = true;
  // }
  // else if (ph1000sh_id == reader->mc_channel_number) { 
  //   if (truth_pt >= 1150 * GeV)
  //     accept = true;
  // }

  if (!accept)
    std::cout << lead_pt << std::endl;
  
  return accept;
}


const xAOD::TruthParticle*  MCFilter::get_mother(const xAOD::TruthParticle* thePart)
{
  const xAOD::TruthVertex* partOriVert = thePart->prodVtx();

  int m_barcodeShift = 1000000;
  int m_barcodeG4Shift = 200001;
  long partPDG       = thePart->pdgId();
  int  partBarcode   = thePart->barcode()%m_barcodeShift;
  long MotherPDG(0);

  const xAOD::TruthVertex*   MothOriVert(0);
  const xAOD::TruthParticle* theMoth(0);

  if(!partOriVert) return theMoth;
  
  int itr=0;
  do {
    if (itr!=0) 
      partOriVert = MothOriVert;
    for (unsigned int ipIn=0;ipIn<partOriVert->nIncomingParticles();ipIn++) {
      theMoth = partOriVert->incomingParticle(ipIn);
      MotherPDG   = theMoth->pdgId();
      MothOriVert = theMoth->prodVtx();
      if (MotherPDG == partPDG) break;
    }
    itr++;
    if (itr>100) { 
      std::cout << "getMother:: infinite while" << std::endl;
      break;
    }
  }  while (MothOriVert !=0 && MotherPDG == partPDG && partBarcode<m_barcodeG4Shift);


  return theMoth;
}

