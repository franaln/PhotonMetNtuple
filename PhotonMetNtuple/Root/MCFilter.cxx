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
  // W/Z + jets
  if ((did >= 361300 && did <= 361371) || (did >= 361374 && did <= 361467)) {
    return accept_vjets_event(event);
  }

  // ttbar
  else if (did == 410000) {
    return accept_ttbar_event(event);
  }

  // ttgamma
  else if (did >= 410082 && did <= 410084) {
    return accept_ttgamma_event(event);
  }


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

Bool_t MCFilter::accept_ttbar_event(xAOD::TEvent& event)
{
  // veto events with at least one ME photon with pt>80 GeV

  const xAOD::TruthParticleContainer* particles = 0;
  RETURN_CHECK("MCFilter", event.retrieve(particles, "TruthParticles"));

  Bool_t accept = true;

  for (auto ip = particles->begin(); ip!=particles->end(); ++ip) {
    
    if ((*ip)->pdgId() != 22 || (*ip)->status() != 1) 
      continue;
    
    Int_t mother_pdg = get_mother(*ip)->absPdgId();
    if (mother_pdg < 100 && mother_pdg != 11) { // ME photon
      
      if ((*ip)->pt() > 80000.) {
        accept = false;
        break;
      }
    }

  }
  
  return accept;
}

Bool_t MCFilter::accept_ttgamma_event(xAOD::TEvent& event)
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
