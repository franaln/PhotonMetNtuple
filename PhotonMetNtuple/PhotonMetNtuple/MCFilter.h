#ifndef MCFilter_h
#define MCFilter_h

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "xAODTruth/TruthParticleContainer.h"

class MCFilter {

public:
  MCFilter() {};
  ~MCFilter() {};

  Bool_t accept_event(uint32_t did, xAOD::TEvent& event);
  Bool_t accept_vjets_event(xAOD::TEvent& event);
  Bool_t accept_ttbar_event(xAOD::TEvent& event);
  Bool_t accept_ttgamma_event(xAOD::TEvent& event);
  Bool_t accept_singletop_event(xAOD::TEvent& event);
  Bool_t accept_photonjet_event(xAOD::TEvent& event, uint32_t did);

  const xAOD::TruthParticle*  get_mother(const xAOD::TruthParticle* thePart);

};

#endif
