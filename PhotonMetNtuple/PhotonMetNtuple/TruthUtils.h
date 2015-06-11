#ifndef TruthUtils_h
#define TruthUtils_h

#include "xAODJet/JetContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include "AthContainers/ConstDataVector.h"
#include <TLorentzVector.h>


class TruthParticle : public TLorentzVector {

 public:
  TruthParticle(const xAOD::TruthParticle *ptcl);
  TruthParticle(const xAOD::Jet *ptcl);

  Int_t index;
  Bool_t good;

  inline bool operator > (const TruthParticle &other) const
  {
    if (Pt() > other.Pt()) return true;
    return false;
  }
  inline bool operator < (const TruthParticle &other) const
  {
    if (Pt() < other.Pt()) return true;
    return false;
  }

};


/// print details about the truth particle to the screen
//void print_truth_ptcl(const xAOD::TruthParticle *ptcl, TString comment="", 
//                      int NchildDepth=0, int NparentDepth=0, int currentDepth=0);
  
Bool_t is_stable(const xAOD::TruthParticle*);

Int_t get_photon_type(const xAOD::TruthParticle *ptcl);

#endif
