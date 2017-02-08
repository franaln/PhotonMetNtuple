#ifndef TruthUtils_h
#define TruthUtils_h

#include "AthContainers/ConstDataVector.h"
#include "xAODJet/JetContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include <TLorentzVector.h>

class TruthParticle : public TLorentzVector {

 public:
  TruthParticle(const xAOD::TruthParticle *ptcl);
  TruthParticle(const xAOD::Jet *ptcl, Bool_t isbjet);

  Int_t index;
  Bool_t good;
  Bool_t isbjet;

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

namespace TruthUtils {
  Double_t GetDeltaPhi(Double_t phi1, Double_t phi2);
  Double_t GetDeltaR(const TruthParticle &p, const TruthParticle &q);
  Bool_t OverlapsOthers(TruthParticle &p, std::vector<TruthParticle> &others, Double_t deltaR_cut);
  void CleanBads(std::vector<TruthParticle> &vector);
}

  
Bool_t is_stable(const xAOD::TruthParticle*);

Int_t get_photon_type(const xAOD::TruthParticle *ptcl);

#endif
