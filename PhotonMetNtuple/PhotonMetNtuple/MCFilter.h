#ifndef MCFilter_h
#define MCFilter_h


class MCFilter {

public:
  MCFilter() {};
  ~MCFilter() {};

  Bool_t accept_event(uint32_t did, xAOD::TEvent& event);
  Bool_t accept_vjets_event(xAOD::TEvent& event);
  Bool_t accept_ttbar_event(xAOD::TEvent& event);
  Bool_t accept_ttgamma_event(xAOD::TEvent& event);

  const xAOD::TruthParticle*  get_mother(const xAOD::TruthParticle* thePart);

};

#endif
