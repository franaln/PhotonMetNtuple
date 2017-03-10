#ifndef Common_H
#define Common_H

#include "xAODBase/IParticleHelpers.h"

static const char *APP_NAME = "PhotonMetNtuple";
static const char *APP_VERSION = "v53";

static bool ptsorter(const xAOD::IParticle* j1, const xAOD::IParticle* j2) {
  return (j1->pt() > j2->pt());
}

#endif
