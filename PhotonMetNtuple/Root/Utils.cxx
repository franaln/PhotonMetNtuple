#include <TMath.h>
#include <TVector2.h>

#include <PhotonMetNtuple/Utils.h>


float get_dphi(float phi1, float phi2)
{
  float  phi = fabs(phi1 - phi2);
  if (phi <= TMath::Pi())  return phi;
  else                     return (2 * TMath::Pi() - phi); 
}

float get_deta(float eta1, float eta2)
{
  return fabs(eta1 - eta2);
}

float get_dr(float eta1, float phi1, float eta2, float phi2)
{
  float dphi = get_dphi(phi1, phi2);
  float deta = get_deta(eta1, eta2);

  return TMath::Hypot(deta, dphi);
}
