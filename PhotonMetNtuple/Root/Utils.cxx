#include <TMath.h>
#include <TVector2.h>

#include <PhotonMetNtuple/Utils.h>


double get_dphi(Double_t phi1, Double_t phi2)
{
  double  phi = fabs(phi1 - phi2);
  if (phi <= TMath::Pi())  return phi;
  else                     return (2 * TMath::Pi() - phi); 
}
