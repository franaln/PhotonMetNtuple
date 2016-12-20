#include <TMath.h>
#include <TVector2.h>

#include <PhotonMetNtuple/Utils.h>


double get_dphi(Double_t phi1, Double_t phi2)
{
  return TVector2::Phi_mpi_pi(phi1 - phi2);
}
