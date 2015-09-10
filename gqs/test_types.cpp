
#include "test_gq_utils.hpp"


int main()
{
  double A,B,C,D,E,F,G,H,J,K;

  A = 1.0; B = 0.0; C = 0.0;
  D = 0.0; E = 0,0; F = 0,0; 
  G = 0.0; H = -1.0; J = 0.0;
  K = 0.0;

  //Start with the Parabolic Cylinder
  GQ_TYPE this_type= characterize_surf(A,B,C,D,E,F,G,H,J,K);
  CHECK_EQUAL(PARABOLIC_CYL, this_type);

  DLIList<RefVolume*> vols;
  gqt->ref_volumes(vols);

  CHECK_EQUAL(1, (int)vols.size());

  //if we get to this point, all of the tests have passed
  std::cout << "PASSED" << std::endl;
  return 0;
}
