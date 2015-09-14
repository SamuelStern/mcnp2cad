
#include "test_gq_utils.hpp"

int main()
{
  double A,B,C,D,E,F,G,H,J,K;
  
  //Start with the Parabolic Cylinder
  A = 1.0; B = 1.0; C = 1.0;
  D = 1.0; E = 0.0; F = 0.0;
  G = 0.0; H = 0.0; J = 0.0;
  K = 9.0;

  double y,z;
  get_rotation(A,B,C,D,E,F,y,z);

  CHECK_REAL_EQUAL( 0, y, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 45, z, CHECK_TOLERANCE);

  //if we get to this point, all of the tests have passed
  std::cout << "PASSED" << std::endl;
  return 0;
}
