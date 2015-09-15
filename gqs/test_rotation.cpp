
#include "test_gq_utils.hpp"

int main()
{

  double A,B,C,D,E,F,G,H,J,K; //coefficients
  double x,y,z; //angles
  
  //TEST 1 - Parabolic Cylinder
  A = 1.0; B = 1.0; C = 1.0;
  D = 1.0; E = 0.0; F = 0.0;
  G = 0.0; H = 0.0; J = 0.0;
  K = 9.0;

  get_rotation(A,B,C,D,E,F,x,y,z);

  CHECK_REAL_EQUAL( 1.5, A, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 0.5, B, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 0.5, C, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 9.0, K, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 0.0, y, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 45.0, z, CHECK_TOLERANCE);

  ///TEST 2 - Ellipsoid
  A = 103.0; B = 125.0; C = 66.0;
  D = -48.0; E = -60.0; F = -12.0;
  G = 0.0; H = 0.0; J = 0.0;
  K = -294.0;

  get_rotation(A,B,C,D,E,F,x,y,z);

  CHECK_REAL_EQUAL( 147 , A, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 98, B, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 49, C, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( -294, K, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 26.56575, x, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 16.6015, y, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 63.434, z, CHECK_TOLERANCE);

  //TEST 3 - No rotation (for robustness)
  A = 1.0; B = 1.0; C = 1.0;
  D = 0.0; E = 0.0; F = 0.0;
  G = 0.0; H = 0.0; J = 0.0;
  K = 9.0;

  get_rotation(A,B,C,D,E,F,x,y,z);
  CHECK_REAL_EQUAL( 1.0, A, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 1.0, B, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 1.0, C, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 0.0, x, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 90.0, y, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 90.0, z, CHECK_TOLERANCE);
  
  
  //TEST 4 - Parabolic Cylinder
  A = 1.0; B = 1.0; C = 1.0;
  D = -1.0; E = 0.0; F = 0.0;
  G = 0.0; H = 0.0; J = 0.0;
  K = 9.0;

  get_rotation(A,B,C,D,E,F,x,y,z);

  CHECK_REAL_EQUAL( 1.5, A, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 0.5, B, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 0.5, C, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 9.0, K, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 0.0, x, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( 0.0, y, CHECK_TOLERANCE);
  CHECK_REAL_EQUAL( -45.0, z, CHECK_TOLERANCE);


  //if we get to this point, all of the tests have passed
  std::cout << "PASSED" << std::endl;
  return 0;
}
