#include <iostream>



#include "iGeom.h"
#include "options.hpp"
#include "geometry.hpp"
#include "MCNPInput.hpp"
#include "volumes.hpp"
#include "testutil.hpp"

struct program_option_struct Gopt;

int main()
{


  //start by creating a Sample GQ and make sure we get the right type
  

  double A,B,C,D,E,F,G,H,J,K;

  A = 1.0; B = 1.0; C = 1.0;
  D = 0.0; E = 0.0; F = 0.0;
  G = 0.0; H = 0.0; J = 0.0;
  K = -1.0;

  GeneralQuadraticSurface GQ(A,B,C,D,E,F,G,H,J,K);

  CHECK_REAL_EQUAL(A,GQ.A,CHECK_TOLERANCE);
  CHECK_REAL_EQUAL(B,GQ.B,CHECK_TOLERANCE);
  CHECK_REAL_EQUAL(C,GQ.C,CHECK_TOLERANCE);
  CHECK_REAL_EQUAL(D,GQ.D,CHECK_TOLERANCE);
  CHECK_REAL_EQUAL(E,GQ.E,CHECK_TOLERANCE);
  CHECK_REAL_EQUAL(F,GQ.F,CHECK_TOLERANCE);
  CHECK_REAL_EQUAL(G,GQ.G,CHECK_TOLERANCE);
  CHECK_REAL_EQUAL(H,GQ.H,CHECK_TOLERANCE);
  CHECK_REAL_EQUAL(J,GQ.J,CHECK_TOLERANCE);
  CHECK_REAL_EQUAL(K,GQ.K,CHECK_TOLERANCE);


  return 0;

}
