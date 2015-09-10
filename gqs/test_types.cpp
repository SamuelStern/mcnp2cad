
#include "test_gq_utils.hpp"


int main()
{
  double A,B,C,D,E,F,G,H,J,K;
  
  //Start with the Parabolic Cylinder
  A = 1.0; B = 0.0; C = 0.0;
  D = 0.0; E = 0,0; F = 0,0; 
  G = 0.0; H = -1.0; J = 0.0;
  K = 0.0;


  GQ_TYPE this_type= characterize_surf(A,B,C,D,E,F,G,H,J,K);
  CHECK_EQUAL(PARABOLIC_CYL, this_type);

  DLIList<RefVolume*> vols;
  gqt->ref_volumes(vols);

  CHECK_EQUAL(1, (int)vols.size());
  
  //Elliptic Cylinder
  clean_cgma();
  A = 1.0; B = 2.0; C = 0.0;
  D = 0.0; E = 0,0; F = 0,0; 
  G = 0.0; H = 0.0; J = 0.0;
  K = -1.0;

  this_type= characterize_surf(A,B,C,D,E,F,G,H,J,K);
  CHECK_EQUAL(ELLIPTIC_CYL,this_type);

  
  //Hyperbolic Cylinder
  clean_cgma();
  A = 1.0; B = -2.0; C = 0.0;
  D = 0.0; E = 0,0; F = 0,0; 
  G = 0.0; H = 0.0; J = 0.0;
  K = -1.0;

  this_type= characterize_surf(A,B,C,D,E,F,G,H,J,K);
  CHECK_EQUAL(HYPERBOLIC_CYL,this_type);


  //Elliptic Paraboloid
  clean_cgma();
  A = 1.0; B = 2.0; C = 0.0;
  D = 0.0; E = 0,0; F = 0,0; 
  G = 0.0; H = 0.0; J = -1.0;
  K = 0.0;

  this_type= characterize_surf(A,B,C,D,E,F,G,H,J,K);
  CHECK_EQUAL(ELLIPTIC_PARABOLOID,this_type);

  //Elliptic Cone
  clean_cgma();
  A = 1.0; B = 2.0; C = -3.0;
  D = 0.0; E = 0,0; F = 0,0; 
  G = 0.0; H = 0.0; J = 0.0;
  K = 0.0;

  this_type= characterize_surf(A,B,C,D,E,F,G,H,J,K);
  CHECK_EQUAL(ELLIPTIC_CONE,this_type);

  //One Sheet Hyperpoloid
  clean_cgma();
  A = 1.0; B = 2.0; C = -3.0;
  D = 0.0; E = 0,0; F = 0,0; 
  G = 0.0; H = 0.0; J = 0.0;
  K = -1.0;

  this_type= characterize_surf(A,B,C,D,E,F,G,H,J,K);
  CHECK_EQUAL(ONE_SHEET_HYPERBOLOID,this_type);


  //Two Sheet Hyperpoloid
  clean_cgma();
  A = 1.0; B = -2.0; C = -3.0;
  D = 0.0; E = 0,0; F = 0,0; 
  G = 0.0; H = 0.0; J = 0.0;
  K = -1.0;

  this_type= characterize_surf(A,B,C,D,E,F,G,H,J,K);
  CHECK_EQUAL(TWO_SHEET_HYPERBOLOID,this_type);

  //Hyperbolic paraboloid will go here someday...


  
  //if we get to this point, all of the tests have passed
  std::cout << "PASSED" << std::endl;
  return 0;
}
