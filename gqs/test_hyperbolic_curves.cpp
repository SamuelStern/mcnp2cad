

#include "gq.hpp"


#define FLAG_ERROR throw 1

inline void flag_error() { FLAG_ERROR; }

#define EQUAL_TEST( TEST, TYPE ) if( !(TEST) ) {	\
  printf( "Equality Test Failed: \n"); \
  printf( "  Expected value: %" #TYPE "\n", A ); \
  printf( "  Actual value:   %" #TYPE "\n", B ); \
  printf( "\n" ); \
  flag_error();	  \
  }

#define CHECK_EQUAL( EXP, ACT ) check_equal ( (EXP), (ACT) ) 


bool check_equal(int A, int B) { EQUAL_TEST( A == B, d) };
bool check_equal(double A, double B) { EQUAL_TEST( A == B, f) };



int main()
{

  
  //hyperbolic function coefficients
  double A = 1;
  double B = -2;

  DLIList<RefEdge*> edges;
  //generate the desired curves
  hyperbolic_curves(A,B,edges);

  
  BasicTopologyEntity* bte1 = dynamic_cast<BasicTopologyEntity*>(edges[0]);
  BasicTopologyEntity* bte2 = dynamic_cast<BasicTopologyEntity*>(edges[1]);
    
  CubitBox edge1_box = bte1->bounding_box();
  CubitBox edge2_box = bte2->bounding_box();
  

  //check the bounds of each curve

  ////// Curve 1 Check ///////
  // CHECK_EQUAL( 2.0 , edge1_box.max_x() );
  // CHECK_EQUAL( 1.0 , edge1_box.min_x() );

  //should be created by now, time to export
  DLIList<RefEntity*> exp_bodies;
  int exp_ents;
  CubitString cubit_version("12.2");
  
  CubitCompat_export_solid_model(exp_bodies, "test_file.sat", "ACIS_SAT", exp_ents, cubit_version);




}
