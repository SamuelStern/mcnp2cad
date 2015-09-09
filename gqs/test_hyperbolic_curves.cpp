

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
#define CHECK_REAL_EQUAL( EXP, ACT, EPS ) check_equal ( (EXP), (ACT), (EPS))

bool check_equal(int A, int B) { EQUAL_TEST( A == B, d) };
bool check_equal(double A, double B) { EQUAL_TEST( A == B, f) };
bool check_equal(double A, double B, double eps) { EQUAL_TEST( (A-B < eps), f ) };


double tol = 1e-6;

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
  CubitVector expected_midpoint(1.0,0.0,0.0);
  CubitVector actual_midpoint;
  
  edges[0]->mid_point(actual_midpoint);
  CHECK_REAL_EQUAL ( 1.0, actual_midpoint[0], tol );
  
  //should be created by now, time to export
  if (false)
    {
      DLIList<RefEntity*> exp_bodies;
      int exp_ents;
      CubitString cubit_version("12.2");
      CubitCompat_export_solid_model(exp_bodies, "test_file.sat", "ACIS_SAT", exp_ents, cubit_version);
    }



}


void check_CubitVector( CubitVector A, CubitVector B )
{

  CHECK_REAL_EQUAL( A[0], B[0], tol);
  CHECK_REAL_EQUAL( A[1], B[1], tol);
  CHECK_REAL_EQUAL( A[2], B[2], tol);


}
