
#include "gq.hpp"


#define CHECK_TOLERANCE 1e-6

#define FLAG_ERROR exit(1)

inline void flag_error() { FLAG_ERROR; }

#define EQUAL_TEST( TEST, TYPE ) if( !(TEST) ) {	\
    printf( "Equality Test Failed at line %d: \n", line);	\
  printf( "  Expected value: %" #TYPE "\n", A ); \
  printf( "  Actual value:   %" #TYPE "\n", B ); \
  printf( "\n" ); \
  flag_error();	  \
  }

#define CHECK_EQUAL( EXP, ACT ) check_equal ( (EXP), (ACT), __LINE__ ) 
#define CHECK_REAL_EQUAL( EXP, ACT, EPS ) check_equal ( (EXP), (ACT), (EPS), __LINE__ )
#define CHECK_CUBITVECTORS_EQUAL( EXP, ACT ) check_CubitVectors_equal( (EXP), (ACT), #EXP, #ACT, __LINE__ )

bool check_equal(int A, int B, int line) { EQUAL_TEST( A == B, d) };
bool check_equal(double A, double B, int line) { EQUAL_TEST( A == B, f) };
bool check_equal(double A, double B, double eps, int line) { EQUAL_TEST( (A-B < eps), f ) };
bool check_CubitVectors_equal( CubitVector A, CubitVector B, const char* Aname, const char* Bname, int line ) 
{
if ( (A[0] - B[0] > CHECK_TOLERANCE) || ( A[1] - B[1] > CHECK_TOLERANCE ) || (A[2] - B[2] > CHECK_TOLERANCE ) )
  {

std::cout << "=====" << std::endl << "ERROR" << std::endl << "=====" << std::endl;
std::cout << "CubitVector " << Aname << " does not equal " << Bname
          << " (from line " << line << ")" << std::endl;
CHECK_REAL_EQUAL( A[0], B[0], CHECK_TOLERANCE);
CHECK_REAL_EQUAL( A[1], B[1], CHECK_TOLERANCE);
CHECK_REAL_EQUAL( A[2], B[2], CHECK_TOLERANCE);

}
}
