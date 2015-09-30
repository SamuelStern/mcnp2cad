
#include <stdio.h>

#define CHECK_TOLERANCE 1e-5

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

void check_equal(int A, int B, int line) { EQUAL_TEST( A == B, d) };
void check_equal(double A, double B, int line) { EQUAL_TEST( A == B, f) };
void check_equal(double A, double B, double eps, int line) { EQUAL_TEST( (fabs(A-B) < eps), f ) };
