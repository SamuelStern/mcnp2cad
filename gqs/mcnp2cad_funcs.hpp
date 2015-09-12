#include "gq.hpp"

enum GQ_TYPE {UNKNOWN = 0,
	      ELLIPSOID,
	      ONE_SHEET_HYPERBOLOID,
	      TWO_SHEET_HYPERBOLOID,
	      ELLIPTIC_CONE,
	      ELLIPTIC_PARABOLOID,
	      HYPERBOLIC_PARABOLOID,
	      ELLIPTIC_CYL,
	      HYPERBOLIC_CYL,
	      PARABOLIC_CYL};


std::map<GQ_TYPE,void (*)(double,double,double,double,double,double,double)>  gq_funcs();
// Function for charaterizing the sub-type of generalized quadratic described by the input coefficients.
GQ_TYPE characterize_surf( double A,
		       double B,
		       double C, 
		       double G, 
		       double H, 
		       double J,
			   double K);


void complete_square ( double &A,
		       double &B,
		       double &C, 
		       double &G, 
		       double &H, 
		       double &J,
		       double &K,
		       double &dx,
		       double &dy,
		       double &dz);


void get_translation( double A,
		      double B,
		      double C, 
		      double D, 
		      double E,
		      double F,
		      double G, 
		      double H, 
		      double J,
		      double K,
		      double &dx,
		      double &dy, 
		      double &dz);
void get_rotation(double &A,
		  double &B,
		  double &C, 
		  double &D, 
		  double &E,
		  double &F,
		  double &alpha,
		  double &beta,
		  double &theta);
