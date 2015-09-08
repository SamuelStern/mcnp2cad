#include "GeometryModifyTool.hpp"
#include "GeometryQueryTool.hpp"
#include "GMem.hpp"
#include "RefVertex.hpp"
#include "RefEdge.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include "Surface.hpp"
#include "InitCGMA.hpp"
#include "CubitCompat.hpp"


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



void elliptic_cone(double a, double b, double c, double g, double h, double j, double k);
void elliptic_paraboloid(double a, double b, double c, double g, double h, double j, double k);
void elliptic_cyl(double a, double b, double c, double g, double h, double j, double k);
void parabolic_cyl(double a, double b, double c, double g, double h, double j, double k);
void hyperbolic_cyl(double a, double b, double c, double g, double h, double j, double k);
void one_sheet_hyperboloid(double a, double b, double c, double g, double h, double j, double k);
void two_sheet_hyperboloid(double a, double b, double c, double g, double h, double j, double k);

void hyperbolic_curves(double a, double b, DLIList<RefEdge*> &edge_list);
void hyperbolic_curves_in_plane( double a, double b, int ax1, int ax2, DLIList<RefEdge*> &edge_list);


GQ_TYPE characterize_surf( double A,
		       double B,
		       double C, 
		       double D, 
		       double E,
		       double F,
		       double G, 
		       double H, 
		       double J,
		       double K);

void complete_square ( double A,
		       double B,
		       double C, 
		       double D, 
		       double E,
		       double F,
		       double G, 
		       double H, 
		       double J,
		       double K,
		       double &a,
		       double &b, 
		       double &c, 
		       double &rhs,
		       double &W);


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
