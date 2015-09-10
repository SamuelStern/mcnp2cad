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

//doing this globally for now to make function signatures easier to write
extern GeometryModifyTool *gmt;
extern GeometryQueryTool *gqt;


void elliptic_cone(double a, double b, double c, double g, double h, double j, double k);
void elliptic_paraboloid(double a, double b, double c, double g, double h, double j, double k);
void elliptic_cyl(double a, double b, double c, double g, double h, double j, double k);
void parabolic_cyl(double a, double b, double c, double g, double h, double j, double k);
void hyperbolic_cyl(double a, double b, double c, double g, double h, double j, double k);
void one_sheet_hyperboloid(double a, double b, double c, double g, double h, double j, double k);
void two_sheet_hyperboloid(double a, double b, double c, double g, double h, double j, double k);



/* Creates two hyperbolic curves in the xy plane using the parameters a & b

Hyperbolic Curve Form:

x^2/a - y^2/b = 1

symmetric axis - x
reflection axis - y

Returns: two RefEdge pointers to the curves

*/
void hyperbolic_curves(double a, double b, DLIList<RefEdge*> &edge_list);

/* this function will return hyperbolic curves in a plane of two principle axes

symmetric_axis - intersecting axis of symmetry for one of the curves

reflecting_axis - reflecting axis for the curves

(for the axes arguments: 0 is x, 1 is y, 2 is z)

*/

void hyperbolic_curves_in_plane( double a, double b, int symmetric_axis, int reflecting_axis, DLIList<RefEdge*> &edge_list);


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
