

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

#include "moab/ProgOptions.hpp"

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


std::ostream& operator<<(std::ostream& out, const GQ_TYPE value){
    static std::map<GQ_TYPE, std::string> strings;
    if (strings.size() == 0){
#define INSERT_ELEMENT(p) strings[p] = #p
      INSERT_ELEMENT(UNKNOWN);
      INSERT_ELEMENT(ELLIPSOID);
      INSERT_ELEMENT(ONE_SHEET_HYPERBOLOID);
      INSERT_ELEMENT(TWO_SHEET_HYPERBOLOID);
      INSERT_ELEMENT(ELLIPTIC_CONE);
      INSERT_ELEMENT(ELLIPTIC_PARABOLOID);
      INSERT_ELEMENT(HYPERBOLIC_PARABOLOID);
      INSERT_ELEMENT(ELLIPTIC_CYL);
      INSERT_ELEMENT(HYPERBOLIC_CYL);
      INSERT_ELEMENT(PARABOLIC_CYL);
#undef INSERT_ELEMENT
    }   
    return out << strings[value];
}


//doing this globally for now to make function signatures easier to write
CubitStatus stat = InitCGMA::initialize_cgma(); 

GeometryModifyTool *gmt = GeometryModifyTool::instance();
GeometryQueryTool *gqt = GeometryQueryTool::instance();


void elliptic_cone(double a, double b, double c, double g, double h, double j, double k);
void elliptic_paraboloid(double a, double b, double c, double k);


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

/*
 Desription: Program for the creation of a Generalized Quadratic (GQ) surface using CGM

 Input: 10 double values indicating the coefficients (A,B,C,D,E,F,G,H,J,K) to be used in defining the GQ equation as follows:

 Ax^2 + By^2 + Cz^2 + Dxy + Eyx + Fzx + Gz + Hy + Jz + K = 0

 Output: .sat file containing the described GQ surface.
*/
int main ( int argc, char** argv ) {

  double A,B,C,D,E,F,G,H,J,K;
  std::string filename = "GQ.sat";

  ProgOptions po( "GQ: A program for generating generalized quadratic surfaces in CGM.");

  po.addRequiredArg<double>("A", "x-squared coefficient", &A);
  po.addRequiredArg<double>("B", "y-squared coefficient", &B);
  po.addRequiredArg<double>("C", "z-squared coefficient", &C);
  po.addRequiredArg<double>("D", "xy coefficient", &D);
  po.addRequiredArg<double>("E", "yz coefficient", &E);
  po.addRequiredArg<double>("F", "xz coefficient", &F);
  po.addRequiredArg<double>("G", "x coefficient", &G);
  po.addRequiredArg<double>("H", "y coefficient", &H);
  po.addRequiredArg<double>("J", "z coefficient", &J);
  po.addRequiredArg<double>("K", "offset", &K);

  po.addOpt<std::string>("o", "Sat File", &filename);

  po.parseCommandLine( argc, argv );

  //make sure this is actually a quadratic surface
  if ( A == 0 && B == 0 && C == 0 )
    {
      std::cout << "All 2nd order coeffs are zero. This is not a GQ. Exiting..." << std::endl;
      return 1;
    }

  //let's rule out rotations (for now)
  if ( D != 0 || E != 0 || F !=0 )
    {
      std::cout << "Rotations are unsupported right now." << std::endl;
      return 1;
    }

  //The first step is to characterize the surface
  GQ_TYPE type = characterize_surf(A,B,C,D,E,F,G,H,J,K);

  if (!type)
    {
      std::cout << "This GQ type is not yet supported. Exiting..." << std::endl;
      return 1;
    }

  std::cout << "This GQ has type: " << type << std::endl; 

  //should be created by now, time to export
  DLIList<RefEntity*> exp_bodies;
  int exp_ents;
  CubitString cubit_version("12.2");
  
  CubitCompat_export_solid_model(exp_bodies, filename.c_str(), "ACIS_SAT", exp_ents, cubit_version);

  double dx, dy, dz; 
  
  get_translation(A,B,C,D,E,F,G,H,J,K,dx,dy,dz);

  return 0;

}


// Function for charaterizing the sub-type of generalized quadratic described by the input coefficients.
GQ_TYPE characterize_surf( double A,
		       double B,
		       double C, 
		       double D, 
		       double E,
		       double F,
		       double G, 
		       double H, 
		       double J,
		       double K)
{

  //first things first, complete the square to get the equation in the correct form
  
  double a,b,c,rhs,W;

  complete_square(A,B,C,D,E,F,G,H,J,K,a,b,c,rhs,W);
  

  std::cout << "a= " << a << std::endl;

  std::cout << "b= " << b << std::endl;
  
  std::cout << "c= " << c << std::endl;

  std::cout << "rhs= " << rhs << std::endl;

  //now start to determine what kind of surface we have 
  int num_zero = 0;
  int num_neg = 0;

  //count negative coefficients
  if ( a < 0 ) num_neg++;
  if ( b < 0 ) num_neg++;
  if ( c < 0 ) num_neg++;

  //count zero coefficients
  if ( a == 0) num_zero++;
  if ( b == 0) num_zero++;
  if ( c == 0) num_zero++;


  double g,h,j;
  g = G; h = H; j = J;
  if (W!=0)
    {
      g/=W;
      h/=W;
      j/=W;
    }

  if ( num_neg == 0 && num_zero == 0 && rhs)
    return ELLIPSOID;
  else if ( num_neg == 1 && num_zero == 0 && rhs )
    return ONE_SHEET_HYPERBOLOID;
  else if ( num_neg == 2 && num_zero == 0 && rhs )
    return TWO_SHEET_HYPERBOLOID;
  else if ( num_neg == 1 && num_zero == 0 && !rhs )
    {
      elliptic_cone( a, b, c, g, h, j, rhs );
      return ELLIPTIC_CONE;
    }
  else if ( num_neg == 0 && num_zero == 1 && !rhs )
    {
      elliptic_paraboloid( a, b, c, rhs );
      return ELLIPTIC_PARABOLOID;
    }
  else if ( num_neg == 1 && num_zero == 1 && !rhs )
    return HYPERBOLIC_PARABOLOID;
  else if ( num_neg == 0 && num_zero == 1 && rhs )
    return ELLIPTIC_CYL;
  else if ( num_neg == 1 && num_zero == 1 && rhs )
    return HYPERBOLIC_CYL;
  else if ( num_zero == 2 && !rhs )
    return PARABOLIC_CYL;

  
  return UNKNOWN;

}


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
		       double &W)
{

  W = -K;
  W += (A == 0) ? 0 : (G*G)/(4*A);
  W += (B == 0) ? 0 : (H*H)/(4*B);
  W += (C == 0) ? 0 : (J*J)/(4*C);


  double dum = ( W == 0 ) ? 1 : W;

  a = A/dum;

  b = B/dum;
  
  c = C/dum;

  rhs = ( W == 0 ) ? 0 : 1;

  return;

}

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
		  double &dz)
{

  dx = (A == 0) ? 0 : G/(2*A);
  dy = (B == 0) ? 0 : H/(2*B);
  dz = (C == 0) ? 0 : J/(2*C);

  if ( G/A < 0 ) dx *= -1;
  if ( H/B < 0 ) dy *= -1;
  if ( J/C < 0 ) dz *= -1;

  return;

}


void elliptic_cone( double a, double b, double c, double g, double h, double j, double k )
{

  //make sure that k is zero for this case
  if ( k != 0 ) 
    {
      std::cout << "Should not be here. Constant of the GQ is non-zero. Exiting..." << std::endl;
      exit(1);
    }

  bool pre_rot;
  double height,mnr,mjr;
  height = 0;

  int axis = 3;
  //figure out which of the coefficients is negative and set params accordingly
  if ( a < 0 ) 
    {
      height = sqrt(fabs(a));
      pre_rot = (b < c); // if the major axis is not along b, an extra rotation is needed
      mnr = pre_rot ? b : c;
      mjr = pre_rot ? c : b;
      axis = 1;
    }
  else if ( b < 0) 
    { 
      height = sqrt(fabs(b));
      pre_rot = (a < c); // if the major axes is not along a, an extra rotation is needed
      mnr = pre_rot ? a : c;
      mjr = pre_rot ? c : a;
      axis = 0;
    }
  else if ( c < 0 ) 
    {
      height = sqrt(fabs(c));
      pre_rot = (a < b); // if the minor axes it not along b, an extra rotation is needed
      mnr = pre_rot ? a : b;
      mjr = pre_rot ? b : a;
      axis = 2;
    }
 
  if ( 0 == height || 3 == axis ) 
    {
      std::cout << "Could not find a negtaive coefficient. Error. Exiting..." << std::endl;
      exit(1);
    }

  CubitStatus result;

  //fortunately there is a direct cgm function for this
  Body* ent = gmt->cylinder( height, mjr, mnr, 0);
  
  //if the minor/major radii need to be switched, do it now
  if (pre_rot)
    {
      result = gqt->rotate(ent,CubitVector(0,0,1),90);
      if (result != CUBIT_SUCCESS)
	std::cout << "Error rotating entity." << std::endl;
    }
  

  if ( 2 != axis ) //only re-orient the body if h is not along z-axis
    {
      double rot_axis[3] = {0,0,0};

      rot_axis[axis] = 1;

      double angle = (axis) ? 90 : -90;

      CubitVector ax(rot_axis[0],rot_axis[1],rot_axis[2]);

      std::cout << rot_axis[0] << rot_axis[1] << rot_axis[2]  << std::endl;

      result = gqt->rotate(ent,ax,angle);
      if ( result != CUBIT_SUCCESS )
	std::cout << "Error rotating entity." << std::endl;
    }


}

void elliptic_paraboloid( double a, double b, double c, double k )
{

  //figure out which direction is zero
  



}
