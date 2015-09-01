

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
	      HYPERBOLOID,
	      CONE,
	      PARABOLOID,
	      CYLINDER};

int characterize_surf( double A,
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
		       double &rhs); 


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


  CubitStatus stat = InitCGMA::initialize_cgma(); 

  GeometryModifyTool *gmt = GeometryModifyTool::instance();

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

  //let's rule out rotations (for now)
  if ( D != 0 || E != 0 || F !=0 )
    {
      std::cout << "Rotations are unsupported right now." << std::endl;
      return 1;
    }

  //The first step is to characterize the surface
  characterize_surf(A,B,C,D,E,F,G,H,J,K);
  


  return 0;

}


// Function for charaterizing the sub-type of generalized quadratic described by the input coefficients.
int characterize_surf( double A,
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
  
  double a,b,c,rhs;

  complete_square( A,B,C,D,E,F,G,H,J,K,a,b,c,rhs);
  

  std::cout << "a= " << a << std::endl;

  std::cout << "b= " << b << std::endl;
  
  std::cout << "c= " << c << std::endl;

  std::cout << "rhs= " << rhs << std::endl;


  return 0;

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
		       double &rhs)
{

  double W = -K;
 
  W += (A == 0) ? 0 : (G*G)/(4*A);
  W += (B == 0) ? 0 : (H*H)/(4*B);
  W += (C == 0) ? 0 : (J*J)/(4*C);

  a = A/W;

  b = B/W;
  
  c = C/W;

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

  return;

}
