
#include "ProgOptions.hpp"
#include "gq.hpp"
#include "mcnp2cad_funcs.hpp"

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

  std::map<GQ_TYPE,void (*)(double,double,double,double,double,double,double)> gqs = gq_funcs();
  

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


  //complete the square on the second and first order terms
  
  
  //The first step is to characterize the surface
  GQ_TYPE type = characterize_surf(A,B,C,D,E,F,G,H,J,K);
  gqs[type](A,B,C,G,H,J,K);
  
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
