
#include "moab/ProgOptions.hpp"
#include "gq.hpp"

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
      elliptic_paraboloid( a, b, c, g, h, j, rhs );
      return ELLIPTIC_PARABOLOID;
    }
  else if ( num_neg == 1 && num_zero == 1 && !rhs )
    return HYPERBOLIC_PARABOLOID;
  else if ( num_neg == 0 && num_zero == 1 && rhs )
    {
      elliptic_cyl(a,b,c,g,h,j,rhs);
      return ELLIPTIC_CYL;
    }
  else if ( num_neg == 1 && num_zero == 1 && rhs )
    return HYPERBOLIC_CYL;
  else if ( num_zero == 2 && !rhs )
    {
      parabolic_cyl(a,b,c,g,h,j,rhs);
      return PARABOLIC_CYL;
    }

  
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

 void elliptic_paraboloid( double a, double b, double c, double g, double h, double j, double k )
{

  //figure out which direction is zero
 
  double r1,r2;
  int axis = 3;
  //figure out which direction is zero
  if (a == 0) 
    {
      axis = 0;
      r1 = b/g;
      r2 = c/g;
    }
  else if (b == 0)
    { 
      axis = 1;
      r1 = a/h;
      r2 = c/h;
    }
  else if (c == 0) 
    {
      axis = 2;
      r1 = a/j;
      r2 = b/j;
    }

if (3 == axis) 
  {
    std::cout << "Could not find a negtaive coefficient. Error. Exiting..." << std::endl;
    exit(1);
  }

//determine the size of the volume
 double height = 1e3; //some arbitrarily large height
 double mag = sqrt(fabs(height/r1)); //value of the trace at that height
 height *= (r1 < 0) ? -1 : 1; //sign adjustment
 

 //create a vertex on the trace of the parabola at this point
 double p1[3] = {0,0,0};
 p1[(a == 0)] = -mag; //<-- didn't know you could do this
 p1[axis] = height;
 RefVertex* v1 = gmt->make_RefVertex(CubitVector(p1[0],p1[1],p1[2]));
 if (!v1) std::cout << "Failed to create the first vertex." << std::endl;

 //create another vertex reflected across the rotation axis
 double pt2[3] = {0,0,0};
 pt2[(a == 0)] = mag; //<-- didn't know you could do this
 pt2[axis] = height;
 CubitVector pt2_pos(pt2[0],pt2[1],pt2[2]);
 RefVertex* v2 = gmt->make_RefVertex(CubitVector(pt2[0],pt2[1],pt2[2]));
 if (!v2) std::cout << "Failed to create the second vertex." << std::endl;

 //mid-point will always be the origin
 CubitVector mdpt(0,0,0);
 RefVertex* mv = gmt->make_RefVertex(mdpt);
 if(!mv) std::cout << "Could not create vertex at the origin." << std::endl;
 
 //we'll need a point above the origin at the height of the parabola
 //to create bounding curves for the surface of rotation
 double pt3_pos[3] = {0,0,0};
 pt3_pos[axis] = height;
 RefVertex* v3 = gmt->make_RefVertex(CubitVector(pt3_pos[0],pt3_pos[1],pt3_pos[2]));
 if(!v3) std::cout << "Error creating the bounding curve vertex." << std::endl;

 //create a vertiical line from the origin to the height
 RefEdge* line1 = gmt->make_RefEdge(STRAIGHT_CURVE_TYPE, v3, mv);
 if(!line1) std::cout << "Error creating the vertical bounding cure." << std::endl;

 //horizontal line to enclose the surface
 RefEdge* line2 = gmt->make_RefEdge(STRAIGHT_CURVE_TYPE, v3, v2);
 if(!line2) std::cout << "Error creating the horizontal bounding curve." << std::endl;

 //now create the parabolic curve 
 RefEdge* parab = gmt->make_RefEdge(PARABOLA_CURVE_TYPE, v1, v2, &mdpt);
 if (!parab) std::cout << "Failed to create the parabolic curve." << std::endl;

 //trim the curve at the origin
 CubitStatus result = gmt->trim_curve(parab,mdpt,pt2_pos);
 if ( result != CUBIT_SUCCESS ) std::cout << "Could not trim the prabolic curve." << std::endl;

 //gather the bounding curves to create the surface 
 // note: trimmed curve is not parab curve
 DLIList<RefEdge*> bounding_curves;
 gqt->ref_edges(bounding_curves);
 
 //create the surface
 RefFace* surf = gmt->make_RefFace(TORUS_SURFACE_TYPE, bounding_curves, true);
 if (!surf) std::cout << "Error creating the surface of rotation." << std::endl;

 //prepare surface entity for sweep
 RefEntity* ent = dynamic_cast<RefEntity*>(surf);
 DLIList<RefEntity*> surf_to_sweep;
 surf_to_sweep.insert(ent);
 DLIList<Body*> new_bodies;
 CubitVector sweep_point(0,0,0);
 double sweep_ax[3] = {0,0,0};
 sweep_ax[axis] = 1;
 CubitVector sweep_axis(sweep_ax[0],sweep_ax[1],sweep_ax[2]);
 
 //sweep the surface about the axis of rotation
 result = gmt->sweep_rotational(surf_to_sweep, sweep_point, sweep_axis, 2*CUBIT_PI, new_bodies, CUBIT_FALSE, CUBIT_FALSE);
 if ( result != CUBIT_SUCCESS ) std::cout << "Failed to create the swept entity." << std::endl;
 
 //need to get rid of the old surface reference
 gqt->delete_RefFace(surf);
 
 //make sure we got the correct number of volumes from the sweep
 assert( 1 == new_bodies.size() );
 
 //perpare to scale the volume along the final axis
 double scale_factor = r2/r1;
 CubitVector scale_pt(0,0,0); 
 std::cout << scale_factor << std::endl;
 double scale_factors[3] = {1,1,1};
 scale_factors[2-(c == 0)] = scale_factor;
 CubitVector factors(scale_factors[0],scale_factors[1],scale_factors[2]);
 
 //scale the volume 
 result = gmt->scale( new_bodies[0], scale_pt, factors);
 if (CUBIT_SUCCESS != result ) std::cout << "Error scaling the volume." << std::endl;

 //all done
 return;
 
}


void elliptic_cyl(double a, double b, double c, double g, double h, double j, double k)
{

  double r1,r2;
  int axis = 3;
  //figure out which direction is zero
  if (a == 0) 
    {
      axis = 0;
      r1 = c;
      r2 = b;
    }
  else if (b == 0)
    { 
      axis = 1;
      r1 = a;
      r2 = c;
    }
  else if (c == 0) 
    {
      axis = 2;
      r1 = a;
      r2 = b;
    }
  
  //use some arbitrary height for now
  double height = 10;

  //need a point on the origin to create the curve
  CubitVector origin(0,0,0);

  //create a surface from these points
  double plane[3] = {0,0,0};
  plane[axis] = 1;
  CubitVector gen_ax(plane[0],plane[1],plane[2]);

  CubitStatus result = gmt->create_ellipse_surface(r1,r2,gen_ax);
  
  //now get ready to sweep the ellipse along the generation axis
  DLIList<RefFace*> surfs;
  gqt->ref_faces(surfs);
  assert(1 == surfs.size());

  DLIList<RefEntity*> sweep_ents;
  RefEntity* ent = dynamic_cast<RefEntity*>(surfs[0]);
  sweep_ents.insert(ent);

  DLIList<Body*> new_bodies;

  result = gmt->sweep_translational(sweep_ents, height*gen_ax, 0, 0, CUBIT_FALSE, CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE, new_bodies);

  
  
  return;
}

void parabolic_cyl(double a, double b, double c, double g, double h, double j, double k)
{
  
  CubitVector ax1,ax2;
  double alpha; //equation of the form y = alpha*x^2

  //find the non-negative 2nd order coeff
  if (a != 0) 
    {
      const double ax[3] = {1,0,0};
      ax1.set(ax);
      alpha = a;
    }
  else if (b != 0)
    { 
      const double ax[3] = {0,1,0};
      ax1.set(ax);
      alpha = b;
    }
  else if (c != 0) 
    {
      const double ax[3] = {0,0,1};
      ax1.set(ax);
      alpha = c;
    }
  
  if (g != 0) 
    {
      const double ax[3] = {1,0,0};
      ax2.set(ax);
      alpha/=g;
    }
  else if (h != 0)
    { 
      const double ax[3] = {0,1,0};
      ax2.set(ax);
      alpha/=h;
    }
  else if (j != 0) 
    {
      const double ax[3] = {0,0,1};
      ax2.set(ax);
      alpha/=j;
    }

  CubitVector gen_ax = ax1*ax2;

  //get value along ax1 at height
  double height = 10;
  double mag = sqrt(fabs(height/alpha));
  DLIList<RefEdge*> bounding_curves;

  //create a trace of the parabola
  CubitVector p1 = mag*ax1 + height*ax2;
  RefVertex* v1 = gmt->make_RefVertex(p1);

  CubitVector p2 = -mag*ax1 + height*ax2;
  RefVertex* v2 = gmt->make_RefVertex(p2);

  CubitVector origin(0,0,0);

  RefEdge* parab = gmt->make_RefEdge(PARABOLA_CURVE_TYPE,v1,v2,&origin);
  bounding_curves.insert(parab);

  RefEdge* line = gmt->make_RefEdge(STRAIGHT_CURVE_TYPE, v1, v2);
  bounding_curves.insert(line);

  RefFace* trace = gmt->make_RefFace(TORUS_SURFACE_TYPE,bounding_curves,true);


  RefEntity* ent = dynamic_cast<RefEntity*>(trace);

  DLIList<RefEntity*> sweep_ents;
  sweep_ents.insert(ent);

  DLIList<Body*> new_bodies;
  

  double length = 10; //arbitrary length for now
  CubitStatus result = gmt->sweep_translational(sweep_ents, length*gen_ax, 0, 0, CUBIT_FALSE, CUBIT_FALSE, CUBIT_TRUE, CUBIT_FALSE, new_bodies);
  
  return;
}

