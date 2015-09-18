#include "moab/CartVect.hpp"
#include "moab/Matrix3.hpp"
#include "mcnp2cad_funcs.hpp"


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


std::map<GQ_TYPE,void (*)(double,double,double,double,double,double,double)>  gq_funcs()
{
std::map<GQ_TYPE, void (*)(double,double,double,double,double,double,double)> gq_map;
 gq_map[UNKNOWN] = &not_supported;
 gq_map[ELLIPSOID] = &not_supported;
 gq_map[ONE_SHEET_HYPERBOLOID]= &one_sheet_hyperboloid;
 gq_map[TWO_SHEET_HYPERBOLOID]= &two_sheet_hyperboloid;
 gq_map[ELLIPTIC_CONE]= &elliptic_cone;
 gq_map[ELLIPTIC_PARABOLOID]= &elliptic_paraboloid;
 gq_map[ELLIPTIC_CYL]= &elliptic_cyl;
 gq_map[HYPERBOLIC_CYL]= &hyperbolic_cyl;
 gq_map[PARABOLIC_CYL]= &parabolic_cyl;

  
  return gq_map;
  }

// Function for charaterizing the sub-type of generalized quadratic described by the input coefficients.
GQ_TYPE characterize_surf( double A,
		       double B,
		       double C, 
		       double G, 
		       double H, 
		       double J,
		       double K)
{

  if (false)
    {
      std::cout << "a= " << A << std::endl;
      
      std::cout << "b= " << B << std::endl;
      
      std::cout << "c= " << C << std::endl;
      
      std::cout << "rhs= " << -K << std::endl;
    }

  //now start to determine what kind of surface we have 
  int num_zero = 0;
  int num_neg = 0;

  //count negative coefficients
  if ( A < 0 ) num_neg++;
  if ( B < 0 ) num_neg++;
  if ( C < 0 ) num_neg++;

  //count zero coefficients
  if ( A == 0) num_zero++;
  if ( B == 0) num_zero++;
  if ( C == 0) num_zero++;

  int rhs = (int)-K;
  
  if ( num_neg == 0 && num_zero == 0 && rhs)
    {
      //we already have a function for this
      return ELLIPSOID;
    }
  else if ( num_neg == 1 && num_zero == 0 && rhs )
    {
      //one_sheet_hyperboloid( a, b, c, g, h, j, rhs );
      return ONE_SHEET_HYPERBOLOID;
    }
  else if ( num_neg == 2 && num_zero == 0 && rhs )
    {
      //two_sheet_hyperboloid( a, b, c, g, h, j, rhs );
      return TWO_SHEET_HYPERBOLOID;
    }
  else if ( num_neg == 1 && num_zero == 0 && !rhs )
    {
      //elliptic_cone( a, b, c, g, h, j, rhs );
      return ELLIPTIC_CONE;
    }
  else if ( num_neg == 0 && num_zero == 1 && !rhs )
    {
      //elliptic_paraboloid( a, b, c, g, h, j, rhs );
      return ELLIPTIC_PARABOLOID;
    }
  else if ( num_neg == 1 && num_zero == 1 && !rhs )
    return HYPERBOLIC_PARABOLOID;
  else if ( num_neg == 0 && num_zero == 1 && rhs )
    {
      //elliptic_cyl(a,b,c,g,h,j,rhs);
      return ELLIPTIC_CYL;
    }
  else if ( num_neg == 1 && num_zero == 1 && rhs )
    {
      //hyperbolic_cyl(a,b,c,g,h,j,rhs);
      return HYPERBOLIC_CYL;
    }
  else if ( num_zero == 2 && !rhs )
    {
      //parabolic_cyl(a,b,c,g,h,j,rhs);
      return PARABOLIC_CYL;
    }

  
  return UNKNOWN;

}


void complete_square ( double &A,
		       double &B,
		       double &C, 
		       double &G, 
		       double &H, 
		       double &J,
		       double &K,
		       double &dx,
		       double &dy,
		       double &dz)
{

  double W;
  
  W = -K;
  W += (A == 0) ? 0 : (G*G)/(4*A);
  W += (B == 0) ? 0 : (H*H)/(4*B);
  W += (C == 0) ? 0 : (J*J)/(4*C);


  double dum = ( W == 0 ) ? 1 : W;

  A/=dum;

  B/=dum;
  
  C/=dum;

  G/=dum;

  H/=dum;

  J/=dum;

  K = ( W == 0 ) ? 0 : -1;

  dx = (A == 0) ? 0 : G/(2*A);
  dy = (B == 0) ? 0 : H/(2*B);
  dz = (C == 0) ? 0 : J/(2*C);

  if ( G/A < 0 ) dx *= -1;
  if ( H/B < 0 ) dy *= -1;
  if ( J/C < 0 ) dz *= -1;

  return;

}



void get_rotation(double &A,
		  double &B,
		  double &C, 
		  double &D, 
		  double &E,
		  double &F,
		  double &alpha,
		  double &theta,
		  double &phi)
{

  //reset angles
  alpha = 0.0;
  theta = 0.0;
  phi = 0.0;

  moab::Matrix3 coeff_mat(A,D/2,F/2,
			  D/2,B,E/2,
			  F/2,E/2,C);

  double eigen_vals[3];

  moab::CartVect eigen_vects[3];

  moab::Matrix::EigenDecomp(coeff_mat,eigen_vals,eigen_vects);

  moab::CartVect x_ax(1,0,0);
  moab::CartVect y_ax(0,1,0);
  moab::CartVect z_ax(0,0,1);
  
  //make sure we have unit vectors
  eigen_vects[0].normalize();
  eigen_vects[1].normalize();
  eigen_vects[2].normalize();


  moab::Matrix3 P( eigen_vects[0], eigen_vects[1], eigen_vects[2] );

  if ( fabs(P.determinant()-1.0) > 1e-6 ) //make sure we have a right-handed system
    {

      moab::Matrix3 new_P( eigen_vects[0], eigen_vects[2], eigen_vects[1]);
      
      //if for some reason we can't achieve a right-handed system, exit
      if ( new_P.determinant()-1.0 > 1e-6 )
	{
	  std::cout << "Could not orient new axes properly" << std::endl;
	  exit(1);

	}
      //update Eigenvector matrix
      P = new_P;
     
      //swap Eigenvalues
      double temp = eigen_vals[2];
      eigen_vals[2] = eigen_vals[1];
      eigen_vals[1] = temp;
      
    }

  A = eigen_vals[0];
  B = eigen_vals[1];
  C = eigen_vals[2];
  D = 0.0;
  E = 0.0; 
  F = 0.0;

  //calculate angles of rotation
  P = P.transpose(); //transpose P to get the correct conversion


  //attempt the simple solution first
  theta = -asin(P[2][0]);
  if ( fabs(cos(theta)) > 1.0e-8 )
    {
      alpha = atan2(P[2][1]/cos(theta),P[2][2]/cos(theta));
      phi = atan2(P[1][0]/cos(theta),P[0][0]/cos(theta));
    }
  else
    {
      phi = 0; //arbitrary value
      theta = ( P[1][1] == P[0][2] ) ? CUBIT_PI/2 : -CUBIT_PI/2;
      alpha = ( P[1][1] == P[0][2] ) ? atan2(P[0][1],P[0][2]) : atan2(-P[0][1],-P[0][2]);
    }

  //convert to degrees
  alpha*=180/CUBIT_PI;
  theta*=180/CUBIT_PI;
  phi *= 180/CUBIT_PI;

  return;
}
