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

 gq_map[ONE_SHEET_HYPERBOLOID]= &one_sheet_hyperboloid;
 gq_map[TWO_SHEET_HYPERBOLOID]= &two_sheet_hyperboloid;
 gq_map[ ELLIPTIC_CONE]= &elliptic_cone;
 gq_map[ ELLIPTIC_PARABOLOID]= &elliptic_paraboloid;
 gq_map[ ELLIPTIC_CYL]= &elliptic_cyl;
 gq_map[ HYPERBOLIC_CYL]= &hyperbolic_cyl;
 gq_map[ PARABOLIC_CYL]= &parabolic_cyl;

  
  return gq_map;
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
  
  if (false)
    {
      std::cout << "a= " << a << std::endl;
      
      std::cout << "b= " << b << std::endl;
      
      std::cout << "c= " << c << std::endl;
      
      std::cout << "rhs= " << rhs << std::endl;
    }

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
    {
      //we already have a function for this
      return ELLIPSOID;
    }
  else if ( num_neg == 1 && num_zero == 0 && rhs )
    {
      one_sheet_hyperboloid( a, b, c, g, h, j, rhs );
      return ONE_SHEET_HYPERBOLOID;
    }
  else if ( num_neg == 2 && num_zero == 0 && rhs )
    {
      two_sheet_hyperboloid( a, b, c, g, h, j, rhs );
      return TWO_SHEET_HYPERBOLOID;
    }
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
    {
      hyperbolic_cyl(a,b,c,g,h,j,rhs);
      return HYPERBOLIC_CYL;
    }
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
