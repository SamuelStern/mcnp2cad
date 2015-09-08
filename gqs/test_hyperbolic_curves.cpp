

#include "gq.hpp"

int main()
{

  
  //hyperbolic function coefficients
  double A = 1;
  double B = -2;

  DLIList<RefEdge*> edges;
  //generate the desired curves
  hyperbolic_curves(A,B,edges);

  
  BasicTopologyEntity* bte = dynamic_cast<BasicTopologyEntity*>(edges[0]);

  CubitBox edge_box = bte->bounding_box();

  //check the bounds of each curve

  ////// Curve 1 Check ///////
  
  




}
