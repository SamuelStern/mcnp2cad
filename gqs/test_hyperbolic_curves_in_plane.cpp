#include "test_gq_utils.hpp"

int main()
{
  
  //hyperbolic function coefficients
  double A = 1;
  double B = -2;
  int reflection_axis = 2;
  int symmetric_axis = 1;
  DLIList<RefEdge*> edges;
  //generate the desired curves
  hyperbolic_curves_in_plane(A,B,symmetric_axis,reflection_axis,edges);

  
  BasicTopologyEntity* bte1 = dynamic_cast<BasicTopologyEntity*>(edges[0]);
  BasicTopologyEntity* bte2 = dynamic_cast<BasicTopologyEntity*>(edges[1]);
    
  CubitBox edge1_box = bte1->bounding_box();
  CubitBox edge2_box = bte2->bounding_box();

  //check the bounds of each curve

  ////// Curve 1 Check ///////
  CubitVector expected_startpoint(0.0,1.4142135,2.0);
  CubitVector actual_startpoint = edges[0]->start_vertex()->center_point();
  CHECK_CUBITVECTORS_EQUAL(expected_startpoint, actual_startpoint);
  
  CubitVector expected_midpoint(0.0,1.0,0.0);
  CubitVector actual_midpoint;
  edges[0]->mid_point(actual_midpoint);
  CHECK_CUBITVECTORS_EQUAL(expected_midpoint, actual_midpoint );

  CubitVector expected_endpoint(0.0,1.4142135,-2.0);
  CubitVector actual_endpoint = edges[0]->end_vertex()->center_point();
  CHECK_CUBITVECTORS_EQUAL(expected_endpoint, actual_endpoint);

  ////// Curve 2 Check ///////
  expected_startpoint.set(0.0,-1.4142135,2.0);
  actual_startpoint = edges[1]->start_vertex()->center_point();
  CHECK_CUBITVECTORS_EQUAL(expected_startpoint, actual_startpoint);
  
  expected_midpoint.set(0.0,1.0,0.0);
  actual_midpoint;
  edges[0]->mid_point(actual_midpoint);
  CHECK_CUBITVECTORS_EQUAL(expected_midpoint, actual_midpoint );

  expected_endpoint.set(0.0,-1.4142135,-2.0);
  actual_endpoint = edges[1]->end_vertex()->center_point();
  CHECK_CUBITVECTORS_EQUAL(expected_endpoint, actual_endpoint);

  //if we get to here, the tests have passed
  std::cout << "PASSED" << std::endl;

}
