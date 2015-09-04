
#include "gq.hpp"



void hyperbolic_curve();

int main()
{

  hyperbolic_curve();

  return 0;
}


void hyperbolic_curve()
{

  double a = 3;
  double b = 5;

  //first create a conic surface
  CubitVector p1(0,0,b);
  CubitVector p2(0,a,0);

  RefVertex* v1 = gmt->make_RefVertex(p1);

  RefVertex* v2 = gmt->make_RefVertex(p2);

  RefEdge* line = gmt->make_RefEdge(STRAIGHT_CURVE_TYPE, v1, v2);

  DLIList<RefEntity*> ents_to_sweep;
  ents_to_sweep.insert(line);

  DLIList<Body*> new_bodies;
  gmt->sweep_rotational(ents_to_sweep, CubitVector(0,0,0), CubitVector(0,0,1), 2*CUBIT_PI, new_bodies, CUBIT_FALSE, CUBIT_FALSE);


  //calculate the offset we want based on a&b
  double offset = a*a/b;
  //now create the plane
  Body* plane = gmt->planar_sheet(CubitVector(offset,-a,-b),CubitVector(offset,-a,b),
				  CubitVector(offset,a,b),CubitVector(offset,a,-b));


  DLIList<RefFace*> surfs;
  gqt->ref_faces(surfs);
  DLIList<RefEdge*> edge_list;

  gmt->surface_intersection(surfs[0],surfs[1],edge_list);


  gqt->delete_single_Body(new_bodies[0]);
  gqt->delete_single_Body(plane);

    //should be created by now, time to export
  DLIList<RefEntity*> exp_bodies;
  int exp_ents;
  CubitString cubit_version("12.2");
  
  CubitCompat_export_solid_model(exp_bodies, "Hyperbola.sat", "ACIS_SAT", exp_ents, cubit_version);


  //now create a curve from the resulting 

  return;
}
