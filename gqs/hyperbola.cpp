
#include "gq.hpp"



void hyperbolic_curves(double a, double b, DLIList<RefEdge*> &edge_list);

int main()
{

  DLIList<RefEdge*> curves;
  hyperbolic_curves(3,5, curves);

  //should be created by now, time to export
  DLIList<RefEntity*> exp_bodies;
  int exp_ents;
  CubitString cubit_version("12.2");
  
  CubitCompat_export_solid_model(exp_bodies, "Hyperbola.sat", "ACIS_SAT", exp_ents, cubit_version);


  return 0;
}

/* Creates two hyperbolic curves in the xy plane using the parameters a & b

Hyperbolic Curve Form:

x^2/a - y^2/b = 1

symmetric axis - x
reflection axis - y

Returns: two RefEdge pointers to the curves

*/

void hyperbolic_curves(double a, double b, DLIList<RefEdge*> &edge_list)
{

  //first create a conic surface
  CubitVector p1(b,0,0);
  CubitVector p2(0,a,0);

  RefVertex* v1 = gmt->make_RefVertex(p1);

  RefVertex* v2 = gmt->make_RefVertex(p2);

  RefEdge* line = gmt->make_RefEdge(STRAIGHT_CURVE_TYPE, v1, v2);

  DLIList<RefEntity*> ents_to_sweep;
  ents_to_sweep.insert(line);

  DLIList<Body*> new_bodies;
  gmt->sweep_rotational(ents_to_sweep, CubitVector(0,0,0), CubitVector(1,0,0), 2*CUBIT_PI, new_bodies, CUBIT_FALSE, CUBIT_FALSE);


  //calculate the offset we want based on a & b to give us a pure hyperbole
  double offset = a*a/b;

  //now create the plane
  Body* plane = gmt->planar_sheet(CubitVector(-b,-a,offset),CubitVector(b,-a,offset),
				  CubitVector(b,a,offset),CubitVector(-b,a,offset));


  DLIList<RefFace*> surfs;
  gqt->ref_faces(surfs);

  gmt->surface_intersection(surfs[0],surfs[1],edge_list);


  gqt->delete_single_Body(new_bodies[0]);
  gqt->delete_single_Body(plane);

  //make a copy of the first curve
  RefEdge *copy = gmt->make_RefEdge( edge_list[0], true);

  //now reflect this curve across the creation plane
  CubitVector reflection_pt(b,0,0);
  CubitVector reflection_ax(1,0,0);

  RefEntity* copy_ent = dynamic_cast<RefEntity*>(copy);
  DLIList<RefEntity*> ents_to_reflect, reflected_ents;
  ents_to_reflect.insert(copy_ent);

  gqt->reflect(ents_to_reflect, reflection_pt, reflection_ax, true, reflected_ents);
  
  //add the reflected curve to the output list
  edge_list.insert(copy);
  
  //now create a curve from the resulting 

  return;
}

/* this function will return hyperbolic curves in a plane of two principle axes



axis 1 (ax1) - axis of symmetry for one of the curves

axis 2 (ax1) - reflecting axis for the curves

(for the axes arguments: 0 is x, 1 is y, 2 is z

*/
void hyperbolic_curves_in_plane( double a, double b, int ax1, int ax2, DLIList<RefEdge*> &edge_list)
{

  //first create our curves in the yz plane
  //symmetric axis - x
  //reflection axis - y
  hyperbolic_curves(a,b, edge_list);


}
