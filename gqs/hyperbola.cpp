
#include "gq.hpp"



void hyperbolic_curves(double a, double b, DLIList<RefEdge*> &edge_list);
void hyperbolic_curves_in_plane( double a, double b, int ax1, int ax2, DLIList<RefEdge*> &edge_list);

int main()
{

  DLIList<RefEdge*> curves;
  hyperbolic_curves_in_plane(3,5,2,0,curves);


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

axis 2 (ax2) - reflecting axis for the curves

(for the axes arguments: 0 is x, 1 is y, 2 is z

*/
void hyperbolic_curves_in_plane( double a, double b, int ax1, int ax2, DLIList<RefEdge*> &edge_list)
{

  //first create our curves in the yz plane
  //symmetric axis - x
  //reflection axis - y
  hyperbolic_curves(a,b, edge_list);

  //get the curves as RefEntities
  DLIList<RefEntity*> edge_ent_list;
  if ( 2 == edge_list.size() )
    {
      edge_ent_list.insert(dynamic_cast<RefEntity*>(edge_list[0]));
      edge_ent_list.insert(dynamic_cast<RefEntity*>(edge_list[1]));
    }
  else
      std::cout << "Expected two curves. Leaving function without action..." << std::endl;

  //make sure our axes aren't the same
  if ( ax1 == ax2 )
    {
      std::cout << "Axes do not define a proper plane. Leaving function without action..." << std::endl;
      return;
    }


  //for tracking the symmetric axis during the first rotation
  int sym_axis = 0;

  //rotate the reflection axis to match the one requested
  if ( ax2 != 1 ) //if this is already the y axis, skip
    {
      //setup the rotation axis
      double ax[3] = {0,0,0};
      ax[( 2 == ax2 ) ? 0 : 2] = 1;
      double degrees = 90;

      DLIList<RefEntity*> ents_rotated;
      gqt->rotate(edge_ent_list, CubitVector(0,0,0), CubitVector(ax[0],ax[1],ax[2]), degrees, true, ents_rotated, false);

      //the symmetric axis may have changed base on this rotation
      if( ax[2] != 0 ) //if we rotated around the z-axis, the symm axis changed
	sym_axis++;
    }

  //if the symmetric axis is incorrect, rotate around reflecting axis
  if ( sym_axis != ax1 )
    {

      //setup the rotation axis
      double ax[3] = {0,0,0};
      ax[ax2] = 1;
      double degrees = 90;

      DLIList<RefEntity*> ents_rotated;
      gqt->rotate(edge_ent_list, CubitVector(0,0,0), CubitVector(ax[0],ax[1],ax[2]), degrees, true, ents_rotated, false);

    }

}
