#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <utility> // defines std::pair
#include <vector>
#include <fstream>
#include <iostream>
#include <CGAL/grid_simplify_point_set.h>

#include <CGAL/Point_set_3.h>
using TScalar = double;
using TKernel = CGAL::Simple_cartesian< TScalar >;
using TPoint  = TKernel::Point_3;
using TPointSet = CGAL::Point_set_3< TPoint >;
using TVector = TKernel::Vector_3;

#include <CGAL/draw_point_set_3.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/compute_average_spacing.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/disable_warnings.h>
typedef std::array<std::size_t,3> Facet;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3  Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
// Point with normal vector stored as a std::pair.
typedef std::pair<Point, Vector> Pwn;

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
namespace SMS = CGAL::Surface_mesh_simplification;

typedef CGAL::Surface_mesh<Point_3>                  Surface_mesh;


struct Construct{
  Mesh& mesh;
  template < typename PointIterator>
  Construct(Mesh& mesh,PointIterator b, PointIterator e)
    : mesh(mesh)
  {
    for(; b!=e; ++b){
      boost::graph_traits<Mesh>::vertex_descriptor v;
      v = add_vertex(mesh);
      mesh.point(v) = *b;
    }
  }
  Construct& operator=(const Facet f)
  {
    typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits<Mesh>::vertices_size_type size_type;
    mesh.add_face(vertex_descriptor(static_cast<size_type>(f[0])),
                  vertex_descriptor(static_cast<size_type>(f[1])),
                  vertex_descriptor(static_cast<size_type>(f[2])));
    return *this;
  }
  Construct&
  operator*() { return *this; }
  Construct&
  operator++() { return *this; }
  Construct
  operator++(int) { return *this; }
};

int main(int argc, char*argv[])
{
  const std::string fname = (argc>1) ? argv[1] : "../depth.ply";
  // Reads a .xyz point set file in points[].
  // Note: read_points() requires an output iterator
  // over points and as well as property maps to access each
  // point position and normal.
  std::vector<Pwn> points;
  std::vector<Pwn> pointsFilter;
  Surface_mesh surface_mesh;
  if(!CGAL::IO::read_PLY(fname,
                         std::back_inserter(points),
                         CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                          .normal_map(CGAL::Second_of_pair_property_map<Pwn>())))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  TPointSet filter_points;
  std::vector<Point_3> pointsTomesh;
  std::vector< std::pair<Point,unsigned> > pointsDeluanay;
  Mesh m;
  for(std::size_t i = 0; i < points.size(); ++ i){
    if(points[i].first[2] != 0){
      TPoint values (points[i].first[0], points[i].first[1], points[i].first[2]);
      Point_3 valuesTwo (points[i].first[0], points[i].first[1], points[i].first[2]);
      filter_points.insert(values);
      pointsTomesh.push_back(valuesTwo);
    }
  }

  /*

  // simplification by clustering using erase-remove idiom
  double cell_size = 0.03;
  unsigned int min_points_per_cell = 3;
  auto iterator_to_first_to_remove
    = CGAL::grid_simplify_point_set
    (pointsTomesh, cell_size); // optional
  pointsTomesh.erase(iterator_to_first_to_remove, pointsTomesh.end());
  // Optional: after erase(), shrink_to_fit to trim excess capacity
  pointsTomesh.shrink_to_fit();

  */

 Construct construct(m,pointsTomesh.begin(),pointsTomesh.end());

 CGAL::advancing_front_surface_reconstruction(pointsTomesh.begin(),
                                               pointsTomesh.end(),
                                               construct);

  CGAL::draw(m);
  return EXIT_SUCCESS;
}