#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <boost/timer/timer.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include "Triangulation_off_ostream_2.h"
#include "utils.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef K::Point_3 point;
typedef std::vector<float> point3D;

int main()
{
    // boundary of fundamental domain:
    float x_min = 0., x_max = 1., y_min = 0., y_max = 1.;

    unsigned int num_interior_points = 10000;
    unsigned int num_boundary_points = 500;

    boost::timer::auto_cpu_timer t;
    for (int n = 4; n < 12; n += 2)
    {
        std::vector<point3D> boundary_points;
        std::vector<point> points;

        trace_regular_polygon(n, x_min, x_max, y_min, y_max, num_boundary_points, num_interior_points,
                              points, boundary_points);
        Delaunay triangulation(points.begin(), points.end());

        std::string f_name = "poly_" + std::to_string(n) + ".off";
        std::ofstream off_stream(f_name);
        export_triangulation_2_to_off(off_stream, triangulation);

        f_name = "poly_connectivity_" + std::to_string(n) + ".txt";
        std::ofstream c_stream(f_name);
        export_connectivity_to_txt(c_stream, boundary_points, n, num_boundary_points);

        f_name = "poly_boundary_" + std::to_string(n) + ".txt";
        std::ofstream bd_stream(f_name);
        export_boundary_points(bd_stream, boundary_points);
    }

    return 0;
}

