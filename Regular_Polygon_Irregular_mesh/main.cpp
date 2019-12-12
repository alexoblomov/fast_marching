#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <iterator>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Surface_mesh.h>
//#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <boost/graph/adjacency_list.hpp>
#include "Triangulation_off_ostream_2.h"

//typedef std::vector<float> point;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef K::Point_3 point;
const float Pi = 3.141592654f;


void generate_points_in_triangle(float Ax, float Bx, float Cx, float Ay, float By, float Cy,
                              unsigned int num_points, std::vector<point> &points)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    float tol = 0.01; // ensures we don't create a point on the boundary
    std::uniform_real_distribution<float> uniformRealDistribution_x(0 + tol, 1);
    std::uniform_real_distribution<float> uniformRealDistribution_y(0 + tol, 1);

    float x, y , r1, r2;
    for (unsigned int i = 0; i < num_points; ++i)
    {
        r1 = uniformRealDistribution_x(rng);
        r2 = uniformRealDistribution_y(rng);
        x = (1- sqrt(r1))*Ax + (sqrt(r1)*(1-r2))*Bx + (r2*sqrt(r1))*Cx;
        y = (1- sqrt(r1))*Ay + (sqrt(r1)*(1-r2))*By + (r2*sqrt(r1))*Cy;
        points.emplace_back(x,y,0);
    }

}

void generate_points_on_segment(float x, float x_prime, float y, float y_prime, unsigned int num_points,
                                std::vector<point> &points)
{
    float px, py;
    float dx = (x_prime - x) / num_points;
    float dy = (y_prime - y) / num_points;

    for (unsigned int i = 0; i <= num_points; ++i)
    {
        px = x + i * dx;
        py = y + i * dy;
        points.emplace_back(px,py,0);
    }
}

inline float DegToRad(float theta)
{
    return theta / 180 * Pi;
}

void rotate_xy(float ix, float iy, float &ox, float &oy, float angle)
{
    float theta = DegToRad(angle);
    ox = ix * std::cos(theta) - iy * std::sin(theta);
    oy = ix * std::sin(theta) + iy * std::cos(theta);
}

void trace_regular_polygon(int degree, float x_min, float x_max, float y_min, float y_max,
                           unsigned int num_boundary_points, unsigned int num_interior_points,
                           std::vector<point> &points)
{
    float x_center = (x_max - x_min) / 2;
    float y_center = (y_max - y_min) / 2;
    float interior_angle = 360./ degree;
    float x = x_max/2;
    float y = y_max;
    float x_prime = 0;
    float y_prime = 0;
    for (int i = 0; i <= degree; i++)
    {
        rotate_xy(x,y,x_prime,y_prime,interior_angle);
        generate_points_on_segment(x, x_prime, y, y_prime, num_boundary_points, points);
        generate_points_in_triangle(x_center, x, x_prime,y_center, y, y_prime,num_interior_points,
                points);
        x = x_prime;
        y = y_prime;
    }

}

int main()
{
    int n = 0;
    std::cout << "Enter degree of polygon (must be even and >= 4): ";
    std::cin >> n;
    if (n % 2 == 1)
    {
        std::cout << "not what we want " << std::endl;
    }
    // boundary of fundamental domain:
    float x_min = 0., x_max = 1., y_min = 0., y_max = 1.;


    unsigned int num_interior_points = 40;
    unsigned int num_boundary_points = 8;
//    add the following as outputs in the future for djikstra
//    std::vector<point> boundary_points;
//    std::vector<point> interior_points;

    // test polygon
    std::vector<point> points;

    trace_regular_polygon(n, x_min, x_max, y_min, y_max, num_boundary_points, num_interior_points, points);
    Delaunay triangulation(points.begin(), points.end());

    std::cout << triangulation.number_of_vertices() << std::endl;
    std::cout << triangulation.number_of_faces() << std::endl;

    std::ofstream off_stream("tst_poly.off");
    export_triangulation_2_to_off(off_stream, triangulation);


    // todo: change connectivity for dijkstra
    return 0;
}
