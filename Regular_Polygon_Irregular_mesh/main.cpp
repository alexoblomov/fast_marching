#include <iostream>
#include <vector>
#include <random>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Surface_mesh.h>


#include "Triangulation_off_ostream_2.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3 point;


void generate_points_in_domain(float x_min, float x_max, float y_min, float y_max,
                               unsigned int num_points, std::vector<point> &points)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> unif_dist_x(x_min, x_max);
    std::uniform_real_distribution<float> unif_dist_y(y_min, y_max);

    float x,y;
    for(unsigned int i = 0; i< num_points; ++i)
    {
        x = unif_dist_x(rng);
        y = unif_dist_y(rng);
        points.push_back(point(x,y,0));
    }

}
/*
 * TO replace by generate points on segment
 */
void generate_points_on_boundary(float x_min, float x_max, float y_min, float y_max,
                               unsigned int num_points,
                               std::vector<point> &points)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> unif_dist_x(x_min, x_max);
    std::uniform_real_distribution<float> unif_dist_y(y_min, y_max);
    std::uniform_real_distribution<float> test(y_max, y_max);
    float x,y;
    float test_y = test(rng);

    points.push_back(point(x_min,y_min,0));
    points.push_back(point(x_max,y_min,0));
    points.push_back(point(x_min,y_max,0));
    points.push_back(point(x_max,y_max,0));
    for(unsigned int i = 0; i< num_points/4; ++i)
    {
        x = unif_dist_x(rng);
        points.push_back(point(x,y_min,0));
    }
    for(unsigned int i = 0; i< num_points/4; ++i)
    {
        x = unif_dist_x(rng);
        points.push_back(point(x,y_max,0));
    }
    for(unsigned int i = 0; i< num_points/4; ++i)
    {
        y = unif_dist_y(rng);
        points.push_back(point(x_min,y,0));
    }
    for(unsigned int i = 0; i< num_points/4; ++i)
    {
        y = unif_dist_y(rng);
        points.push_back(point(x_max,y,0));
    }

}

void generate_points_on_segment(float x_a, float x_b, float y_a, float y_b,
                                 unsigned int num_points,
                                std::vector<point> &points)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> unif_dist_x(x_a,x_b);
    std::uniform_real_distribution<float> unif_dist_y(y_a,y_b);

    float x,y;
    for(unsigned int i = 0; i< num_points; ++i)
    {
        x = unif_dist_x(rng);
        y = unif_dist_y(rng);
        points.push_back(point(x,y,0));
    }
}
//
//void trace_ray(point &ipoint_a, float itheta, float &opointb)
//{
//
//}

void print_points(std::vector<point> &points)
{
    for(auto v: points)
    {
        std::cout << v[0] << " " << v[1] << std::endl;
    }
}

int main() {
    int n = 0;
    std::cout << "Enter degree of polygon (must be even): ";
    std::cin >> n;
    if(n%2 == 1)
    {
        std::cout << "not what we want " << std::endl;
    }
    // boundary of fundamental domain:
    float x_min = 0., x_max = 1., y_min = 0., y_max = 1.;
//    std::cout << "Enter boundary of fundamental domain in the following order: x_min, x_max, y_min, y_max ";
//    std::cin >> x_min >> x_max >> y_min >> y_max;
    std::vector<point> points;

    unsigned int num_interior_points = 40;
    generate_points_in_domain(x_min, x_max, y_min, y_max, num_interior_points,points);

    unsigned int num_boundary_points = 30;
    // will be used to construct rays.
    //std::vector<point> boundary_points;
    generate_points_on_boundary(x_min, x_max, y_min, y_max, num_boundary_points, points);
    std::cout << "with boundary points: " << std::endl;
    Delaunay triangulation(points.begin(),points.end());

    std::cout << triangulation.number_of_vertices() << std::endl;
    std::cout << triangulation.number_of_faces() << std::endl;

    std::ofstream off_stream("tst.off");
    export_triangulation_2_to_off(off_stream, triangulation);

    return 0;

//    std::cout << "Constructing polygon of degree " << n << std::endl;
//
//    float x_center = (x_max - x_min)/2;
//    float y_center = (y_max - y_min)/2;
//    // the ray between the bounding box and the edge of the regular 2n-gon
//    float theta = 180. - (180.-360./n); // see notes for more detailed explanation

    return 0;
}

