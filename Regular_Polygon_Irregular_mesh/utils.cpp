//
// Created by alex on 02/02/2020.
//
#include "utils.h"

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <iterator>

#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K> Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;

typedef K::Point_3 point;
typedef std::vector<float> point3D;
const float Pi = 3.141592654f;


void generate_points_in_triangle(float Ax, float Bx, float Cx, float Ay, float By, float Cy,
                                 unsigned int num_points, std::vector<point> &points)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    float tol = 0.0001; // ensures we don't duplicate points in the center
    std::uniform_real_distribution<float> uniformRealDistribution_x(0 + tol, 1);
    std::uniform_real_distribution<float> uniformRealDistribution_y(0 + tol, 1);

    float x, y, r1, r2;
    for (unsigned int i = 0; i < num_points; ++i)
    {
        r1 = uniformRealDistribution_x(rng);
        r2 = uniformRealDistribution_y(rng);
        x = (1 - sqrt(r1)) * Ax + (sqrt(r1) * (1 - r2)) * Bx + (r2 * sqrt(r1)) * Cx;
        y = (1 - sqrt(r1)) * Ay + (sqrt(r1) * (1 - r2)) * By + (r2 * sqrt(r1)) * Cy;
        points.emplace_back(x, y, 0);
    }

}

void generate_points_on_segment(float x, float x_prime, float y, float y_prime, unsigned int num_points,
                                std::vector<point> &points, std::vector<point3D> &boundary_points)
{
    float px, py;
    float dx = (x_prime - x) / num_points;
    float dy = (y_prime - y) / num_points;

    for (unsigned int i = 0; i <= num_points; ++i)
    {
        px = x + i * dx;
        py = y + i * dy;
        points.emplace_back(px, py, 0);
        point3D p{px, py, 0};
        boundary_points.push_back(p);
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
                           std::vector<point> &points, std::vector<point3D> &boundary_points)
{
    float x_center = (x_max - x_min) / 2;
    float y_center = (y_max - y_min) / 2;
    float interior_angle = 360. / degree;
    float x = x_max / 2;
    float y = y_max;
    float x_prime = x_min;
    float y_prime = y_max;
    for (int i = 0; i < degree; i++) // should be <
    {
        rotate_xy(x, y, x_prime, y_prime, interior_angle);
        generate_points_on_segment(x, x_prime, y, y_prime, num_boundary_points, points, boundary_points);
        generate_points_in_triangle(x_center, x, x_prime, y_center, y, y_prime, num_interior_points, points);
        x = x_prime;
        y = y_prime;
    }
}


void export_connectivity_to_txt(std::ofstream &stream, std::vector<point3D> &boundary_points, int degree,
                                unsigned int num_boundary_points)
{
    std::vector<point3D> connect;
    int count = 0;
    for (unsigned int d = 0; d < degree / 2; d++)
    {
        for (unsigned int p_idx = 0; p_idx < num_boundary_points; p_idx++)
        {
            connect.push_back(boundary_points[p_idx + d * (num_boundary_points + 1)]);
            connect.push_back(boundary_points[(d + 1 + degree / 2) * (num_boundary_points + 1) - p_idx - 1]);
        }
        count++;

    }
    for (const auto &vt : connect)
    {
        std::copy(vt.cbegin(), vt.cend(),
                  std::ostream_iterator<float>(stream, " "));
        stream << '\n';
    }
}

void export_boundary_points(std::ofstream &stream, std::vector<point3D> &points)
{
    for (const auto &vt : points)
    {
        std::copy(vt.cbegin(), vt.cend(),
                  std::ostream_iterator<float>(stream, " "));
        stream << '\n';
    }
}
