//
// Created by alex on 02/02/2020.
//

#pragma once
#include <vector>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 point;
typedef std::vector<float> point3D;

void generate_points_in_triangle(float Ax, float Bx, float Cx, float Ay, float By, float Cy,
                                 unsigned int num_points, std::vector<point> &points);
void generate_points_on_segment(float x, float x_prime, float y, float y_prime, unsigned int num_points,
                                std::vector<point> &points, std::vector<point3D> &boundary_points);
inline float DegToRad(float theta);

void rotate_xy(float ix, float iy, float &ox, float &oy, float angle);

void trace_regular_polygon(int degree, float x_min, float x_max, float y_min, float y_max,
                           unsigned int num_boundary_points, unsigned int num_interior_points,
                           std::vector<point> &points, std::vector<point3D> &boundary_points);

void export_connectivity_to_txt(std::ofstream &stream, std::vector<point3D> &boundary_points, int degree,
                                unsigned int num_boundary_points);
void export_boundary_points(std::ofstream &stream, std::vector<point3D> &points);
