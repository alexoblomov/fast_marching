cmake_minimum_required(VERSION 3.15)
project(RegularPolygons)
set(CMAKE_CXX_STANDARD 14)

find_package(CGAL)

add_executable(RegularPolygons main.cpp)
target_link_libraries(RegularPolygons CGAL::CGAL)
