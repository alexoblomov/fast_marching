cmake_minimum_required(VERSION 3.15)
project(delaunary2D)
set(CMAKE_CXX_STANDARD 14)

find_package(CGAL)

add_executable(delaunary2D main.cpp)
target_link_libraries(delaunary2D CGAL::CGAL)
