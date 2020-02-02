cmake_minimum_required(VERSION 3.15)
project(RegularPolygons)
set(CMAKE_CXX_STANDARD 14)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_package(Boost REQUIRED COMPONENTS timer)
find_package(CGAL)

add_executable(RegularPolygons main.cpp utils.h utils.cpp)
target_link_libraries(RegularPolygons CGAL::CGAL)
target_link_libraries(RegularPolygons ${Boost_TIMER_LIBRARY} ${Boost_SYSTEM_LIBRARY})
