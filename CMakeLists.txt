## =========================================================================
## @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
## =========================================================================

## -- Base configuration
cmake_minimum_required(VERSION 3.18)
project(SimplexMesh VERSION 1.0.0.1)

## -- CMake directives
set(CMAKE_CXX_STANDARD 17)
include(GenerateExportHeader)

## -- Find Eigen3
find_package(Eigen3 REQUIRED)

## -- Find CGAL
find_package(CGAL REQUIRED COMPONENTS Core Qt5)
if(CGAL_FOUND AND CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
else()
  message(
    FATAL_ERROR
    "ERROR: this program requires CGAL and CGAL_Qt5 and will not be compiled."
    )
endif()

add_executable(build_simplex_mesh build_simplex_mesh.cxx)
target_link_libraries(
  build_simplex_mesh Eigen3::Eigen CGAL::CGAL CGAL::CGAL_Core CGAL::CGAL_Qt5
  )

## eof - $RCSfile$
