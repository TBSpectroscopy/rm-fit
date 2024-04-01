cmake_minimum_required(VERSION 3.5)
project(calc_g)

#set(pybind11_DIR /usr/include/pybind11)
find_package(pybind11 REQUIRED)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR})

add_compile_options(
    -O2
    -Wall
    -Wformat=0
    -flto=auto
)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
pybind11_add_module(calc_g calc_g.cpp)
target_link_libraries(calc_g PRIVATE Eigen3::Eigen)
