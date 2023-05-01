cmake_minimum_required(VERSION 3.4)
project(calc_g)

set(pybind11_DIR /usr/include/pybind11)
find_package(pybind11 REQUIRED)

add_compile_options(
    -O2
    -Wall
    -Wformat=0
    -flto=auto
)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
pybind11_add_module(calc_g calc_g.cpp)
