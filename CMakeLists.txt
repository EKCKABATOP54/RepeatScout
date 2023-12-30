cmake_minimum_required(VERSION 3.27)
project(bioinformatics)

set(CMAKE_C_COMPILER clang)
set(CMAKE_CPP_COMPILER clang++)

set(CMAKE_CXX_STANDARD 20)

add_executable(bioinformatics main.cpp
        repeat_scout.cpp
        repeat_scout.h)
