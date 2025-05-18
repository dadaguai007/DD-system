cmake_minimum_required(VERSION 3.8)

project(ML)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 99)

set(RANDOM_Src
    "ML.cpp"
)

set(RANDOM_Inc
    "ML.hpp"
)
add_definitions("-O3")

add_library(ML_lib STATIC 
    ${RANDOM_Src} 
    ${RANDOM_Inc} 
)


