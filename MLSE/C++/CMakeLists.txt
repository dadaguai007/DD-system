cmake_minimum_required(VERSION 3.10)
project(ML)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_definitions("-O3 -ffast-math")
set(CMAKE MODULE_PATH "/usr/share/cmake/Modules/")
 

include_directories(
    src
)
add_subdirectory("src")
add_executable(test_ML
    "test_ML.cpp"
)
target_link_libraries(test_ML PUBLIC ML_lib)
