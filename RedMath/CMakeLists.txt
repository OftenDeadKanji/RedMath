cmake_minimum_required(VERSION 3.23)
project(RedMathLib VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory("RedMath/" "${CMAKE_SOURCE_DIR}/out/RedMath/")
add_subdirectory("Testing/" "${CMAKE_SOURCE_DIR}/out/Testing/")