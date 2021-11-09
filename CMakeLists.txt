cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(sacryst LANGUAGES C CXX)
include_directories(src)

set(sacryst_SRCS 
  src/main.cpp
  src/mt19937ar_c.cpp
)

add_executable(sacryst ${sacryst_SRCS})
set_property(TARGET sacryst PROPERTY CXX_STANDARD 17)
set_property(TARGET sacryst PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET sacryst PROPERTY CUDA_SEPARABLE_COMPILATION ON)