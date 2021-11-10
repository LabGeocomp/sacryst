cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(sacryst LANGUAGES C CXX)
include_directories(src)

set(sacryst_SRCS 
  src/main.cpp
)

add_executable(sacryst ${sacryst_SRCS})
target_compile_definitions(sacryst PRIVATE PRINT_INTERMEDIATE_RESULTS)
target_compile_definitions(sacryst PRIVATE WRITE_STATISTICS)