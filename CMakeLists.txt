cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(sacryst LANGUAGES C CXX)
include_directories(src)

set(sacryst_SRCS 
  src/main.cpp
  src/cParameters.cpp
)

add_executable(sacryst ${sacryst_SRCS})
target_compile_definitions(sacryst PRIVATE PRINT_INTERMEDIATE_RESULTS)
target_compile_definitions(sacryst PRIVATE WRITE_STATISTICS)
target_compile_definitions(sacryst PRIVATE CUSTOM_ENERGY)

add_executable(sacryst_ackley ${sacryst_SRCS} src/energies/ackley_energy.cpp)
target_compile_definitions(sacryst_ackley PRIVATE PRINT_INTERMEDIATE_RESULTS)
target_compile_definitions(sacryst_ackley PRIVATE WRITE_STATISTICS)

add_executable(sacryst_griewangk ${sacryst_SRCS} src/energies/griewangk_energy.cpp)
target_compile_definitions(sacryst_griewangk PRIVATE PRINT_INTERMEDIATE_RESULTS)
target_compile_definitions(sacryst_griewangk PRIVATE WRITE_STATISTICS)

add_executable(sacryst_rastrigin ${sacryst_SRCS} src/energies/rastrigin_energy.cpp)
target_compile_definitions(sacryst_rastrigin PRIVATE PRINT_INTERMEDIATE_RESULTS)
target_compile_definitions(sacryst_rastrigin PRIVATE WRITE_STATISTICS)

add_executable(sacryst_rosembrock ${sacryst_SRCS} src/energies/rosenbrock_energy.cpp)
target_compile_definitions(sacryst_rosembrock PRIVATE PRINT_INTERMEDIATE_RESULTS)
target_compile_definitions(sacryst_rosembrock PRIVATE WRITE_STATISTICS)

add_executable(sacryst_sphere ${sacryst_SRCS} src/energies/sphere_energy.cpp)
target_compile_definitions(sacryst_sphere PRIVATE PRINT_INTERMEDIATE_RESULTS)
target_compile_definitions(sacryst_sphere PRIVATE WRITE_STATISTICS)

add_executable(sacryst_weierstrass ${sacryst_SRCS} src/energies/weierstrass_energy.cpp)
target_compile_definitions(sacryst_weierstrass PRIVATE PRINT_INTERMEDIATE_RESULTS)
target_compile_definitions(sacryst_weierstrass PRIVATE WRITE_STATISTICS)

add_executable(sacryst_zakharov ${sacryst_SRCS} src/energies/zakharov_energy.cpp)
target_compile_definitions(sacryst_zakharov PRIVATE PRINT_INTERMEDIATE_RESULTS)
target_compile_definitions(sacryst_zakharov PRIVATE WRITE_STATISTICS)
