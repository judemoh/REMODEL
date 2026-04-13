cmake_minimum_required(VERSION 3.20)
project(WP1_NURBS_Partition)

enable_language(CXX Fortran)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

# Point cmake directly at the conda OpenCASCADE installation
set(CMAKE_PREFIX_PATH "/Users/jude/anaconda3/envs/cad_hpc/lib/cmake/opencascade")
find_package(OpenCASCADE REQUIRED COMPONENTS FoundationClasses)
link_directories(/Users/jude/anaconda3/envs/cad_hpc/lib)
add_definitions(-DUSE_OpenCASCADE)


set(INCLUDE_DIRS ${OpenCASCADE_INCLUDE_DIRS})
set(LINK_LIBRARIES ${OpenCASCADE_LIBRARIES})

set(NURBS_ROOT "${CMAKE_SOURCE_DIR}/../NURBS_Fortran")

set(LIB_SRC
    ${NURBS_ROOT}/occ_interface.hpp
    ${NURBS_ROOT}/occ_interface.cpp
    ${NURBS_ROOT}/utils.f90
    ${NURBS_ROOT}/nurbsCurves.f90
    ${NURBS_ROOT}/nurbsSurfaces.f90
    ${NURBS_ROOT}/occ_interface.f90
    ${NURBS_ROOT}/occ_reader.f90
)

add_library(occ_interface ${LIB_SRC})
target_include_directories(occ_interface PRIVATE ${INCLUDE_DIRS} ${NURBS_ROOT})
target_link_libraries(occ_interface PRIVATE ${LINK_LIBRARIES})

add_executable(main
    main.f90
    nurbsPartition_JH.f90
)
target_link_libraries(main PRIVATE ${LINK_LIBRARIES} occ_interface)