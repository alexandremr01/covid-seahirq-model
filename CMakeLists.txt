cmake_minimum_required(VERSION 3.10)

# set the project name
project(covid-seahirq-model)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set( CMAKE_CXX_FLAGS "-fPIC -Wall" )

# includes EIGEN
SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )


# library
add_library(modellib STATIC src/seahirq_lib/solvers.cpp src/seahirq_lib/scenario.cpp src/seahirq_lib/models.cpp)

# executable
add_executable(spatial_covid0d_estrat src/cmd/spatial_covid0d_estrat.cpp src/cmd/real_data_reader.cpp)
target_link_libraries (spatial_covid0d_estrat modellib)

# python implementation
find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module(cmodels MODULE src/python_interface/bind.cpp src/python_interface/calculate_r0.cpp
        src/python_interface/matrix.cpp src/python_interface/run_model.cpp)
target_link_libraries(cmodels PRIVATE modellib)