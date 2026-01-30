# CMakeLists for the radiative ep generator under the peaking approximation

# Set up
cmake_minimum_required (VERSION 3.0)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

# Project name
project (generators)

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for running the radiative generator
add_executable(generator generator.cpp)
target_link_libraries(generator ${ROOT_LIBRARIES} )

# Add executable for running the born generator
add_executable(generator_born generator_born.cpp)
target_link_libraries(generator_born ${ROOT_LIBRARIES} )

# Add subdirectory for lund file conversion
add_subdirectory(lund_conversion)

# Add subdirectory for elastic generators on nuclei
add_subdirectory(eA)
