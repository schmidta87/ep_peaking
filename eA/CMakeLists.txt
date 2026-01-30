# Add executable for converting generator output to lund format.
add_executable(generator_eHe4_born generator_eHe4_born.cpp )
target_link_libraries(generator_eHe4_born ${ROOT_LIBRARIES} )
