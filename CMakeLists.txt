cmake_minimum_required(VERSION 3.6)
project(rb130102d)
#-std=c++11
#add_definitions(-fopenmp -lm)

# Add openmp flags
FIND_PACKAGE( OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

# Setup dirs
set(CMAKE_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Compile every .c file to single target
file( GLOB APP_SOURCES *.c )
foreach( testsourcefile ${APP_SOURCES} )
    # I used a simple string replace, to cut off .cpp.
    string( REPLACE ".c" "" testname ${testsourcefile} )
    get_filename_component( testname_strip ${testname} NAME )
    get_filename_component( testsourcefile_strip ${testsourcefile} NAME )
    add_executable( ${BIN_PREFIX}${testname_strip} ${testsourcefile_strip} )
    # Make sure YourLib is linked to each app
    # target_link_libraries( ${testname} YourLib )
endforeach( testsourcefile ${APP_SOURCES} )