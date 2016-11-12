cmake_minimum_required(VERSION 3.6)
project(rb130102d)
add_definitions(-std=c++11 -fopenmp -lm)

set(CMAKE_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file( GLOB APP_SOURCES *.cpp )
foreach( testsourcefile ${APP_SOURCES} )
    # I used a simple string replace, to cut off .cpp.
    string( REPLACE ".cpp" "" testname ${testsourcefile} )
    get_filename_component( testname_strip ${testname} NAME )
    get_filename_component( testsourcefile_strip ${testsourcefile} NAME )
    add_executable( ${BIN_PREFIX}${testname_strip} ${testsourcefile_strip} )
    # Make sure YourLib is linked to each app
    # target_link_libraries( ${testname} YourLib )
endforeach( testsourcefile ${APP_SOURCES} )