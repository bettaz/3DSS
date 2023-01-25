cmake_minimum_required(VERSION 2.9)

project( Third_Assignment )
set(CMAKE_SOURCE_DIR "./src")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

find_package(OpenMP)

add_executable ( Third_Assignment src/main.cpp src/pcManager.h)
target_link_libraries ( Third_Assignment ${OpenCV_LIBS} )
cmake_path(GET CMAKE_SOURCE_DIR PARENT_PATH root)
file(COPY ${root}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})