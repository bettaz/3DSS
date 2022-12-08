cmake_minimum_required(VERSION 2.9)

project( Second_Assignment )
set(CMAKE_SOURCE_DIR "./src")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

set(OpenCV_DIR "~/openCVbuild")
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
find_package(OpenMP)

add_executable ( Second_Assignment src/main.cpp)
target_link_libraries ( Second_Assignment ${OpenCV_LIBS} )
