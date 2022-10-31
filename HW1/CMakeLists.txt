cmake_minimum_required(VERSION 2.9)

project( First_Assignment )
set(CMAKE_SOURCE_DIR "./build")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

#set(OpenCV_DIR "~/openCVbuild")
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
find_package(OpenMP)

add_executable ( First_Assignment src/main.cpp src/StereoMatcher.cpp src/StereoMatcher.h src/PointCloud.cpp src/PointCloud.h)
target_link_libraries ( First_Assignment ${OpenCV_LIBS} )
