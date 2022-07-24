cmake_minimum_required(VERSION 3.15)
project(libdai)

set(CMAKE_CXX_STANDARD 11)

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS program_options)
message(STATUS  "Boost version: ${Boost_VERSION}")

include_directories(${Boost_INCLUDE_DIR})

find_library(GMP gmp REQUIRED)
message(STATUS "GMP found: ${GMP}")
find_library(GMPXX gmpxx REQUIRED)
message(STATUS "GMPXX found: ${GMPXX}")
list(APPEND LIBS ${GMPXX} ${GMP})

include_directories(include)

file(GLOB SOURCES src/*.cpp)
add_library(dai SHARED ${SOURCES})
target_link_libraries(dai ${LIBS})

add_executable(multi-predictor wrapper/multipredict-wrapper.cpp)
target_link_libraries(multi-predictor dai)

add_executable(predictor wrapper/predict-wrapper.cpp)
target_link_libraries(predictor dai)

add_executable(updater wrapper/update-wrapper.cpp)
target_link_libraries(updater dai)