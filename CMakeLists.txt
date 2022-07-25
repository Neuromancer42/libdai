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
add_library(dai_static STATIC ${SOURCES})
target_link_libraries(dai ${LIBS})

add_executable(swigwrapper wrapper/libdai-swig-wrapper-test.cpp)
target_link_libraries(swigwrapper dai)

add_library(daifg SHARED wrapper/libdai-swig-wrapper_wrap.cxx)
target_include_directories(daifg  SYSTEM
        PRIVATE "/Library/Java/JavaVirtualMachines/jdk-11.0.14.jdk/Contents/Home/include"
        PRIVATE "/Library/Java/JavaVirtualMachines/jdk-11.0.14.jdk/Contents/Home/include/darwin"
        )
target_link_libraries(daifg dai_static ${LIBS})
