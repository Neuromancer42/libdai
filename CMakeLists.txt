cmake_minimum_required(VERSION 3.20)
project(libdai)

set(CMAKE_CXX_STANDARD 11)

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS program_options)
message(STATUS  "Boost found: ${Boost_INCLUDE_DIR} (version: ${Boost_VERSION})")

include_directories(${Boost_INCLUDE_DIR})

find_library(GMP gmp REQUIRED)
message(STATUS "GMP found: ${GMP}")
find_library(GMPXX gmpxx REQUIRED)
message(STATUS "GMPXX found: ${GMPXX}")
list(APPEND LIBS ${GMPXX} ${GMP})
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    list(APPEND LIBS OpenMP::OpenMP_CXX)
endif()
include_directories(include)

file(GLOB SOURCES src/*.cpp)
add_library(dai ${SOURCES})
set_property(TARGET dai PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(dai ${LIBS})
install(TARGETS dai DESTINATION lib)

if(NOT WRAPPER_ONLY)
    add_subdirectory(utils)
endif()

add_subdirectory(ext)

find_package(SWIG 4.0 COMPONENTS java)
if(SWIG_FOUND)
    message(STATUS "SWIG found: ${SWIG_EXECUTABLE} (version: ${SWIG_VERSION})")
    add_subdirectory(wrapper)
endif()
