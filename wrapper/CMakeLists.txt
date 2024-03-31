add_executable(swigwrapper-test-bp libdai-swig-wrapper-test-bp.cpp)
target_link_libraries(swigwrapper-test-bp dai dai_ext)

add_executable(swigwrapper-test-em libdai-swig-wrapper-test-em.cpp)
target_link_libraries(swigwrapper-test-em dai dai_ext)

add_executable(fgcheck fgcheck.cpp)
target_link_libraries(fgcheck dai)

if(NOT DEFINED JAVA_HOME)
    message(STATUS "Finding java automatically")
    find_package(Java 11 COMPONENTS Development REQUIRED)
    message(STATUS "Java found: ${Java_JAVA_EXECUTABLE} (version: ${Java_VERSION})")
    cmake_path(GET Java_JAVA_EXECUTABLE PARENT_PATH java_bin)
    cmake_path(GET java_bin PARENT_PATH JAVA_HOME)
endif()
message(STATUS "Java home: ${JAVA_HOME}")
find_package(JNI 11 REQUIRED)
message(STATUS "JNI found: ${JNI_INCLUDE_DIRS}")

include(UseSWIG)
set_property(SOURCE libdai-swig-wrapper.i PROPERTY CPLUSPLUS ON)
set_property(SOURCE libdai-swig-wrapper.i PROPERTY COMPILE_OPTIONS -package ${JAVA_PACKAGE})
swig_add_library(daifg
        TYPE SHARED
        LANGUAGE java
        SOURCES libdai-swig-wrapper.i)
target_include_directories(daifg SYSTEM
        PRIVATE ${JNI_INCLUDE_DIRS}
        PRIVATE .)
target_link_libraries(daifg PRIVATE dai dai_ext)

install(TARGETS daifg DESTINATION jnilib)
