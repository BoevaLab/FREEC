cmake_minimum_required(VERSION 3.21)
project(ploidy-purity-cpp)

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

endif(NOT CMAKE_BUILD_TYPE)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -pedantic")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
option(test "Build tests." ON)

find_package(NLopt)
add_definitions(-DUSE_NLOPT)
add_definitions(-DPROFILE_TRACE)
add_definitions(-DCLUSTER)

file(GLOB ppc_src
        "src/*.h"
        "src/*.cpp"
        )

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_executable(PPC src/main.cpp ${ppc_src})

target_include_directories(PPC PUBLIC limbo/src)
target_link_libraries(PPC Eigen3::Eigen)

get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Boost 1.74.0 COMPONENTS date_time filesystem system)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(PPC Eigen3::Eigen ${Boost_LIBRARIES})
target_include_directories(PPC PUBLIC ${NLOPT_INCLUDE_DIRS})
target_link_libraries(PPC ${NLOPT_LIBRARIES})
