cmake_minimum_required(VERSION 3.17)
project(cplex_test)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

######## Cplex Dependency
find_package(CPLEX REQUIRED)
if (CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIR})
    link_directories(${CPLEX_LIBRARY})
endif ()

if (CPLEX_CONCERT_FOUND)
    include_directories(${CPLEX_CONCERT_INCLUDE_DIR})
    link_directories(${CPLEX_CONCERT_LIBRARY})
endif ()

if (CPLEX_ILOCPLEX_FOUND)
    include_directories(${CPLEX_ILOCPLEX_INCLUDE_DIR})
    link_directories(${CPLEX_ILOCPLEX_LIBRARY})
endif ()

# http lib
include(FetchContent)
FetchContent_Declare(cpr GIT_REPOSITORY https://github.com/libcpr/cpr.git
        GIT_TAG 871ed52d350214a034f6ef8a3b8f51c5ce1bd400) # The commit hash for 1.9.0. Replace with the latest from: https://github.com/libcpr/cpr/releases
FetchContent_MakeAvailable(cpr)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(json)

add_subdirectory(boost-cmake)

######## Complier message
SET(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
message(STATUS "System: ${CMAKE_SYSTEM}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

######## Set output directories.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../build/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../build/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../build/lib)

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ldl")
endif ()

######## find all the source files.
file(GLOB SOURCE_FILES *.cpp)

######## Set executable file name, and add the source files for it.
add_executable(tap ${SOURCE_FILES} utils.h SolverVPLSHammingSX.cpp SolverVPLSHammingSX.h Solution.h pricingSolver.cpp pricingSolver.h CGTAPInstance.cpp CGTAPInstance.h Query.cpp Query.h JVMAdapter.cpp JVMAdapter.h ActiveDomains.h ActiveDomains.cpp DiversificationInit.cpp DiversificationInit.h RandomInit.h IntensificationInit.cpp IntensificationInit.h randThenSortInit.cpp randThenSortInit.h RandomInitSkew.h KnapsackSolver.cpp KnapsackSolver.h KSInit.h KMppInit.h InitTargets.h princingCPSolver.cpp princingCPSolver.h UserProfile.cpp UserProfile.h FakeTimeStats.cpp FakeTimeStats.h)

######## Add Dependency Library
target_link_libraries(tap Boost::filesystem)
target_link_libraries(tap Boost::boost)
target_link_libraries(tap ilocplex cplex-library cplex-concert cplex-cp)
target_link_libraries(tap cpr::cpr)
target_link_libraries(tap nlohmann_json::nlohmann_json)
#add_dependencies(tap eigen)

####### Fix linker issue on ubuntu
set(CMAKE_CXX_LINK_EXECUTABLE  "${CMAKE_CXX_LINK_EXECUTABLE} -ldl")

###### Verbose ON/OFF
set(CMAKE_VERBOSE_MAKEFILE on)


