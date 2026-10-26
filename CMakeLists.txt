cmake_minimum_required(VERSION 3.17)
project(cplex_test)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

if (MSVC)
    add_compile_definitions(_HAS_STD_BYTE=0)
endif ()

######## Set output directories.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../build/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../build/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../build/lib)

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
        GIT_TAG b3755347973a23567e1f786615ebb480e1f704c4) # The commit hash for 1.9.0. Replace with the latest from: https://github.com/libcpr/cpr/releases
FetchContent_MakeAvailable(cpr)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz)
FetchContent_MakeAvailable(json)

set(BOOST_URL "https://archives.boost.io/release/1.71.0/source/boost_1_71_0.tar.bz2" CACHE STRING "Boost download URL" FORCE)
add_subdirectory(boost-cmake)

######## Complier message
SET(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
message(STATUS "System: ${CMAKE_SYSTEM}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ldl")
endif ()

######## find all the source files.
file(GLOB SOURCE_FILES *.cpp)

######## Set executable file name, and add the source files for it.
add_executable(tap ${SOURCE_FILES} utils.h SolverVPLSHammingSX.cpp SolverVPLSHammingSX.h Solution.h pricingSolver.cpp pricingSolver.h CGTAPInstance.cpp CGTAPInstance.h Query.cpp Query.h JVMAdapter.cpp JVMAdapter.h ActiveDomains.h ActiveDomains.cpp RandomInit.h randThenSortInit.cpp randThenSortInit.h KnapsackSolver.cpp KnapsackSolver.h KSInit.h KMppInit.h InitTargets.h princingCPSolver.cpp princingCPSolver.h UserProfile.cpp UserProfile.h DiversificationInit.h IntensificationInit.h
        CheckInsight.h)

######## Add Dependency Library
target_link_libraries(tap Boost::filesystem)
target_link_libraries(tap Boost::boost)
target_link_libraries(tap ilocplex cplex-library cplex-concert cplex-cp)
target_link_libraries(tap cpr::cpr)
target_link_libraries(tap nlohmann_json::nlohmann_json)
#add_dependencies(tap eigen)

####### Copy required runtime DLLs on Windows (cpr.dll, libcurl.dll, cplex2220.dll)
if (WIN32)
    add_custom_command(TARGET tap POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E $<IF:$<BOOL:$<TARGET_RUNTIME_DLLS:tap>>,copy_if_different,true>
            $<TARGET_RUNTIME_DLLS:tap>
            $<TARGET_FILE_DIR:tap>
        COMMAND_EXPAND_LISTS
    )
    if (CPLEX_DLL)
        add_custom_command(TARGET tap POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CPLEX_DLL}"
                $<TARGET_FILE_DIR:tap>
        )
    endif ()
endif ()

####### Fix linker issue on ubuntu
if (UNIX)
    set(CMAKE_CXX_LINK_EXECUTABLE  "${CMAKE_CXX_LINK_EXECUTABLE} -ldl")
endif ()

###### Verbose ON/OFF
set(CMAKE_VERBOSE_MAKEFILE on)

message("Sortie bin: ${CMAKE_BINARY_DIR}")
