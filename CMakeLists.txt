# Project
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(atlas VERSION 1.0.0)

cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# compile commands: by default, compile in Release mode
# type 'cmake -DCMAKE_BUILD_TYPE=Debug ..' in terminal for compiling in Debug mode
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# locate external libraries
find_package(ZLIB REQUIRED)
find_package(Armadillo REQUIRED)

# Download coretools and stattools
include(FetchContent)
FetchContent_Declare(coretools
        GIT_REPOSITORY    https://bitbucket.org/wegmannlab/coretools.git
        GIT_TAG           b58485467705901fde55ab39b70f0068025bbc3d
        SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/coretools"
        )
FetchContent_MakeAvailable(coretools)

FetchContent_Declare(stattools
        GIT_REPOSITORY    https://bitbucket.org/wegmannlab/stattools.git
        GIT_TAG           e1660d0aabccc1c1af409612aa447b13a8256dad
        SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/stattools"
        )
FetchContent_MakeAvailable(stattools)

FetchContent_Declare(genometools
        GIT_REPOSITORY    https://bitbucket.org/wegmannlab/genometools.git
        GIT_TAG           8c528ec36755b592fa6f65b2bacb74cceaea9a83
        SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/genometools"
        )
FetchContent_MakeAvailable(genometools)

set(LINK_LIBRARIES
        coreLib
        coretools
        genometools
        stattools
        ${ZLIB_LIBRARIES}
        ${ARMADILLO_LIBRARIES}
        )

include_directories(
        BAM
        bamtools
        GenomeTasks
        GenotypeLikelihoods
        GLF
        PopulationTools
        Simulations
        TestUtilities
        VCF
        coretools/core
        coretools/core/Types
        coretools/core/IntegrationTests
        genometools/core        
        stattools/core        
        ${ARMADILLO_INCLUDE_DIRS}
)

#######################
# Create core library
#######################
# assemble list of cpp and h files in all source folders
file(GLOB_RECURSE BAMTOOLS bamtools/*)
file(GLOB CORE_SOURCES BAM/*.h BAM/*.cpp
                               GenomeTasks/*.h GenomeTasks/*.cpp
                               GenotypeLikelihoods/*.h GenotypeLikelihoods/*.cpp
                               GLF/*.h GLF/*.cpp
                               PopulationTools/*.h PopulationTools/*.cpp
                               Simulations/*.h Simulations/*.cpp
                               TestUtilities/*.h TestUtilities/*.cpp
                               VCF/*.h VCF/*.cpp
        )
list(APPEND CORE_SOURCES ${BAMTOOLS})

#Core library
add_library(coreLib STATIC ${CORE_SOURCES})

#######################
# Executable
#######################

# add executable to be built
add_executable(${PROJECT_NAME}
        atlas.cpp
        )

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC
        ${LINK_LIBRARIES}
        )

# add the binary tree to the search path for include files
target_include_directories(${PROJECT_NAME} PUBLIC)

#######################
# Tests
#######################
#assemble list of test cpp and h files
file(GLOB_RECURSE TEST_SOURCES
        BAM/tests/*.h BAM/tests/*.cpp
        bamtools/tests/*
        GenomeTasks/tests/*.h GenomeTasks/tests/*.cpp
        GenotypeLikelihoods/tests/*.h GenotypeLikelihoods/tests/*.cpp
        GLF/tests/*.h GLF/tests/*.cpp
        PopulationTools/tests/*.h PopulationTools/tests/*.cpp
        Simulations/tests/*.h Simulations/tests/*.cpp
        TestUtilities/tests/*.h TestUtilities/tests/*.cpp
        VCF/tests/*.h VCF/tests/*.cpp
        )

add_executable("${PROJECT_NAME}_unitTests" ${TEST_SOURCES})

target_link_libraries("${PROJECT_NAME}_unitTests" PUBLIC
        gtest_main
        gmock_main
        ${LINK_LIBRARIES}
        )

target_include_directories("${PROJECT_NAME}_unitTests" PUBLIC
        "${PROJECT_BINARY_DIR}" #google tests will be cloned in there
        BAM/tests
        bamtools/tests
        GenomeTasks/tests
        GenotypeLikelihoods/tests
        GLF/tests
        PopulationTools/tests
        Simulations/tests
        TestUtilities/tests
        VCF/tests
        coretools/tests
        genometools/tests
        stattools/tests
        )
