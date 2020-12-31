# Project
project(atlas)

cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_CXX_STANDARD 14)

# add libraries
add_subdirectory(commonutilities)

# locate external libraries
find_package(ZLIB REQUIRED)
find_package(Armadillo REQUIRED)

set(LINK_LIBRARIES
        coreLib
        commonutilities
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
        commonutilities/core
        commonutilities/core/IntegrationTests
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

# compile commands: by default, compile in Release mode
# type 'cmake -DCMAKE_BUILD_TYPE=Debug ..' in terminal for compiling in Debug mode
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

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
        commonutilities/tests
        )