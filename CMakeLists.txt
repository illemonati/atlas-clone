cmake_minimum_required(VERSION 3.14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

# Project
project(atlas VERSION 1.0.0 LANGUAGES CXX)

# Atlas
add_executable(${PROJECT_NAME} atlas.cpp)

file(GLOB CORE_SOURCES CONFIGURE_DEPENDS
  BAM/*.h BAM/*.cpp
  GenomeTasks/*.h GenomeTasks/*.cpp
  GenotypeLikelihoods/*.h GenotypeLikelihoods/*.cpp
  GLF/*.h GLF/*.cpp
  PopulationTools/*.h PopulationTools/*.cpp
  Simulations/*.h Simulations/*.cpp
  #  TestUtilities/*.h TestUtilities/*.cpp
  VCF/*.h VCF/*.cpp)

find_package(ZLIB REQUIRED)
find_package(Armadillo REQUIRED)

# Download coretools and stattools
include(FetchContent)
FetchContent_Declare(coretools
        GIT_REPOSITORY    https://bitbucket.org/wegmannlab/coretools.git
        GIT_TAG           856ee71
        SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/coretools"
        )

FetchContent_Declare(stattools
        GIT_REPOSITORY    https://bitbucket.org/wegmannlab/stattools.git
        GIT_TAG           9e38e5d
        SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/stattools"
        )

FetchContent_Declare(genometools
        GIT_REPOSITORY    https://bitbucket.org/wegmannlab/genometools.git
        GIT_TAG           8c46607
        SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/genometools"
        )

FetchContent_Declare(bamtools
        GIT_REPOSITORY    https://github.com/pezmaster31/bamtools.git
        GIT_TAG           master
        SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/bamtools"
        )

FetchContent_MakeAvailable(coretools stattools genometools bamtools)

target_sources(${PROJECT_NAME} PRIVATE ${CORE_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE 
  ${CMAKE_CURRENT_BINARY_DIR}/coretools/core/IntegrationTests
  BAM
  GenomeTasks
  GenotypeLikelihoods
  GLF
  PopulationTools
  Simulations
#  TestUtilities
  VCF)

target_link_libraries(${PROJECT_NAME} PRIVATE BamTools BamTools-utils coretools stattools genometools ${ZLIB_LIBRARIES} ${ARMADILLO_LIBRARIES})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)

add_executable("${PROJECT_NAME}_unitTests" EXCLUDE_FROM_ALL)

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS 
        BAM/tests/*.h BAM/tests/*.cpp
        bamtools/tests/*
        GenomeTasks/tests/*.h GenomeTasks/tests/*.cpp
        GenotypeLikelihoods/tests/*.h GenotypeLikelihoods/tests/*.cpp
        GLF/tests/*.h GLF/tests/*.cpp
        PopulationTools/tests/*.h PopulationTools/tests/*.cpp
        Simulations/tests/*.h Simulations/tests/*.cpp
#        TestUtilities/tests/*.h TestUtilities/tests/*.cpp
        VCF/tests/*.h VCF/tests/*.cpp)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        main
)
#FetchContent_MakeAvailable(googletest)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
option(gtest_disable_pthreads "Disable uses of pthreads in gtest." ON)

target_sources("${PROJECT_NAME}_unitTests" PRIVATE ${CORE_SOURCES} ${TEST_SOURCES})

target_include_directories("${PROJECT_NAME}_unitTests" PRIVATE 
  coretools/core/IntegrationTests
  coretools/tests
        BAM
        GenomeTasks
        GenotypeLikelihoods
        GLF
        PopulationTools
        Simulations
#        TestUtilities
        VCF)

target_link_libraries("${PROJECT_NAME}_unitTests" PRIVATE BamTools BamTools-utils coretools stattools genometools ${ZLIB_LIBRARIES} ${ARMADILLO_LIBRARIES} gtest_main gmock_main)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options("${PROJECT_NAME}_unitTests" PRIVATE -Wall -Wextra)
