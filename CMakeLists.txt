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
  GenotypeLikelihoods/SequencingError/*.h GenotypeLikelihoods/SequencingError/*.cpp
  GenotypeLikelihoods/PMD/*.h GenotypeLikelihoods/PMD/*.cpp
  GLF/*.h GLF/*.cpp
  PopulationTools/*.h PopulationTools/*.cpp
  Simulations/*.h Simulations/*.cpp
  VCF/*.h VCF/*.cpp)

include(FetchContent)
FetchContent_Declare(coretools
  GIT_REPOSITORY    https://bitbucket.org/wegmannlab/coretools.git
  GIT_TAG           4c04b9e
  SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/coretools"
)

include(FetchContent)
FetchContent_Declare(genometools
  GIT_REPOSITORY    https://bitbucket.org/wegmannlab/genometools.git
  GIT_TAG           dc258ac
  SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/genometools"
)

include(FetchContent)
FetchContent_Declare(stattools
  GIT_REPOSITORY    https://bitbucket.org/wegmannlab/stattools.git
  GIT_TAG           81631d1
  SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/stattools"
)

FetchContent_Declare(bamtools
  GIT_REPOSITORY    https://github.com/glis-glis/bamtools.git
  GIT_TAG           7c3abb4
  SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/bamtools"
)

FetchContent_MakeAvailable(coretools genometools stattools bamtools)

target_sources(${PROJECT_NAME} PRIVATE ${CORE_SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE 
  ${CMAKE_CURRENT_BINARY_DIR}/coretools/core/IntegrationTests
  BAM
  GenomeTasks
  GenotypeLikelihoods
  GLF
  PopulationTools
  Simulations
  VCF)

target_link_libraries(${PROJECT_NAME} PRIVATE BamTools BamTools-utils stattools)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra ${SAN})
target_compile_definitions(${PROJECT_NAME} PRIVATE WITH_ARMADILLO ARMA_WARN_LEVEL=0 DEVTOOLS)

add_executable("${PROJECT_NAME}_unitTests" EXCLUDE_FROM_ALL)

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS 
  BAM/tests/*.h BAM/tests/*.cpp
  bamtools/tests/*
  GenomeTasks/tests/*.h GenomeTasks/tests/*.cpp
  GenotypeLikelihoods/tests/*.h GenotypeLikelihoods/tests/*.cpp
  GLF/tests/*.h GLF/tests/*.cpp
  PopulationTools/tests/*.h PopulationTools/tests/*.cpp
  Simulations/tests/*.h Simulations/tests/*.cpp
  VCF/tests/*.h VCF/tests/*.cpp)

target_sources("${PROJECT_NAME}_unitTests" PRIVATE ${CORE_SOURCES} ${TEST_SOURCES})

target_include_directories("${PROJECT_NAME}_unitTests" PRIVATE 
  BAM
  GenomeTasks
  GenotypeLikelihoods
  GLF
  PopulationTools
  Simulations
  VCF)

target_link_libraries("${PROJECT_NAME}_unitTests" PRIVATE BamTools BamTools-utils stattools gtest_main gmock_main)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options("${PROJECT_NAME}_unitTests" PRIVATE -Wall -Wextra)
target_compile_definitions("${PROJECT_NAME}_unitTests" PRIVATE WITH_ARMADILLO ARMA_WARN_LEVEL=0 DEVTOOLS)
