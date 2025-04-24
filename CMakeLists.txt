cmake_minimum_required(VERSION 3.14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Don't warn about missing policies!
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif (POLICY CMP0135)

# Project
project(atlas VERSION 2.0.0 LANGUAGES CXX)

# Atlas
add_executable(${PROJECT_NAME} atlas.cpp)

file(GLOB CORE_SOURCES CONFIGURE_DEPENDS
  BAM/*.h BAM/*.cpp
  GenomeTasks/*.h GenomeTasks/*.cpp
  GenotypeLikelihoods/*.h GenotypeLikelihoods/*.cpp
  GenotypeLikelihoods/SequencingError/*.h GenotypeLikelihoods/SequencingError/*.cpp
  GenotypeLikelihoods/PMD/*.h GenotypeLikelihoods/PMD/*.cpp
  GenotypeLikelihoods/oldPMD/*.cpp
  GLF/*.h GLF/*.cpp
  PopulationTools/*.h PopulationTools/*.cpp
  Simulations/*.h Simulations/*.cpp
  VCF/*.h VCF/*.cpp)

include(FetchContent)
FetchContent_Declare(coretools
  URL        https://bitbucket.org/wegmannlab/coretools/get/b9df2aa57f68ef84a8795d258a391007a0b4bed2.tar.gz
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/coretools"
)

FetchContent_Declare(genometools
  URL        https://bitbucket.org/wegmannlab/genometools/get/a54a0952572656680af4da6e5c416f1208d03570.tar.gz
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/genometools"
)

FetchContent_Declare(stattools
  URL        https://bitbucket.org/wegmannlab/stattools/get/fb4f5232198ea712b6dc9b31fe5057041dd102c6.tar.gz
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/stattools"
)

FetchContent_Declare(bamtools
  URL        https://bitbucket.org/wegmannlab/bamtools/get/7c3abb4e716e3646ad0f50f1ef5d04b13b3f2877.tar.gz
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bamtools"
)

FetchContent_MakeAvailable(coretools genometools stattools bamtools)

# htslib
if(CONDA)
  find_package(nlohmann_json REQUIRED)
  find_library(HTSLIB NAMES hts)
else()
  include(htslib.cmake)
  FetchContent_Declare(json
	URL        https://github.com/nlohmann/json/archive/refs/tags/v3.11.3.tar.gz
	SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/json"
  )
  FetchContent_MakeAvailable(json)
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

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

target_link_libraries(${PROJECT_NAME} PRIVATE coretools genometools stattools BamTools BamTools-utils hts nlohmann_json::nlohmann_json)
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra ${SAN})
target_compile_definitions(${PROJECT_NAME} PRIVATE WITH_ARMADILLO ARMA_WARN_LEVEL=0 DEVTOOLS)

# Test
if(BUILD_TESTING)
  add_executable("${PROJECT_NAME}_unitTests")

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

  target_link_libraries("${PROJECT_NAME}_unitTests" PRIVATE coretools genometools stattools BamTools BamTools-utils hts nlohmann_json::nlohmann_json gtest_main )
  target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
  set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
  target_compile_options("${PROJECT_NAME}_unitTests" PRIVATE -Wall -Wextra)
  target_compile_definitions("${PROJECT_NAME}_unitTests" PRIVATE WITH_ARMADILLO ARMA_WARN_LEVEL=0 DEVTOOLS)
endif()
