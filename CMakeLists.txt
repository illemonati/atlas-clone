cmake_minimum_required(VERSION 3.14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Asan")
  set(SAN "-O1;-g;-fsanitize=address;-fno-omit-frame-pointer;-DNDEBUG")
  set(SAN_LINK "-fsanitize=address")
  message("Compiling and linking with adress sanitizer")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Msan")
  set(SAN "-O1;-g;-fsanitize=memory;-fno-omit-frame-pointer;-fsanitize-memory-track-origins")
  set(SAN_LINK "-fsanitize=memory")
  message("Compiling and linking with memory sanitizer")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "UBsan")
  set(SAN "-O1;-g;-fsanitize=undefined;-fno-omit-frame-pointer")
  set(SAN_LINK "-fsanitize=undefined")
  message("Compiling and linking with undefined behavior sanitizer")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Perf")
  set(SAN "-O3;-DNDEBUG;-fno-omit-frame-pointer")
  set(SAN_LINK "")
  message("Compiling and linking for performance measurements")
else()
  set(SAN "")
  set(SAN_LINK "")
endif()

# Project
project(atlas VERSION 1.0.0 LANGUAGES CXX)

# Atlas
add_executable(${PROJECT_NAME} atlas.cpp)

file(GLOB CORE_SOURCES CONFIGURE_DEPENDS
  BAM/*.h BAM/*.cpp
  GenomeTasks/*.h GenomeTasks/*.cpp
  GenotypeLikelihoods/*.h GenotypeLikelihoods/*.cpp
  GenotypeLikelihoods/SequencingError/*.h GenotypeLikelihoods/SequencingError/*.cpp
  GLF/*.h GLF/*.cpp
  PopulationTools/*.h PopulationTools/*.cpp
  Simulations/*.h Simulations/*.cpp
  VCF/*.h VCF/*.cpp)

find_package(ZLIB REQUIRED)

# Download coretools and stattools
include(FetchContent)
FetchContent_Declare(coretools
        GIT_REPOSITORY    https://bitbucket.org/wegmannlab/coretools.git
        GIT_TAG           17af274
        SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/coretools"
        )

FetchContent_Declare(stattools
        GIT_REPOSITORY    https://bitbucket.org/wegmannlab/stattools.git
        GIT_TAG           fcfa7ff
        SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/stattools"
        )

FetchContent_Declare(genometools
        GIT_REPOSITORY    https://bitbucket.org/wegmannlab/genometools.git
        GIT_TAG			  d0167fd
        SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/genometools"
        )

FetchContent_Declare(bamtools
        GIT_REPOSITORY    https://github.com/glis-glis/bamtools.git
        GIT_TAG           master
        SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/bamtools"
        )

FetchContent_Declare(armadillo
		  GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code.git
		  GIT_TAG        11.2.x
		)

FetchContent_MakeAvailable(coretools stattools genometools bamtools armadillo)

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

target_link_libraries(${PROJECT_NAME} PRIVATE BamTools BamTools-utils coretools stattools genometools ${ZLIB_LIBRARIES} armadillo ${SAN_LINK})
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

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        main
)
FetchContent_MakeAvailable(googletest)

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
        VCF)

target_link_libraries("${PROJECT_NAME}_unitTests" PRIVATE BamTools BamTools-utils coretools stattools genometools ${ZLIB_LIBRARIES} armadillo gtest_main gmock_main ${SAN_LINK})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options("${PROJECT_NAME}_unitTests" PRIVATE -Wall -Wextra ${SAN})
