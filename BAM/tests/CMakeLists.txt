#[[
cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(BAMTestsLib CXX C)

add_library(${PROJECT_NAME}
        TBed_Test.cpp
        TGenomePosition_Test.cpp
        TBamFile_Test.cpp
        )

include_directories("../..")
include_directories("${gtest_SOURCE_DIR}/include"
        "${gmock_SOURCE_DIR}/include"
        )
include_directories(../../bamtools/)
include_directories(../../BAM/)
include_directories(../../commonutilities/)
include_directories(../../GenomeTasks/)
include_directories(../../GenotypeLikelihoods/)
include_directories(../../GLF/)
include_directories(../../PopulationTools/)
include_directories(../../Simulations/)
include_directories(../../TestUtilities/)
include_directories(../../VCF/)
]]
