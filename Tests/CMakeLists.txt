cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(TestsLib CXX C)

add_library(${PROJECT_NAME}
        TAtlasTest.cpp TAtlasTest.h
        TAtlasTestFilter.cpp TAtlasTestFilter.h
        TAtlasTestMergePairs.cpp TAtlasTestMergePairs.h
        TAtlasTestPMD.cpp TAtlasTestPMD.h
        TAtlasTestRecalibration.cpp TAtlasTestRecalibration.h
        TVcfTest.cpp TVcfTest.h
        )

include_directories(".")
include_directories(../BAM/)
include_directories(../bamtools/)
include_directories(../commonutilities/)
include_directories(../GenomeTasks/)
include_directories(../GenotypeLikelihoods/)
include_directories(../GLF/)
include_directories(../PopulationTools/)
include_directories(../Simulations/)
include_directories(../VCF/)