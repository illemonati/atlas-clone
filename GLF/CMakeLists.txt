
cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(GLFLib CXX C)

add_library(${PROJECT_NAME}
        TGLF.cpp TGLF.h
        TGlfMultiReader.cpp TGlfMultiReader.h
        )

include_directories(".")
include_directories(../BAM/)
include_directories(../bamtools/)
include_directories(../commonutilities/)
include_directories(../GenomeTasks/)
include_directories(../GenotypeLikelihoods/)
include_directories(../PopulationTools/)
include_directories(../Simulations/)
include_directories(../Tests/)
include_directories(../VCF/)
