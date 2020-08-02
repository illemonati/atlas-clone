cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(VCFLib CXX C)

add_library(${PROJECT_NAME}
        TVcfCompare.cpp TVcfCompare.h
        TVcfConverter.cpp TVcfConverter.h
        TVcfDiagnostics.cpp TVcfDiagnostics.h
        TVCFFields.h TVcfFile.cpp
        TVcfFile.h
        TVcfParser.cpp TVcfParser.h
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
include_directories(../Tests/)
