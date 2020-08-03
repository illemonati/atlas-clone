cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(SimulationsLib CXX C)

add_library(${PROJECT_NAME}
        SFS.cpp SFS.h
        TSimulatorAuxiliaryTools.cpp TSimulatorAuxiliaryTools.h
        TSimulator.cpp TSimulator.h
        TSimulatorQualityTransformation.cpp TSimulatorQualityTransformation.h
        TSimulatorRead.cpp TSimulatorRead.h
        TSimulatorReadLength.cpp TSimulatorReadLength.h
        )

include_directories(".")
include_directories(../BAM/)
include_directories(../bamtools/)
include_directories(../commonutilities/)
include_directories(../GenomeTasks/)
include_directories(../GenotypeLikelihoods/)
include_directories(../GLF/)
include_directories(../PopulationTools/)
include_directories(../Tests/)
include_directories(../VCF/)