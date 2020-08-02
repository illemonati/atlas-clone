cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(BAMLib CXX C)

add_library(${PROJECT_NAME}
        TAlignment.cpp TAlignment.h
        TBamFile.cpp TBamFile.h
        TBamFilter.cpp TBamFilter.h
        TBase.cpp TBase.h
        TBed.cpp TBed.h
        TBedReaderWindows.cpp TBedReaderWindows.h
        TChromosomes.cpp TChromosomes.h
        TCigar.cpp TCigar.h
        TFastaBuffer.cpp TFastaBuffer.h
        TGenomePosition.cpp TGenomePosition.h
        TQualityMap.cpp TQualityMap.h
        TReadGroups.cpp TReadGroups.h
        TSamFlags.h
        TSamHeader.cpp TSamHeader.h
        )

include_directories(".")
include_directories(../bamtools/)
include_directories(../commonutilities/)
include_directories(../GenomeTasks/)
include_directories(../GenotypeLikelihoods/)
include_directories(../GLF/)
include_directories(../PopulationTools/)
include_directories(../Simulations/)
include_directories(../Tests/)
include_directories(../VCF/)