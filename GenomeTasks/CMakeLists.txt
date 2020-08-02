cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(GenomeTasksLib CXX C)

add_library(${PROJECT_NAME}
        TAllelicDepthCounts.cpp TAllelicDepthCounts.h
        TBamDiagnoser.cpp TBamDiagnoser.h
        TBamDownsampler.cpp TBamDownsampler.h
        TBamFilter.cpp TBamFilter.h
        TCaller.cpp TCaller.h
        TContextQuantifer.cpp TContextQuantifer.h
        TCreateBedMask.cpp TCreateBedMask.h
        TDepthWriter.cpp TDepthWriter.h
        TDuplicateQuantifier.cpp TDuplicateQuantifier.h
        TEstimateRecalibration.cpp TEstimateRecalibration.h
        TEstimateTheta.cpp TEstimateTheta.h
        TGenome.cpp TGenome.h
        TPileup.cpp TPileup.h
        TPMDEstimator.cpp TPMDEstimator.h
        TPMDSCalculator.cpp TPMDSCalculator.h
        TPSMCInput.cpp TPSMCInput.h
        TQualityDistribution.cpp TQualityDistribution.h
        TReadGroupMerger.cpp TReadGroupMerger.h
        TSoftClipping.cpp TSoftClipping.h
        TWriteGLF.cpp TWriteGLF.h
        )

include_directories(".")
include_directories(../BAM/)
include_directories(../bamtools/)
include_directories(../commonutilities/)
include_directories(../GenotypeLikelihoods/)
include_directories(../GLF/)
include_directories(../PopulationTools/)
include_directories(../Simulations/)
include_directories(../Tests/)
include_directories(../VCF/)