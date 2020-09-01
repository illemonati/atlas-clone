cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(PopulationToolsLib CXX C)

add_library(${PROJECT_NAME}
        TAlleleCountEstimator.cpp TAlleleCountEstimator.h
        TAlleleCountFileFormat.cpp TAlleleCountFileFormat.h
        TAlleleFrequencyEstimator.cpp TAlleleFrequencyEstimator.h
        TDistanceEstimator.cpp TDistanceEstimator.h
        TGenotypeFrequencies.cpp TGenotypeFrequencies.h
        THardyWeinbergTest.cpp THardyWeinbergTest.h
        TInbreedingEstimator.cpp TInbreedingEstimator.h
        TMajorMinor.cpp TMajorMinor.h
        TPolymorhicWindowIdentifier.cpp TPolymorhicWindowIdentifier.h
        TPopulationLikelihoodLocus.cpp TPopulationLikelihoodLocus.h
        TPopulationLikelihoods.cpp TPopulationLikelihoods.h
        )

include_directories(".")
include_directories(../BAM/)
include_directories(../bamtools/)
include_directories(../commonutilities/)
include_directories(../GenomeTasks/)
include_directories(../GenotypeLikelihoods/)
include_directories(../GLF/)
include_directories(../Simulations/)
include_directories(../TestUtilities/)
include_directories(../VCF/)

