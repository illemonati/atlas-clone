cmake_minimum_required(VERSION 3.4.1 FATAL_ERROR)

project(GenotypeLikelihoodsLib CXX C)

add_library(${PROJECT_NAME}
        auxiliaryTools.cpp auxiliaryTools.h
        TGenotypeData.cpp TGenotypeData.h
        TGenotypeDistribution.cpp TGenotypeDistribution.h
        TGenotypeLikelihoodCalculator.cpp TGenotypeLikelihoodCalculator.h
        TGenotypeMap.h TGenotypePrior.h
        TPostMortemDamage.cpp TPostMortemDamage.h
        TRecalibrationEMEstimator.cpp TRecalibrationEMEstimator.h
        TSequencingErrorCovariate.cpp TSequencingErrorCovariate.h
        TSequencingErrorCovariateFunction.cpp TSequencingErrorCovariateFunction.h
        TSequencingErrorModel.cpp TSequencingErrorModel.h
        TSequencingErrorModels.cpp TSequencingErrorModels.h
        TSite.cpp TSite.h
        TSiteSubset.cpp TSiteSubset.h
        TThetaEstimator.cpp TThetaEstimator.h
        TThetaEstimatorData.cpp TThetaEstimatorData.h
        TWindow.cpp TWindow.h
        )

include_directories(".")
include_directories(../BAM/)
include_directories(../bamtools/)
include_directories(../commonutilities/)
include_directories(../GenomeTasks/)
include_directories(../GLF/)
include_directories(../PopulationTools/)
include_directories(../Simulations/)
include_directories(../TestUtilities/)
include_directories(../VCF/)
