//
// Created by linkv on 8/17/20.
//

#include "TGenotypeLikelihoodCalculator.h"
#include "TGenotypeMap.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using namespace GenotypeLikelihoods;

TEST(TGenotypeLikelihoodCalculator_test, calculateGenotypeLikelihoods_emptySite){
    TSite site;
    TGenotypeLikelihoods genotypeLikelihoods;
    TGenotypeLikelihoodCalculator calculator;

    calculator.calculateGenotypeLikelihoods(site, genotypeLikelihoods);

    EXPECT_EQ(genotypeLikelihoods[AA],1);
    EXPECT_EQ(genotypeLikelihoods[AC],1);
    EXPECT_EQ(genotypeLikelihoods[AG],1);
    EXPECT_EQ(genotypeLikelihoods[AT],1);
    EXPECT_EQ(genotypeLikelihoods[CC],1);
    EXPECT_EQ(genotypeLikelihoods[CG],1);
    EXPECT_EQ(genotypeLikelihoods[CT],1);
    EXPECT_EQ(genotypeLikelihoods[GG],1);
    EXPECT_EQ(genotypeLikelihoods[GG],1);
    EXPECT_EQ(genotypeLikelihoods[TT],1);
};

TEST(TGenotypeLikelihoodCalculator_test, calculateGenotypeLikelihoods_noPMDnoRecal){
    TSite site;
    TGenotypeLikelihoods genotypeLikelihoods;
    TGenotypeLikelihoodCalculator calculator;

    BAM::TBase base;
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    base.base = A;
    site.add(base);

    calculator.calculateGenotypeLikelihoods(site, genotypeLikelihoods);

    EXPECT_FLOAT_EQ(genotypeLikelihoods[AA],oneMinusError);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[AC],0.5 - errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[AG], (oneMinusError + errorOneThird) / 2.0);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[AT],0.5 - errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[CC], errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[CG],(2.0  * errorOneThird) / 2.0);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[CT], errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[GG], errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[GT],(2.0 * errorOneThird) / 2.0);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[TT], errorOneThird);
}

//vec[AA] *= oneMinusError;
//vec[CC] *= errorOneThird;
//vec[GG] *= (1.0 - PMD_GA) * errorOneThird + PMD_GA * oneMinusError;
//vec[TT] *= errorOneThird;
//vec[AC] *= 0.5 - errorOneThird;
//vec[AG] *= ((1.0 + PMD_GA) * oneMinusError + (1.0 - PMD_GA) * errorOneThird) / 2.0;
//vec[AT] *= 0.5 - errorOneThird;
//vec[CT] *= errorOneThird;
//double tmp = (PMD_GA * oneMinusError + (2.0 - PMD_GA) * errorOneThird) / 2.0;
//vec[CG] *= tmp;
//vec[GT] *= tmp;


TEST(TGenotypeLikelihoodCalculator_test, calculateGenotypeLikelihoods_PMD){

}

TEST(TGenotypeLikelihoodCalculator_test, calculateGenotypeLikelihoods_recal){

}