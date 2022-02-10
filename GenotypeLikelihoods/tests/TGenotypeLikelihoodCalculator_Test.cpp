//
// Created by linkv on 8/17/20.
//

#include "gtest/gtest.h"

#include "GenotypeTypes.h"
#include "TGenotypeLikelihoodCalculator.h"

using namespace GenotypeLikelihoods;
using genometools::Base;
using genometools::Genotype;

TEST(TGenotypeLikelihoodCalculator_test, calculateGenotypeLikelihoods_emptySite){
    TSite site;
    TGenotypeLikelihoods genotypeLikelihoods;
    TGenotypeLikelihoodCalculator calculator;

    calculator.calculateGenotypeLikelihoods(site, genotypeLikelihoods);

    EXPECT_EQ(genotypeLikelihoods[Genotype::AA],1);
    EXPECT_EQ(genotypeLikelihoods[Genotype::AC],1);
    EXPECT_EQ(genotypeLikelihoods[Genotype::AG],1);
    EXPECT_EQ(genotypeLikelihoods[Genotype::AT],1);
    EXPECT_EQ(genotypeLikelihoods[Genotype::CC],1);
    EXPECT_EQ(genotypeLikelihoods[Genotype::CG],1);
    EXPECT_EQ(genotypeLikelihoods[Genotype::CT],1);
    EXPECT_EQ(genotypeLikelihoods[Genotype::GG],1);
    EXPECT_EQ(genotypeLikelihoods[Genotype::GG],1);
    EXPECT_EQ(genotypeLikelihoods[Genotype::TT],1);
};

TEST(TGenotypeLikelihoodCalculator_test, calculateGenotypeLikelihoods_noPMDnoRecal){
    TSite site;
    TGenotypeLikelihoods genotypeLikelihoods;
    TGenotypeLikelihoodCalculator calculator;

    BAM::TSequencedBase base;
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    base.base = Base::A;
    site.add(base);

    calculator.calculateGenotypeLikelihoods(site, genotypeLikelihoods);

    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::AA],oneMinusError);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::AC],0.5 - errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::AG], (oneMinusError + errorOneThird) / 2.0);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::AT],0.5 - errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::CC], errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::CG],(2.0  * errorOneThird) / 2.0);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::CT], errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::GG], errorOneThird);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::GT],(2.0 * errorOneThird) / 2.0);
    EXPECT_FLOAT_EQ(genotypeLikelihoods[Genotype::TT], errorOneThird);
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
