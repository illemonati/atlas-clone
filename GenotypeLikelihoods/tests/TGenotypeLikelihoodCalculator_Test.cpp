//
// Created by linkv on 8/17/20.
//
#include "TErrorModels.h"
#include "gtest/gtest.h"

#include <memory>

#include "TReadGroupInfo.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TGenotypeData.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "TSite.h"
#include "coretools/Types/probability.h"

using namespace GenotypeLikelihoods;
using genometools::Base;
using genometools::Genotype;

TEST(TGenotypeLikelihoodCalculator_test, calculateGenotypeLikelihoods_emptySite){
    TSite site;
	BAM::RGInfo::TReadGroupInfo rgi;
    TErrorModels calculator(rgi);

	const auto genotypeLikelihoods = calculator.calculateGenotypeLikelihoods(site);

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
	BAM::TReadGroups rg;
	rg.add("test");
	BAM::RGInfo::TReadGroupInfo rgi(rg);
    TErrorModels calculator(rgi);

    BAM::TSequencedBase base;
    base.quality_orig = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    base.base = Base::A;
	base.readGroupID = 0;
    site.add(base);

    const auto genotypeLikelihoods = calculator.calculateGenotypeLikelihoods(site);

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
