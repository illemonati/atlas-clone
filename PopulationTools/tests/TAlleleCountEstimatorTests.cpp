//
// Created by caduffm on 8/19/21.
//
#include "../TAlleleCountEstimator.h"
#include "gtest/gtest.h"

#include <stddef.h>
#include <cstdint>
#include <memory>
#include <vector>

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TRandomGenerator.h"
#include "TSampleLikelihoods.h"
#include "devtools.h"
#include "probability.h"

using namespace testing;
using namespace genometools;

class BridgeTSiteAlleleFrequencyLikelihoods : public PopulationTools::TSiteAlleleFrequencyLikelihoods {
    // bridge to protected functions of TSiteAlleleFrequencyLikelihoods
public:
    BridgeTSiteAlleleFrequencyLikelihoods(int numIndividuals) : PopulationTools::TSiteAlleleFrequencyLikelihoods(numIndividuals){};
    void _fillLog(PopulationTools::TSampleLikelihoods* data, uint32_t numSamples){
        PopulationTools::TSiteAlleleFrequencyLikelihoods::_fillLog(data, numSamples);
    };
    void _fillNatural(PopulationTools::TSampleLikelihoods* data, uint32_t numSamples){
        PopulationTools::TSiteAlleleFrequencyLikelihoods::_fillNatural(data, numSamples);
    };
};

class TSiteAlleleFrequencyLikelihoods_Test : public Test {
public:
    coretools::TRandomGenerator randomGenerator;
    size_t N = 5;

    void fillSampleLikelihoodsOnlyDiploid(genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> * SampleLikelihoods,
                                          const std::vector<uint32_t> & PhredScores){
        // fill phred scores
        size_t sample = 0;
        for (size_t i = 0; i < N*3; i += 3, sample++){
            SampleLikelihoods[sample].setDiploid(HighPrecisionPhredIntProbability(PhredScores[i]),
                                                 HighPrecisionPhredIntProbability(PhredScores[i+1]),
                                                 HighPrecisionPhredIntProbability(PhredScores[i+2]));
        }
    }

    void fillSampleLikelihoodsOnlyHaploid(genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> * SampleLikelihoods,
                                          const std::vector<uint32_t> & PhredScores){
        // fill phred scores
        size_t sample = 0;
        for (size_t i = 0; i < N*2; i += 2, sample++){
            SampleLikelihoods[sample].setHaploid(HighPrecisionPhredIntProbability(PhredScores[i]),
                                                 HighPrecisionPhredIntProbability(PhredScores[i+1]));
        }
    }

    void fillSampleLikelihoodsHaploidDiploid(genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> * SampleLikelihoods,
                                             const std::vector<uint32_t> & PhredScores,
                                             const std::vector<bool> & isDiploid){
        // fill phred scores
        size_t sample = 0;
        for (size_t i = 0; i < N*3; i += 3, sample++){
            if (isDiploid[sample]){
                SampleLikelihoods[sample].setDiploid(HighPrecisionPhredIntProbability(PhredScores[i]),
                                                     HighPrecisionPhredIntProbability(PhredScores[i+1]),
                                                     HighPrecisionPhredIntProbability(PhredScores[i+2]));
            } else {
                SampleLikelihoods[sample].setHaploid(HighPrecisionPhredIntProbability(PhredScores[i]),
                                                     HighPrecisionPhredIntProbability(PhredScores[i+1]));
            }
        }
    }
};

//----------------------------------
// Only diploid samples
//----------------------------------

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyDiploid){
    // simulated genotypes 2,0,0,1,1 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
    std::vector<uint32_t> phredScores = {1000,1000,0,
                                         0,1000,1000,
                                         0,1000,1000,
                                         1000,0,1000,
                                         1000,0,1000};
    genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
    fillSampleLikelihoodsOnlyDiploid(sampleLikelihoods, phredScores);

    // estimate allele counts
    BridgeTSiteAlleleFrequencyLikelihoods estimator(N);
    for (size_t r = 0; r < 2; r++){ // test both _fillNatural and _fillLog
        if (r == 0){
			OUT(r);
            estimator._fillNatural(sampleLikelihoods, N);
			ECHO("");
        } else {
			OUT(r);
            estimator._fillLog(sampleLikelihoods, N);
			ECHO("");
        }
        auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

        // check against results from R
        // 1) log allele frequency likelihoods
        std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {-3.11211114923584,-1.66754787999198,-0.686247040081626,-0.619649775043733,0,-0.874674287024046,-1.36107447309946,-2.27572980058999,-2.83296540780635,-3.92855654614028,-5.41469624222989};
        EXPECT_EQ(estimatedLogAlleleFrequencyLikelihoods.size(), estimatedLogAlleleFrequencyLikelihoods_fromR.size());
        for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++){
            EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
        }

        // 2) ML allele count
        int MLE_AlleleCount = estimator.getMLAlleleCount();
        EXPECT_EQ(MLE_AlleleCount, 4);
        EXPECT_EQ(estimator.getNumAlleles(), 10);
    }
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyDiploid_fixedRef){
    // simulated genotypes 0,0,0,0,0 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
    std::vector<uint32_t> phredScores = {0,1000,1000,
                                         0,1000,1000,
                                         0,1000,1000,
                                         0,1000,1000,
                                         0,1000,1000};
    genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
    fillSampleLikelihoodsOnlyDiploid(sampleLikelihoods, phredScores);

    // estimate allele counts
    BridgeTSiteAlleleFrequencyLikelihoods estimator(N);
    for (size_t r = 0; r < 2; r++){ // test both _fillNatural and _fillLog
        if (r == 0){
            estimator._fillNatural(sampleLikelihoods, N);
        } else {
            estimator._fillLog(sampleLikelihoods, N);
        }
        auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

        // check against results from R
        // 1) log allele frequency likelihoods
        std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {0,-2.30258509299405,-3.91202300542815,-5.52146091786225,-6.82551718074516,-8.10224933935353,-9.12810227373921,-10.1266311038503,-10.8197782844103,-11.5129254649702,-11.5129254649702};
        EXPECT_EQ(estimatedLogAlleleFrequencyLikelihoods.size(), estimatedLogAlleleFrequencyLikelihoods_fromR.size());
        for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++){
            EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
        }

        // 2) ML allele count
        int MLE_AlleleCount = estimator.getMLAlleleCount();
        EXPECT_EQ(MLE_AlleleCount, 0);
        EXPECT_EQ(estimator.getNumAlleles(), 10);
    }
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyDiploid_fixedAlt){
    // simulated genotypes 2,2,2,2,2 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
    std::vector<uint32_t> phredScores = {1000,1000,0,
                                         1000,1000,0,
                                         1000,1000,0,
                                         1000,1000,0,
                                         1000,1000,0};
    genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
    fillSampleLikelihoodsOnlyDiploid(sampleLikelihoods, phredScores);

    // estimate allele counts
    BridgeTSiteAlleleFrequencyLikelihoods estimator(N);
    for (size_t r = 0; r < 2; r++){ // test both _fillNatural and _fillLog
        if (r == 0){
            estimator._fillNatural(sampleLikelihoods, N);
        } else {
            estimator._fillLog(sampleLikelihoods, N);
        }
        auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

        // check against results from R
        // 1) log allele frequency likelihoods
        std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {-11.5129254649702,-11.5129254649702,-10.8197782844103,-10.1266311038503,-9.12810227373921,-8.10224933935353,-6.82551718074516,-5.52146091786225,-3.91202300542815,-2.30258509299405,0};
        EXPECT_EQ(estimatedLogAlleleFrequencyLikelihoods.size(), estimatedLogAlleleFrequencyLikelihoods_fromR.size());
        for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++){
            EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
        }

        // 2) ML allele count
        int MLE_AlleleCount = estimator.getMLAlleleCount();
        EXPECT_EQ(MLE_AlleleCount, 10);
        EXPECT_EQ(estimator.getNumAlleles(), 10);
    }
}

//----------------------------------
// Only haploid samples
//----------------------------------

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyHaploid){
    // simulated genotypes 1,0,0,0,1 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
    std::vector<uint32_t> phredScores = {1000,0,
                                         0,1000,
                                         0,1000,
                                         0,1000,
                                         1000,0};
    genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
    fillSampleLikelihoodsOnlyHaploid(sampleLikelihoods, phredScores);

    // estimate allele counts
    BridgeTSiteAlleleFrequencyLikelihoods estimator(N);
    for (size_t r = 0; r < 2; r++){ // test both _fillNatural and _fillLog
        if (r == 0){
            estimator._fillNatural(sampleLikelihoods, N);
        } else {
            estimator._fillLog(sampleLikelihoods, N);
        }
        auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

        // check against results from R
        // 1) log allele frequency likelihoods
        std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {-2.36113697994366,-0.95995400633002,0,-1.24268938477478,-2.86531806099098,-4.66372207293771};
        for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++){
            EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
        }

        // 2) ML allele count
        int MLE_AlleleCount = estimator.getMLAlleleCount();
        EXPECT_EQ(MLE_AlleleCount, 2);
        EXPECT_EQ(estimator.getNumAlleles(), 5);
    }
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyHaploid_refFixed){
    // simulated genotypes 0,0,0,0,0 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
    std::vector<uint32_t> phredScores = {0,1000,
                                         0,1000,
                                         0,1000,
                                         0,1000,
                                         0,1000};
    genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
    fillSampleLikelihoodsOnlyHaploid(sampleLikelihoods, phredScores);

    // estimate allele counts
    BridgeTSiteAlleleFrequencyLikelihoods estimator(N);
    for (size_t r = 0; r < 2; r++){ // test both _fillNatural and _fillLog
        if (r == 0){
            estimator._fillNatural(sampleLikelihoods, N);
        } else {
            estimator._fillLog(sampleLikelihoods, N);
        }
        auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

        // check against results from R
        // 1) log allele frequency likelihoods
        std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {0,-2.30258509299405,-4.60517018598809,-6.90775527898214,-9.21034037197618,-11.5129254649702};
        for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++){
            EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
        }

        // 2) ML allele count
        int MLE_AlleleCount = estimator.getMLAlleleCount();
        EXPECT_EQ(MLE_AlleleCount, 0);
        EXPECT_EQ(estimator.getNumAlleles(), 5);
    }
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyHaploid_altFixed){
    // simulated genotypes 1,1,1,1,1 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
    std::vector<uint32_t> phredScores = {1000,0,
                                         1000,0,
                                         1000,0,
                                         1000,0,
                                         1000,0};
    genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
    fillSampleLikelihoodsOnlyHaploid(sampleLikelihoods, phredScores);

    // estimate allele counts
    BridgeTSiteAlleleFrequencyLikelihoods estimator(N);
    for (size_t r = 0; r < 2; r++){ // test both _fillNatural and _fillLog
        if (r == 0){
            estimator._fillNatural(sampleLikelihoods, N);
        } else {
            estimator._fillLog(sampleLikelihoods, N);
        }
        auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

        // check against results from R
        // 1) log allele frequency likelihoods
        std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {-11.5129254649702,-9.21034037197618,-6.90775527898214,-4.60517018598809,-2.30258509299405,0};
        for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++){
            EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
        }

        // 2) ML allele count
        int MLE_AlleleCount = estimator.getMLAlleleCount();
        EXPECT_EQ(MLE_AlleleCount, 5);
        EXPECT_EQ(estimator.getNumAlleles(), 5);
    }
}

//----------------------------------
// Mixed diploid and haploid samples
//----------------------------------

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, haploidDiploid){
    // simulate diploid, haploid, diploid, haploid, diploid
    // simulated genotypes 1,0,0,0,1 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
    std::vector<uint32_t> phredScores = {0,1000,1000,
                                         0,1000,1000,
                                         1000,0,1000,
                                         1000,0,1000,
                                         1000,1000,0};
    genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
    fillSampleLikelihoodsHaploidDiploid(sampleLikelihoods, phredScores, {true, false, true, false, true});

    // estimate allele counts
    BridgeTSiteAlleleFrequencyLikelihoods estimator(N);
    for (size_t r = 0; r < 2; r++){ // test both _fillNatural and _fillLog
        if (r == 0){
            estimator._fillNatural(sampleLikelihoods, N);
        } else {
            estimator._fillLog(sampleLikelihoods, N);
        }
        auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

        // check against results from R
        // 1) log allele frequency likelihoods
        std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {-3.48187285612952,-2.08624716758074,-1.17671963502976,-0.826181680011293,0,-0.826181680011293,-1.17671963502976,-2.08624716758074,-3.48187285612952};
        for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++){
            EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
        }

        // 2) ML allele count
        int MLE_AlleleCount = estimator.getMLAlleleCount();
        EXPECT_EQ(MLE_AlleleCount, 4);
        EXPECT_EQ(estimator.getNumAlleles(), 8);
    }
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, haploidDiploid_RefFixed){
    // simulate diploid, haploid, diploid, haploid, diploid
    // simulated genotypes 0,0,0,0,0 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
    std::vector<uint32_t> phredScores = {0,1000,1000,
                                         0,1000,1000,
                                         0,1000,1000,
                                         0,1000,1000,
                                         0,1000,1000};
    genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
    fillSampleLikelihoodsHaploidDiploid(sampleLikelihoods, phredScores, {true, false, true, false, true});

    // estimate allele counts
    BridgeTSiteAlleleFrequencyLikelihoods estimator(N);
    for (size_t r = 0; r < 2; r++){ // test both _fillNatural and _fillLog
        if (r == 0){
            estimator._fillNatural(sampleLikelihoods, N);
        } else {
            estimator._fillLog(sampleLikelihoods, N);
        }
        auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

        // check against results from R
        // 1) log allele frequency likelihoods
        std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {0,-2.30258509299405,-3.93004151093082,-5.5486119069282,-6.88236604497732,-8.18072095479502,-9.27678547138433,-10.3342704686286,-11.5129254649702};
        for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++){
            EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
        }

        // 2) ML allele count
        int MLE_AlleleCount = estimator.getMLAlleleCount();
        EXPECT_EQ(MLE_AlleleCount, 0);
        EXPECT_EQ(estimator.getNumAlleles(), 8);
    }
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, haploidDiploid_AltFixed){
    // simulate diploid, haploid, diploid, haploid, diploid
    // simulated genotypes 2,1,2,1,2 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
    std::vector<uint32_t> phredScores = {1000,1000,0,
                                         1000,0,1000,
                                         1000,1000,0,
                                         1000,0,1000,
                                         1000,1000,0};
    genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
    fillSampleLikelihoodsHaploidDiploid(sampleLikelihoods, phredScores, {true, false, true, false, true});

    // estimate allele counts
    BridgeTSiteAlleleFrequencyLikelihoods estimator(N);
    for (size_t r = 0; r < 2; r++){ // test both _fillNatural and _fillLog
        if (r == 0){
            estimator._fillNatural(sampleLikelihoods, N);
        } else {
            estimator._fillLog(sampleLikelihoods, N);
        }
        auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

        // check against results from R
        // 1) log allele frequency likelihoods
        std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {-11.5129254649702,-10.3342704686286,-9.27678547138433,-8.18072095479502,-6.88236604497732,-5.5486119069282,-3.93004151093082,-2.30258509299405,0};
        for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++){
            EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
        }

        // 2) ML allele count
        int MLE_AlleleCount = estimator.getMLAlleleCount();
        EXPECT_EQ(MLE_AlleleCount, 8);
        EXPECT_EQ(estimator.getNumAlleles(), 8);
    }
}
