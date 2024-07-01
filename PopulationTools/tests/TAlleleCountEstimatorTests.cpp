//
// Created by caduffm on 8/19/21.
//
#include "TAlleleCountEstimator.h"
#include "coretools/Main/TRandomGenerator.h"
#include "gtest/gtest.h"


using namespace testing;
using namespace genometools;
using namespace PopulationTools;
using coretools::logP;

class TSiteAlleleFrequencyLikelihoods_Test : public Test {
public:
    coretools::TRandomGenerator randomGenerator;
    static constexpr size_t N = 5;

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

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyDiploid) {
	// simulated genotypes 2,0,0,1,1 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
	std::vector<uint32_t> phredScores = {1000, 1000, 0, 0, 1000, 1000, 0, 1000, 1000, 1000, 0, 1000, 1000, 0, 1000};
	genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
	fillSampleLikelihoodsOnlyDiploid(sampleLikelihoods, phredScores);

	// estimate allele counts
	TSiteAlleleFrequencyLikelihoods estimator;
	estimator.fill(sampleLikelihoods, N);
	auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();
	const auto max =
		*std::max_element(estimatedLogAlleleFrequencyLikelihoods.begin(), estimatedLogAlleleFrequencyLikelihoods.end());

	// check against results from R
	// 1) log allele frequency likelihoods
	std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {
		logP(-3.11211114923584), logP(-1.66754787999198), logP(-0.686247040081626), logP(-0.619649775043733), logP(0),
		logP(-0.874674287024046),      logP(-1.36107447309946), logP(-2.27572980058999),  logP(-2.83296540780635),  logP(-3.92855654614028),
		logP(-5.41469624222989)};
	EXPECT_EQ(estimatedLogAlleleFrequencyLikelihoods.size(), estimatedLogAlleleFrequencyLikelihoods_fromR.size());
	for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++) {
		EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i].scale(max),
						estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
	}

	// 2) ML allele count
	int MLE_AlleleCount = estimator.MLAlleleCount();
	EXPECT_EQ(MLE_AlleleCount, 4);
	EXPECT_EQ(estimator.Nalleles(), 10);
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyDiploid_fixedRef) {
	// simulated genotypes 0,0,0,0,0 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
	std::vector<uint32_t> phredScores = {0, 1000, 1000, 0, 1000, 1000, 0, 1000, 1000, 0, 1000, 1000, 0, 1000, 1000};
	genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
	fillSampleLikelihoodsOnlyDiploid(sampleLikelihoods, phredScores);

	// estimate allele counts
	TSiteAlleleFrequencyLikelihoods estimator;
	estimator.fill(sampleLikelihoods, N);
	auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

	// check against results from R
	// 1) log allele frequency likelihoods
	std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {logP(0),
		logP(-2.30258509299405),
		logP(-3.91202300542815),
		logP(-5.52146091786225),
		logP(-6.82551718074516),
		logP(-8.10224933935353),
		logP(-9.12810227373921),
		logP(-10.1266311038503),
		logP(-10.8197782844103),
		logP(-11.5129254649702),
		logP(-11.5129254649702)};
	EXPECT_EQ(estimatedLogAlleleFrequencyLikelihoods.size(), estimatedLogAlleleFrequencyLikelihoods_fromR.size());
	for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++) {
		EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
	}

	// 2) ML allele count
	int MLE_AlleleCount = estimator.MLAlleleCount();
	EXPECT_EQ(MLE_AlleleCount, 0);
	EXPECT_EQ(estimator.Nalleles(), 10);
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyDiploid_fixedAlt) {
	// simulated genotypes 2,2,2,2,2 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
	std::vector<uint32_t> phredScores = {1000, 1000, 0, 1000, 1000, 0, 1000, 1000, 0, 1000, 1000, 0, 1000, 1000, 0};
	genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
	fillSampleLikelihoodsOnlyDiploid(sampleLikelihoods, phredScores);

	// estimate allele counts
	TSiteAlleleFrequencyLikelihoods estimator;
	estimator.fill(sampleLikelihoods, N);
	auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

	// check against results from R
	// 1) log allele frequency likelihoods
	std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {logP(-11.5129254649702),
																						   logP(-11.5129254649702),
																						   logP(-10.8197782844103),
																						   logP(-10.1266311038503),
																						   logP(-9.12810227373921),
																						   logP(-8.10224933935353),
																						   logP(-6.82551718074516),
																						   logP(-5.52146091786225),
																						   logP(-3.91202300542815),
																						   logP(-2.30258509299405),
																						   logP(0)};
	EXPECT_EQ(estimatedLogAlleleFrequencyLikelihoods.size(), estimatedLogAlleleFrequencyLikelihoods_fromR.size());
	for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++) {
		EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
	}

	// 2) ML allele count
	int MLE_AlleleCount = estimator.MLAlleleCount();
	EXPECT_EQ(MLE_AlleleCount, 10);
	EXPECT_EQ(estimator.Nalleles(), 10);
}

//----------------------------------
// Only haploid samples
//----------------------------------

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyHaploid) {
	// simulated genotypes 1,0,0,0,1 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
	std::vector<uint32_t> phredScores = {1000, 0, 0, 1000, 0, 1000, 0, 1000, 1000, 0};
	genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
	fillSampleLikelihoodsOnlyHaploid(sampleLikelihoods, phredScores);

	// estimate allele counts
	TSiteAlleleFrequencyLikelihoods estimator;
	estimator.fill(sampleLikelihoods, N);
	auto lkh = estimator.getLogAlleleFrequencyLikelihoods();
	const auto max = *std::max_element(lkh.begin(), lkh.end());

	// check against results from R
	// 1) log allele frequency likelihoods
	std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {
		logP(-2.36113697994366), logP(-0.95995400633002), logP(0), logP(-1.24268938477478), logP(-2.86531806099098), logP(-4.66372207293771)};
	for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++) {
		EXPECT_FLOAT_EQ(lkh[i].scale(max), estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
	}

	// 2) ML allele count
	int MLE_AlleleCount = estimator.MLAlleleCount();
	EXPECT_EQ(MLE_AlleleCount, 2);
	EXPECT_EQ(estimator.Nalleles(), 5);
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyHaploid_refFixed) {
	// simulated genotypes 0,0,0,0,0 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
	std::vector<uint32_t> phredScores = {0, 1000, 0, 1000, 0, 1000, 0, 1000, 0, 1000};
	genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
	fillSampleLikelihoodsOnlyHaploid(sampleLikelihoods, phredScores);

	// estimate allele counts
	TSiteAlleleFrequencyLikelihoods estimator;
	estimator.fill(sampleLikelihoods, N);
	auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

	// check against results from R
	// 1) log allele frequency likelihoods
	std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {
		logP(0), logP(-2.30258509299405), logP(-4.60517018598809), logP(-6.90775527898214), logP(-9.21034037197618), logP(-11.5129254649702)};
	for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++) {
		EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
	}

	// 2) ML allele count
	int MLE_AlleleCount = estimator.MLAlleleCount();
	EXPECT_EQ(MLE_AlleleCount, 0);
	EXPECT_EQ(estimator.Nalleles(), 5);
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, onlyHaploid_altFixed) {
	// simulated genotypes 1,1,1,1,1 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
	std::vector<uint32_t> phredScores = {1000, 0, 1000, 0, 1000, 0, 1000, 0, 1000, 0};
	genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
	fillSampleLikelihoodsOnlyHaploid(sampleLikelihoods, phredScores);

	// estimate allele counts
	TSiteAlleleFrequencyLikelihoods estimator;
	estimator.fill(sampleLikelihoods, N);
	auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

	// check against results from R
	// 1) log allele frequency likelihoods
	std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {
		logP(-11.5129254649702), logP(-9.21034037197618), logP(-6.90775527898214), logP(-4.60517018598809), logP(-2.30258509299405), logP(0)};
	for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++) {
		EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
	}

	// 2) ML allele count
	int MLE_AlleleCount = estimator.MLAlleleCount();
	EXPECT_EQ(MLE_AlleleCount, 5);
	EXPECT_EQ(estimator.Nalleles(), 5);
}

//----------------------------------
// Mixed diploid and haploid samples
//----------------------------------

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, haploidDiploid) {
	// simulate diploid, haploid, diploid, haploid, diploid
	// simulated genotypes 1,0,0,0,1 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
	std::vector<uint32_t> phredScores = {0, 1000, 1000, 0, 1000, 1000, 1000, 0, 1000, 1000, 0, 1000, 1000, 1000, 0};
	genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
	fillSampleLikelihoodsHaploidDiploid(sampleLikelihoods, phredScores, {true, false, true, false, true});

	// estimate allele counts
	TSiteAlleleFrequencyLikelihoods estimator;
	estimator.fill(sampleLikelihoods, N);
	auto lkh = estimator.getLogAlleleFrequencyLikelihoods();
	const auto max = *std::max_element(lkh.begin(), lkh.end());

	// check against results from R
	// 1) log allele frequency likelihoods
	std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {
		logP(-3.48187285612952),  logP(-2.08624716758074), logP(-1.17671963502976), logP(-0.826181680011293), logP(0),
		logP(-0.826181680011293), logP(-1.17671963502976), logP(-2.08624716758074), logP(-3.48187285612952)};
	for (size_t i = 0; i < lkh.size(); i++) {
		EXPECT_FLOAT_EQ(lkh[i].scale(max), estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
	}

	// 2) ML allele count
	int MLE_AlleleCount = estimator.MLAlleleCount();
	EXPECT_EQ(MLE_AlleleCount, 4);
	EXPECT_EQ(estimator.Nalleles(), 8);
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, haploidDiploid_RefFixed) {
	// simulate diploid, haploid, diploid, haploid, diploid
	// simulated genotypes 0,0,0,0,0 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
	std::vector<uint32_t> phredScores = {0, 1000, 1000, 0, 1000, 1000, 0, 1000, 1000, 0, 1000, 1000, 0, 1000, 1000};
	genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
	fillSampleLikelihoodsHaploidDiploid(sampleLikelihoods, phredScores, {true, false, true, false, true});

	// estimate allele counts
	TSiteAlleleFrequencyLikelihoods estimator;
	estimator.fill(sampleLikelihoods, N);
	auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

	// check against results from R
	// 1) log allele frequency likelihoods
	std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {logP(0),
		logP(-2.30258509299405),
		logP(-3.93004151093082),
		logP(-5.5486119069282),
		logP(-6.88236604497732),
		logP(-8.18072095479502),
		logP(-9.27678547138433),
		logP(-10.3342704686286),
		logP(-11.5129254649702)};
	for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++) {
		EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
	}

	// 2) ML allele count
	int MLE_AlleleCount = estimator.MLAlleleCount();
	EXPECT_EQ(MLE_AlleleCount, 0);
	EXPECT_EQ(estimator.Nalleles(), 8);
}

TEST_F(TSiteAlleleFrequencyLikelihoods_Test, haploidDiploid_AltFixed) {
	// simulate diploid, haploid, diploid, haploid, diploid
	// simulated genotypes 2,1,2,1,2 in RScript Dropbox/PhD/atlas/TAlleleCountsUnitTests.R
	std::vector<uint32_t> phredScores = {1000, 1000, 0, 1000, 0, 1000, 1000, 1000, 0, 1000, 0, 1000, 1000, 1000, 0};
	genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> sampleLikelihoods[N];
	fillSampleLikelihoodsHaploidDiploid(sampleLikelihoods, phredScores, {true, false, true, false, true});

	// estimate allele counts
	TSiteAlleleFrequencyLikelihoods estimator;
	estimator.fill(sampleLikelihoods, N);
	auto estimatedLogAlleleFrequencyLikelihoods = estimator.getLogAlleleFrequencyLikelihoods();

	// check against results from R
	// 1) log allele frequency likelihoods
	std::vector<coretools::LogProbability> estimatedLogAlleleFrequencyLikelihoods_fromR = {
	    logP(-11.5129254649702), logP(-10.3342704686286), logP(-9.27678547138433),
	    logP(-8.18072095479502), logP(-6.88236604497732), logP(-5.5486119069282),
	    logP(-3.93004151093082), logP(-2.30258509299405), logP(0)};
	for (size_t i = 0; i < estimatedLogAlleleFrequencyLikelihoods_fromR.size(); i++) {
		EXPECT_FLOAT_EQ(estimatedLogAlleleFrequencyLikelihoods[i], estimatedLogAlleleFrequencyLikelihoods_fromR[i]);
	}

	// 2) ML allele count
	int MLE_AlleleCount = estimator.MLAlleleCount();
	EXPECT_EQ(MLE_AlleleCount, 8);
	EXPECT_EQ(estimator.Nalleles(), 8);
}
