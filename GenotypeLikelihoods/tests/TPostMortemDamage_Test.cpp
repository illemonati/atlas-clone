//
// Created by vivian on 18.08.20.
//

#include "TPostMortemDamage.h"
#include <gtest/gtest.h>

#include <stddef.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TGenotypeData.h"
#include "TLog.h"
#include "TPMDTables.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "TSequencingErrorModels.h"
#include "gtest/gtest.h"
#include "probability.h"
#include "weakTypes.h"

using namespace GenotypeLikelihoods;
using genometools::Base;

TEST(TPostMortemDamage_test, PMDTable) {
	using namespace genometools;
	TPMDTable t1(10);
	EXPECT_EQ(t1.size(), 10);
	t1.resize(100);
	EXPECT_EQ(t1.size(), 100);

	EXPECT_EQ(t1[Base::G][Base::C][99], 0);
	for (size_t _ = 0; _ < 3; ++_) t1.add(99, Base::G, Base::C);

	t1.add(99, Base::G, Base::A);
	t1.add(99, Base::G, Base::G);
	t1.add(99, Base::G, Base::T);
	EXPECT_EQ(t1[Base::G][Base::C][99], 3);
	EXPECT_EQ(t1[Base::G][Base::A][99], 1);
	EXPECT_EQ(t1[Base::G][Base::G][99], 1);
	EXPECT_EQ(t1[Base::G][Base::T][99], 1);
	EXPECT_EQ(t1.sums(Base::G)[99], 6);
	t1.add(101, Base::A, Base::A);
	t1.add(102, Base::A, Base::A);
	t1.add(103, Base::A, Base::A);
	t1.add(104, Base::A, Base::A);
	t1.add(105, Base::A, Base::A);
	EXPECT_EQ(t1[Base::A][Base::A][t1.size()], 5);

	TPMDTable t2(100);
	EXPECT_EQ(t2.size(), 100);
	t2.add(99, Base::G, Base::C);
	t2.add(0,Base::C, Base::T);
	EXPECT_EQ(t2[Base::G][Base::C][99], 1);
	EXPECT_EQ(t2[Base::C][Base::T][0], 1);
	EXPECT_EQ(t2.sums(Base::G)[99], 1);
	EXPECT_EQ(t2.sums(Base::C)[0], 1);

	t1.add(t2);
	EXPECT_EQ(t1[Base::G][Base::C][99], 4);
	EXPECT_EQ(t1[Base::G][Base::A][99], 1);
	EXPECT_EQ(t1[Base::G][Base::G][99], 1);
	EXPECT_EQ(t1[Base::G][Base::T][99], 1);
	EXPECT_EQ(t1[Base::C][Base::T][0], 1);
	EXPECT_EQ(t1.sums(Base::G)[99], 7);
	EXPECT_EQ(t1.sums(Base::C)[0], 1);

	t1.empty();
	EXPECT_EQ(t1[Base::G][Base::C][99], 0);
	EXPECT_EQ(t1.sums(Base::G)[99], 0);
}

TEST(TPostMortemDamage_test, PMDTables) {
	using namespace genometools;
	BAM::TReadGroups rg;
	rg.add("1");
	BAM::TReadGroupMap rgm(rg);
	TPMDTables ts(&rg, 10, &rgm);
	EXPECT_EQ(ts[0][forward3].size(), 10);

	BAM::TSequencedBase sqbase;
	sqbase.base           = Base::A;
	sqbase.distFrom3Prime = 1;
	sqbase.distFrom5Prime = 33;
	sqbase.readGroupID    = 0;
	sqbase.setReverseStrand(false);

	// forward3
	EXPECT_EQ(ts[0][forward3].sums(Base::C)[sqbase.distFrom3Prime], 0);
	EXPECT_EQ(ts[0][forward3][Base::C][sqbase.base][sqbase.distFrom3Prime], 0);
	ts.add(sqbase, Base::C);
	EXPECT_EQ(ts[0][forward3].sums(Base::C)[sqbase.distFrom3Prime], 1);
	EXPECT_EQ(ts[0][forward3][Base::C][sqbase.base][sqbase.distFrom3Prime], 1);

	sqbase.base = Base::C;
	ts.add(sqbase, Base::C);
	EXPECT_EQ(ts[0][forward3].sums(Base::C)[sqbase.distFrom3Prime], 2);
	EXPECT_EQ(ts[0][forward3][Base::C][sqbase.base][sqbase.distFrom3Prime], 1);

	// forward5
	EXPECT_EQ(ts[0][forward5].sums(Base::A)[sqbase.distFrom3Prime], 0);
	EXPECT_EQ(ts[0][forward5][Base::A][sqbase.base][sqbase.distFrom3Prime], 0);
	sqbase.distFrom5Prime = sqbase.distFrom3Prime;
	ts.add(sqbase, Base::A);
	EXPECT_EQ(ts[0][forward5].sums(Base::A)[sqbase.distFrom3Prime], 1);
	EXPECT_EQ(ts[0][forward5][Base::A][sqbase.base][sqbase.distFrom3Prime], 1);

	// reverse3
	EXPECT_EQ(ts[0][reverse3].sums(Base::C)[sqbase.distFrom3Prime], 0);
	EXPECT_EQ(ts[0][reverse3][Base::C][sqbase.base][sqbase.distFrom3Prime], 0);
	sqbase.distFrom5Prime++;
	sqbase.setReverseStrand(true);
	ts.add(sqbase, Base::G);
	EXPECT_EQ(ts[0][reverse3].sums(Base::C)[sqbase.distFrom3Prime], 1);
	EXPECT_EQ(ts[0][reverse3][Base::C][flipped(sqbase.base)][sqbase.distFrom3Prime], 1);

	// reverse5
	EXPECT_EQ(ts[0][reverse5].sums(Base::A)[sqbase.distFrom3Prime], 0);
	EXPECT_EQ(ts[0][reverse5][Base::A][sqbase.base][sqbase.distFrom3Prime], 0);
	sqbase.distFrom5Prime--;
	sqbase.setReverseStrand(true);
	ts.add(sqbase, Base::T);
	EXPECT_EQ(ts[0][reverse5].sums(Base::A)[sqbase.distFrom3Prime], 1);
	EXPECT_EQ(ts[0][reverse5][Base::A][flipped(sqbase.base)][sqbase.distFrom3Prime], 1);
}

TEST(TPostMortemDamage_test, noPMD) {
	EXPECT_ANY_THROW(TPMDFunctionNoPMD("[3, 4]"));
	const TPMDFunctionNoPMD fn("[]");
	EXPECT_FALSE(fn.hasDamage());
	EXPECT_FLOAT_EQ(fn.prob(0), 0.);
	EXPECT_FLOAT_EQ(fn.prob(122), 0.);
	EXPECT_FLOAT_EQ(fn.prob(static_cast<uint16_t>(-1)), 0.);
}

TEST(TPostMortemDamage_test, exponential) {
	EXPECT_ANY_THROW(TPMDFunctionExponential("[1]"));
	EXPECT_ANY_THROW(TPMDFunctionExponential("[3, 4]"));
	EXPECT_ANY_THROW(TPMDFunctionExponential("[3, 4, 5]"));

	const TPMDFunctionExponential fn0("[]");
	EXPECT_FALSE(fn0.hasDamage());
	EXPECT_FLOAT_EQ(fn0.prob(0), 0.);
	EXPECT_FLOAT_EQ(fn0.prob(33), 0.);
	EXPECT_FLOAT_EQ(fn0.prob(static_cast<uint16_t>(-1)), 0.);

	constexpr auto N = 8;
	constexpr auto a = 1.;
	constexpr auto b = 1.;
	constexpr auto c = 1.;
	const TPMDFunctionExponential fn1("[8, 1., 1., 1.]");
	for (size_t p = 0; p < N + 1; ++p) {
		EXPECT_TRUE(fn1.hasDamage());
		EXPECT_FLOAT_EQ(fn1.prob(p), a*std::exp(-b * p) + c);
	}
	EXPECT_FLOAT_EQ(fn1.prob(N + 1), a*std::exp(-b * N) + c);
	EXPECT_FLOAT_EQ(fn1.prob(N + 13), a*std::exp(-b * N) + c);
}

TEST(TPostMortemDamage_test, empiric) {
	EXPECT_ANY_THROW(TPMDFunctionEmpiric("[1.3]"));
	EXPECT_ANY_THROW(TPMDFunctionEmpiric("[0.5, 0.3, 3, 0.4]"));
	EXPECT_ANY_THROW(TPMDFunctionEmpiric("[-0.3]"));

	const TPMDFunctionEmpiric fn0("[]");
	EXPECT_FALSE(fn0.hasDamage());
	EXPECT_FLOAT_EQ(fn0.prob(0), 0.);
	EXPECT_FLOAT_EQ(fn0.prob(33), 0.);
	EXPECT_FLOAT_EQ(fn0.prob(static_cast<uint16_t>(-1)), 0.);

	std::vector<double> params(10);
	std::iota(params.begin(), params.end(), 0);
	std::reverse(params.begin(), params.end());
	const double fac = params.front()*2;
	std::string s = "[";
	for (auto &p : params) {
		p /= fac;
		s += std::to_string(p) + ",";
	}
	s += "]";
	const TPMDFunctionEmpiric fn1(s);
	EXPECT_TRUE(fn1.hasDamage());
	for (size_t p = 0; p < params.size(); ++p) {
		EXPECT_NEAR(fn1.prob(p), params[p], 1e-5);
	}
	EXPECT_NEAR(fn1.prob(113), params.back(), 1e-5);
}

TEST(TPostMortemDamage_test, empiric_learn) {
	using namespace genometools;
	constexpr size_t N = 10;
	TPMDTable t1(N);

	for (size_t i = 0; i < N; ++i) {
		for (size_t _ = 0; _ < N - i; ++_) t1.add(i, Base::G, Base::A);
		for (size_t _ = 0; _ < i; ++_) t1.add(i, Base::G, Base::G);

		EXPECT_EQ(t1[Base::G][Base::A][i], N - i);
		EXPECT_EQ(t1.sums(Base::G)[i], 10);
	}

	TPMDFunctionEmpiric fne("[]");
	fne.learn(t1, Base::G, Base::A, TPMDEstimationParameters{});
	EXPECT_EQ(fne.string(), "Empiric[1.000000,0.900000,0.800000,0.700000,0.600000,0.500000,0.400000,0.300000,0.200000,0.100000]");
}

TEST(TPostMortemDamage_test, exp_learn) {
	using namespace genometools;

	constexpr auto a    = 0.5;
	constexpr auto b    = 0.7;
	constexpr auto c    = 0.3;
	constexpr size_t N  = 20;
	constexpr size_t Nb = 1000;

	TPMDTable t1(N);

	for (size_t i = 0; i < N; ++i) {
		const size_t nPMD = Nb*(a*std::exp(-b*i) + c);
		for (size_t _ = 0; _ < nPMD; ++_) t1.add(i, Base::G, Base::A);
		for (size_t _ = 0; _ < Nb - nPMD; ++_) t1.add(i, Base::G, Base::G);
	}

	TPMDEstimationParameters eparams;
	eparams.emplace(TPMDFunctionExponential::epsilon, 0.001);
	eparams.emplace(TPMDFunctionExponential::numNR, 1000);
	TPMDFunctionExponential fne("[]");
	fne.learn(t1, Base::G, Base::A, eparams);
	EXPECT_EQ(fne.string(), "Exponential[19,0.500266,0.700669,0.299745]");
}

TEST(TPostMortemDamage_test, baseANoPMD) {
		constexpr auto err = 0.01;

		SequencingError::TModels sem;
		TBaseLikelihoods sem_likelihoods;
		TPostMortemDamage pmd;
		TBaseLikelihoods pmd_likelihoods;

		BAM::TSequencedBase base;
		base.originalQuality_phredInt = 20;
		base.base                     = Base::A;

		sem.fillBaseLikelihoods(base, sem_likelihoods);
		pmd.fillBaseLikelihoods(base, sem_likelihoods, pmd_likelihoods);

		for (Base b = Base::min; b < Base::max; ++b) {
			base.base = b;
			sem.fillBaseLikelihoods(base, sem_likelihoods);
			pmd.fillBaseLikelihoods(base, sem_likelihoods, pmd_likelihoods);

			for (Base trueBase = Base::min; trueBase < Base::max; ++trueBase) {
				if (trueBase == b) {
					EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], 1. - err);
					EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase], 1. - err);
				} else {
					EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], err / 3);
					EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase], err / 3);
				}
			}
		}
}

TEST(TPostMortemDamage_test, baseAWithPMD) {
	using namespace genometools;
	constexpr auto err = 0.01;

	SequencingError::TModels sem;
	TBaseLikelihoods sem_likelihoods;
	TBaseLikelihoods pmd_likelihoods;

	BAM::TSequencedBase base;
	base.originalQuality_phredInt = 20;
	base.base                     = Base::A;
	base.readGroupID              = 0;
	base.distFrom5Prime           = 2;
	base.setReverseStrand(false);

	// initialize RG
	BAM::TReadGroups ReadGroups;
	ReadGroups.add("testRG");

	// initialize PMD
	coretools::TLog logfile;
	std::vector<uint16_t> ReadGroupsWithoutPMD;
	TPostMortemDamage pmd("doubleStrand:Empiric[0.3]:Empiric[0.1]", ReadGroups, ReadGroupsWithoutPMD);

	for (uint16_t dfrom3 = 0; dfrom3 < 3; dfrom3 += 2) {
		base.distFrom3Prime = dfrom3;
		const auto from3    = base.distFrom3Prime < base.distFrom5Prime;
		for (Base b = Base::min; b < Base::max; ++b) {
			base.base = b;

			sem.fillBaseLikelihoods(base, sem_likelihoods);
			pmd.fillBaseLikelihoods(base, sem_likelihoods, pmd_likelihoods);

			for (Base trueBase = Base::min; trueBase < Base::max; ++trueBase) {

				if (trueBase == b) {
					EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], 1. - err);
					if (from3 && trueBase == Base::G) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
								(1.0 - 0.1) * sem_likelihoods[Base::G] + 0.1 * sem_likelihoods[Base::A]);
					} else if (!from3 && trueBase == Base::C) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
								(1.0 - 0.3) * sem_likelihoods[Base::C] + 0.3 * sem_likelihoods[Base::T]);
					} else {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase], 1. - err);
					}
				} else {
					EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], err/3);
					if (from3 && trueBase == Base::G) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
								(1.0 - 0.1) * sem_likelihoods[Base::G] + 0.1 * sem_likelihoods[Base::A]);
					} else if (!from3 && trueBase == Base::C) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
								(1.0 - 0.3) * sem_likelihoods[Base::C] + 0.3 * sem_likelihoods[Base::T]);
					} else {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase], err/3);
					}
				}
			}
		}
	}
}
