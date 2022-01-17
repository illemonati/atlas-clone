//
// Created by vivian on 18.08.20.
//

#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TestCase.h"
#include <algorithm>
#include <cstdint>
#include <gtest/gtest.h>

using namespace GenotypeLikelihoods;
using genometools::Base;

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
	EXPECT_TRUE(fn0.hasDamage());
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
	EXPECT_TRUE(fn0.hasDamage());
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

TEST(TPostMortemDamage_test, baseANoPMD) {
	constexpr auto err = 0.01;

	TSequencingErrorModels sem;
	TBaseLikelihoods sem_likelihoods;
	TPostMortemDamage pmd;
	TBaseLikelihoods pmd_likelihoods;

	BAM::TSequencedBase base;
	base.originalQuality_phredInt = 20;
	base.base                     = genometools::A;

	sem.fillBaseLikelihoods(base, sem_likelihoods);
	pmd.fillBaseLikelihoods(base, sem_likelihoods, pmd_likelihoods);

	for (Base b = Base::min(); b < Base::max(); ++b) {
		base.base = b;
		sem.fillBaseLikelihoods(base, sem_likelihoods);
		pmd.fillBaseLikelihoods(base, sem_likelihoods, pmd_likelihoods);

		for (Base trueBase = Base::min(); trueBase < Base::max(); ++trueBase) {
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

	TSequencingErrorModels sem;
	TBaseLikelihoods sem_likelihoods;
	TBaseLikelihoods pmd_likelihoods;

	BAM::TSequencedBase base;
	base.originalQuality_phredInt = 20;
	base.base                     = A;
	base.readGroupID              = 0;
	base.distFrom5Prime           = 2;
	base.setReverseStrand(false);

	// initialize RG
	BAM::TReadGroups ReadGroups;
	ReadGroups.add("testRG");

	// initialize PMD
	TLog logfile;
	std::vector<uint16_t> ReadGroupsWithoutPMD;
	TPostMortemDamage pmd("doubleStrand:Empiric[0.3]:Empiric[0.1]", ReadGroups, &logfile, ReadGroupsWithoutPMD);

	for (uint16_t dfrom3 = 0; dfrom3 < 3; dfrom3 += 2) {
		base.distFrom3Prime = dfrom3;
		const auto from3    = base.distFrom3Prime < base.distFrom5Prime;
		for (Base b = Base::min(); b < Base::max(); ++b) {
			base.base = b;

			sem.fillBaseLikelihoods(base, sem_likelihoods);
			pmd.fillBaseLikelihoods(base, sem_likelihoods, pmd_likelihoods);

			for (Base trueBase = Base::min(); trueBase < Base::max(); ++trueBase) {

				if (trueBase == b) {
					EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], 1. - err);
					if (from3 && trueBase == G) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
								(1.0 - 0.1) * sem_likelihoods[G] + 0.1 * sem_likelihoods[A]);
					} else if (!from3 && trueBase == C) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
								(1.0 - 0.3) * sem_likelihoods[C] + 0.3 * sem_likelihoods[T]);
					} else {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase], 1. - err);
					}
				} else {
					EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], err/3);
					if (from3 && trueBase == G) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
								(1.0 - 0.1) * sem_likelihoods[G] + 0.1 * sem_likelihoods[A]);
					} else if (!from3 && trueBase == C) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
								(1.0 - 0.3) * sem_likelihoods[C] + 0.3 * sem_likelihoods[T]);
					} else {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase], err/3);
					}
				}
			}
		}
	}
}
