//
// Created by vivian on 18.08.20.
//

#include "PMD/TFunction.h"
#include "PMD/TModel.h"
#include "PMD/TModels.h"
#include <gtest/gtest.h>

#include <stddef.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TGenotypeData.h"
#include "coretools/Main/TLog.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"
#include "SequencingError/TModels.h"
#include "gtest/gtest.h"
#include "coretools/Types/probability.h"
#include "coretools/Types/weakTypes.h"

using namespace GenotypeLikelihoods;
using genometools::Base;


TEST(TPostMortemDamage_test, noPMD) {
	EXPECT_ANY_THROW(PMD::TNo("[3, 4]"));
	const PMD::TNo fn("[]");
	EXPECT_FALSE(fn.hasDamage());
	EXPECT_FLOAT_EQ(fn.prob(0), 0.);
	EXPECT_FLOAT_EQ(fn.prob(122), 0.);
	EXPECT_FLOAT_EQ(fn.prob(static_cast<uint16_t>(-1)), 0.);
}

TEST(TPostMortemDamage_test, exponential) {
	EXPECT_ANY_THROW(PMD::TExponential("[1]"));
	EXPECT_ANY_THROW(PMD::TExponential("[3, 4]"));
	EXPECT_ANY_THROW(PMD::TExponential("[3, 4, 5]"));

	const PMD::TExponential fn0("[]");
	EXPECT_FALSE(fn0.hasDamage());
	EXPECT_FLOAT_EQ(fn0.prob(0), 0.);
	EXPECT_FLOAT_EQ(fn0.prob(33), 0.);
	EXPECT_FLOAT_EQ(fn0.prob(static_cast<uint16_t>(-1)), 0.);

	constexpr auto N = 8;
	constexpr auto a = .3;
	constexpr auto b = .2;
	constexpr auto c = .1;
	const PMD::TExponential fn1("[8, .3, .2, .1]");
	for (size_t p = 0; p < N + 1; ++p) {
		EXPECT_TRUE(fn1.hasDamage());
		EXPECT_FLOAT_EQ(fn1.prob(p), a*std::exp(-b * p) + c);
	}
	EXPECT_FLOAT_EQ(fn1.prob(N + 1), a*std::exp(-b * N) + c);
	EXPECT_FLOAT_EQ(fn1.prob(N + 13), a*std::exp(-b * N) + c);
}

TEST(TPostMortemDamage_test, empiric) {
	EXPECT_ANY_THROW(PMD::TEmpiric("[1.3]"));
	EXPECT_ANY_THROW(PMD::TEmpiric("[0.5, 0.3, 3, 0.4]"));
	EXPECT_ANY_THROW(PMD::TEmpiric("[-0.3]"));

	const PMD::TEmpiric fn0("[]");
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
	const PMD::TEmpiric fn1(s);
	EXPECT_TRUE(fn1.hasDamage());
	for (size_t p = 0; p < params.size(); ++p) {
		EXPECT_NEAR(fn1.prob(p), params[p], 1e-5);
	}
	EXPECT_NEAR(fn1.prob(113), params.back(), 1e-5);
}


TEST(TPostMortemDamage_test, baseANoPMD) {
	constexpr auto err = 0.01;

	// initialize RG
	BAM::TReadGroups rgs;
	rgs.add("testRG");

	SequencingError::TModels sem;
	PMD::TModels pmd;

	sem.initializeNoRecal(rgs);

	BAM::TSequencedBase base;
	base.originalQuality_phredInt = 20;
	base.readGroupID = 0;

	for (Base b = Base::min; b < Base::max; ++b) {
		base.base                  = b;
		const auto sem_likelihoods = sem.baseLikelihoods(base);
		const auto pmd_likelihoods = pmd.baseLikelihoods(base, sem_likelihoods);

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

	// initialize RG
	BAM::TReadGroups rgs;
	rgs.add("testRG");

	SequencingError::TModels sem;
	sem.initializeNoRecal(rgs);

	BAM::TSequencedBase base;
	base.originalQuality_phredInt = 20;
	base.base                     = Base::A;
	base.readGroupID              = 0;
	base.distFrom5Prime           = 2;
	base.setReverseStrand(false);


	// initialize PMD
	coretools::TLog logfile;
	std::vector<size_t> ReadGroupsWithoutPMD;
	PMD::TModels pmd("doubleStrand:Empiric[0.3]:Empiric[0.1]", rgs, ReadGroupsWithoutPMD);

	for (size_t dfrom3 = 0; dfrom3 < 3; dfrom3 += 2) {
		base.distFrom3Prime = dfrom3;
		const auto from3    = base.distFrom3Prime < base.distFrom5Prime;
		for (Base b = Base::min; b < Base::max; ++b) {
			base.base = b;

			const auto sem_likelihoods = sem.baseLikelihoods(base);
			const auto pmd_likelihoods = pmd.baseLikelihoods(base, sem_likelihoods);

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
					EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], err / 3);
					if (from3 && trueBase == Base::G) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
										(1.0 - 0.1) * sem_likelihoods[Base::G] + 0.1 * sem_likelihoods[Base::A]);
					} else if (!from3 && trueBase == Base::C) {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
										(1.0 - 0.3) * sem_likelihoods[Base::C] + 0.3 * sem_likelihoods[Base::T]);
					} else {
						EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase], err / 3);
					}
				}
			}
		}
	}
}
