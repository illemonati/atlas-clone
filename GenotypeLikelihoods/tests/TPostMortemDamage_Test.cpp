//
// Created by vivian on 18.08.20.
//

#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TestCase.h"
#include <cstdint>

using namespace GenotypeLikelihoods;
using genometools::Base;

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
