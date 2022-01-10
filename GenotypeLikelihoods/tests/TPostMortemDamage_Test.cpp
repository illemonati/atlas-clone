//
// Created by vivian on 18.08.20.
//

#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TestCase.h"


using namespace GenotypeLikelihoods;
using genometools::Base;

TEST(TPostMortemDamage_test, baseANoPMD){
    TSequencingErrorModels sem;
    TBaseLikelihoods sem_likelihoods;

    TPostMortemDamage pmd;
    TBaseLikelihoods pmd_likelihoods;
    BAM::TSequencedBase base;
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    base.base = genometools::A;

    sem.fillBaseLikelihoods(base, sem_likelihoods);
    pmd.fillBaseLikelihoods(base, sem_likelihoods, pmd_likelihoods);

    for(Base b = Base::min(); b < Base::max(); ++b){
        base.base = (b);
        sem.fillBaseLikelihoods(base, sem_likelihoods);

        for(Base trueBase = Base::min(); trueBase < Base::max(); ++trueBase){
            //true base is A
            if(trueBase == b){
                EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], oneMinusError);
            } else {
                EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], errorOneThird);
            }
        }
    }
}

TEST(TPostMortemDamage_test, baseAWithPMD){
    TSequencingErrorModels sem;
    TBaseLikelihoods sem_likelihoods;
    TPostMortemDamage pmd;
    TBaseLikelihoods pmd_likelihoods;
    BAM::TSequencedBase base;

    //initialize RG
    BAM::TReadGroups ReadGroups;
    ReadGroups.add("testRG");

    //initialize PMD
    TLog logfile;
    std::vector<uint16_t> ReadGroupsWithoutPMD;
    pmd.initialize("doubleStrand:Empiric[0.5]:Empiric[0.5]", ReadGroups, &logfile, ReadGroupsWithoutPMD);

    //initialize base
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    base.base = genometools::A;

    for(Base b = Base::min(); b < Base::max(); ++b){
        //calculate base likelihoods
	    std::cout << "b " << b << '\n';
        base.base = b;
        base.readGroupID = 0;
        base.setReverseStrand(false);

        sem.fillBaseLikelihoods(base, sem_likelihoods);
        pmd.fillBaseLikelihoods(base, sem_likelihoods, pmd_likelihoods);

		for (Base trueBase = Base::min(); trueBase < Base::max(); ++trueBase) {
			std::cout <<"trueBase " <<  trueBase << '\n';
			std::cout << pmd_likelihoods << '\n';
			std::cout << sem_likelihoods << '\n';
			// true base is A
			if (trueBase == b) {
				if (trueBase == genometools::A || trueBase == genometools::T) {
					EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], oneMinusError);
				}
				else if (trueBase == genometools::C) {
					EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
							(1.0 - 0.5) * 0.5 + 0.5 * sem_likelihoods[genometools::T]);
				} else if (trueBase == genometools::G) {
					EXPECT_FLOAT_EQ(pmd_likelihoods[trueBase],
							(1.0 - 0.5) * 0.5 + 0.5 * sem_likelihoods[genometools::A]);
				}
			} else {
				EXPECT_FLOAT_EQ(sem_likelihoods[trueBase], errorOneThird);
			}
			std::cout << '\n';
		}
	}
}
