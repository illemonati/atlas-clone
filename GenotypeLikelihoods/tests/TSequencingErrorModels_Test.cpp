//
// Created by vivian on 17.08.20.
//

#include "TSequencingErrorModels.h"
#include "TGenotypeMap.h"
#include "TPostMortemDamage.h"
#include "TestCase.h"

using namespace GenotypeLikelihoods;
using genometools::Base;

TEST(TSequencingErrorModels_test, oneBase){
    TSequencingErrorModels _sequencingErrorModels;
    TBaseLikelihoods _baseLikelihoods;

    BAM::TSequencedBase base;
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3.0;

    for(Base b = Base::min(); b < Base::max(); ++b){
        base.base = b;
        _sequencingErrorModels.fillBaseLikelihoods(base, _baseLikelihoods);

        for(Base trueBase = Base::min(); trueBase < Base::max(); ++trueBase){
            //true base is A
            if(trueBase == b){
                EXPECT_FLOAT_EQ(_baseLikelihoods[trueBase], oneMinusError);
            } else {
                EXPECT_FLOAT_EQ(_baseLikelihoods[trueBase], errorOneThird);
            }
        }
    }
}
