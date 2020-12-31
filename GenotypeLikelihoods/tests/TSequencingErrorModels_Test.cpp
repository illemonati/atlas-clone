//
// Created by vivian on 17.08.20.
//

#include "TSequencingErrorModels.h"
#include "TGenotypeMap.h"
#include "TPostMortemDamage.h"
#include "TestCase.h"

using namespace GenotypeLikelihoods;

TEST(TSequencingErrorModels_test, oneBase){
    TSequencingErrorModels _sequencingErrorModels;
    TBaseData _baseLikelihoods;
    TGenotypeMap genoMap;

    BAM::TBase base;
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3.0;
    std::vector<char> bases = {'A', 'C', 'G', 'T'};

    for(int b = 0; b < bases.size(); ++b){
        base.base = genoMap.getBaseOnlyCapitals(bases[b]);
        _sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoods);

        for(int trueBase = 0; trueBase < bases.size(); ++trueBase){
            //true base is A
            if(trueBase == b){
                EXPECT_FLOAT_EQ(_baseLikelihoods[trueBase], oneMinusError);
            } else {
                EXPECT_FLOAT_EQ(_baseLikelihoods[trueBase], errorOneThird);
            }
        }
    }
}