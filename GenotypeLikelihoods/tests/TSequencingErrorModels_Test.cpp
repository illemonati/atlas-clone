//
// Created by vivian on 17.08.20.
//

#include "TSequencingErrorModels.h"
#include "TGenotypeMap.h"
#include "TPostMortemDamage.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using namespace GenotypeLikelihoods;

TEST(TSequencingErrorModels_test, baseA){
    //TPostMortemDamage _pmd;
    TSequencingErrorModels _sequencingErrorModels;
    TBaseData _baseLikelihoods;
    TGenotypeMap genoMap;

    BAM::TBase base;
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    base.base = A;

    _sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoods);

    for(int first = 0; first < bases.size(); ++first){
        for(int second = 0; second < bases.size(); ++second){
            if(first == second){
                EXPECT_EQ(_baseLikelihoods[genoMap.toGenotype(bases[first], bases[second])], oneMinusError);
            } else {
                EXPECT_EQ(_baseLikelihoods[genoMap.toGenotype(bases[first], bases[second])], errorOneThird);
            }
        }
    }
}