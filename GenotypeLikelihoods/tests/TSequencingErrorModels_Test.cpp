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
    float errorOneThird = 0.01 / 3.0;
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    base.base = A;

    _sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoods);

    for(int first = 0; first < bases.size(); ++first){
        for(int second = 0; second < bases.size(); ++second) {
            //genotypes containing A
            if(first == 0 || second == 0){
                if(first == second){
                    //std::cout << "equal: first " << first << " second " << second << " oneMinusError " << oneMinusError << std::endl;
                    EXPECT_FLOAT_EQ(_baseLikelihoods[genoMap.toGenotype(bases[first], bases[second])], oneMinusError);
                } else {
                    EXPECT_FLOAT_EQ(_baseLikelihoods[genoMap.toGenotype(bases[first], bases[second])], 0.5 * (errorOneThird + oneMinusError));
                    std::cout << "NOT equal: first " << first << " second " << second << " error " << 0.5 * (errorOneThird + oneMinusError) << std::endl;
                }
            }

            //genotypes not containing A
            else {
                EXPECT_FLOAT_EQ(_baseLikelihoods[genoMap.toGenotype(bases[first], bases[second])], errorOneThird);
            }
        }
    }
}