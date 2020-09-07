//
// Created by vivian on 18.08.20.
//

#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TGenotypeMap.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

using namespace GenotypeLikelihoods;

TEST(TPostMortemDamage_test, baseANoPMD){
    TSequencingErrorModels _sequencingErrorModels;
    std::vector<TBaseData> _baseLikelihoods;
    TBaseData _baseLikelihoodsNoPMD;
    TPostMortemDamage _pmd;
    BAM::TBase base;
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    base.base = A;
    TGenotypeMap genoMap;

    _baseLikelihoods.resize(1);
    _sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD);
    _pmd.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

    for(int b = 0; b < bases.size(); ++b){
        base.base = genoMap.getBaseOnlyCapitals(bases[b]);
        _sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD);

        for(int trueBase = 0; trueBase < bases.size(); ++trueBase){
            //true base is A
            if(trueBase == b){
                EXPECT_FLOAT_EQ(_baseLikelihoodsNoPMD[trueBase], oneMinusError);
            } else {
                EXPECT_FLOAT_EQ(_baseLikelihoodsNoPMD[trueBase], errorOneThird);
            }
        }
    }
}

TEST(TPostMortemDamage_test, baseAWithPMD){
    TSequencingErrorModels _sequencingErrorModels;
    std::vector<TBaseData> _baseLikelihoods;
    TBaseData _baseLikelihoodsNoPMD;
    TPostMortemDamage _pmd;
    BAM::TBase base;
    TGenotypeMap genoMap;

    //initialize RG
    BAM::TReadGroups ReadGroups;
    ReadGroups.add("testRG");

    //initialize PMD
    TParameters params;
    TLog logfile;
    params.addParameter("pmd", "Empiric[0.5]");
    _pmd.initialize(params, ReadGroups, &logfile);

    std::cout << "here!!" << std::endl;

    //initialize base
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    base.base = A;

//
//    //calculate likelihoods
//    _baseLikelihoods.resize(1);
//    _sequencingErrorModels.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD);
//    _pmd.calculateBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

//    for(int first = 0; first < bases.size(); ++first){
//        for(int second = 0; second < bases.size(); ++second){
//            if(first == second){
//                EXPECT_EQ(_baseLikelihoods[genoMap.toGenotype(bases[first], bases[second])], oneMinusError);
//            } else {
//                EXPECT_EQ(_baseLikelihoods[genoMap.toGenotype(bases[first], bases[second])], errorOneThird);
//            }
//        }
//    }
}
