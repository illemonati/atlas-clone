//
// Created by vivian on 18.08.20.
//

#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TestCase.h"


using namespace GenotypeLikelihoods;
using genometools::Base;

TEST(TPostMortemDamage_test, baseANoPMD){
    TSequencingErrorModels _sequencingErrorModels;
    std::vector<TBaseLikelihoods> _baseLikelihoods;
    TBaseLikelihoods _baseLikelihoodsNoPMD;
    TPostMortemDamage _pmd;
    BAM::TSequencedBase base;
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    base.base = genometools::A;

    _baseLikelihoods.resize(1);
    _sequencingErrorModels.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD);
    _pmd.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

    for(Base b = Base::min(); b < Base::max(); ++b){
        base.base = (b);
        _sequencingErrorModels.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD);

        for(Base trueBase = Base::min(); trueBase < Base::max(); ++trueBase){
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
    std::vector<TBaseLikelihoods> _baseLikelihoods;
    TBaseLikelihoods _baseLikelihoodsNoPMD;
    TPostMortemDamage _pmd;
    BAM::TSequencedBase base;

    //initialize RG
    BAM::TReadGroups ReadGroups;
    ReadGroups.add("testRG");

    //initialize PMD
    TLog logfile;
    std::vector<uint16_t> ReadGroupsWithoutPMD;
    _pmd.initialize("Empiric[0.5]", ReadGroups, &logfile, ReadGroupsWithoutPMD);

    //initialize base
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    base.base = genometools::A;
    _baseLikelihoods.resize(1);

    for(Base b = Base::min(); b < Base::max(); ++b){
        //calculate base likelihoods
        base.base = b;
        base.readGroupID = 0;
        base.setReverseStrand(false);

        _sequencingErrorModels.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD);
        _pmd.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

        for(Base trueBase = Base::min(); trueBase < Base::max(); ++trueBase){
            //true base is A
            if(trueBase == b && (trueBase == genometools::A || trueBase == genometools::T)) {
                EXPECT_FLOAT_EQ(_baseLikelihoods[0][trueBase], oneMinusError);
            } else if(trueBase == b && (trueBase == genometools::C)){
                EXPECT_FLOAT_EQ(_baseLikelihoods[0][trueBase], (1.0 - 0.5) * 0.5 + 0.5 * _baseLikelihoodsNoPMD[genometools::T]);
            } else if(trueBase == b && (trueBase == genometools::G)){
                EXPECT_FLOAT_EQ(_baseLikelihoods[0][trueBase], (1.0 - 0.5) * 0.5 + 0.5 * _baseLikelihoodsNoPMD[genometools::A]);
//            } else {
//                EXPECT_FLOAT_EQ(_baseLikelihoods[0][trueBase], errorOneThird);
            }
        }
    }
}
