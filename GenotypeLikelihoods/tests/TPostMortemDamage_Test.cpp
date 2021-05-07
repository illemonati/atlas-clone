//
// Created by vivian on 18.08.20.
//

#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"
#include "TGenotypeMap.h"
#include "TestCase.h"


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
    _sequencingErrorModels.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD);
    _pmd.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

    for(int b = 0; b < bases.size(); ++b){
        base.base = genoMap.getBaseOnlyCapitals(bases[b]);
        _sequencingErrorModels.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD);

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

    //initialize base
    base.originalQuality_phredInt = 20;
    float oneMinusError = 0.99;
    float errorOneThird = 0.01 / 3;
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    base.base = A;
    _baseLikelihoods.resize(1);

    for(int b = 0; b < bases.size(); ++b){
        //calculate base likelihoods
        base.base = genoMap.getBaseOnlyCapitals(bases[b]);
        base.readGroupID = 0;
        base.setReverseStrand(false);

        _sequencingErrorModels.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD);
        _pmd.fillBaseLikelihoods(base, _baseLikelihoodsNoPMD, _baseLikelihoods[0]);

        for(int trueBase = 0; trueBase < bases.size(); ++trueBase){
            //true base is A
            if(trueBase == b && (trueBase == A || trueBase == T)) {
                EXPECT_FLOAT_EQ(_baseLikelihoods[0][trueBase], oneMinusError);
            } else if(trueBase == b && (trueBase == C)){
                EXPECT_FLOAT_EQ(_baseLikelihoods[0][trueBase], (1.0 - 0.5) * 0.5 + 0.5 * _baseLikelihoodsNoPMD[T]);
            } else if(trueBase == b && (trueBase == G)){
                EXPECT_FLOAT_EQ(_baseLikelihoods[0][trueBase], (1.0 - 0.5) * 0.5 + 0.5 * _baseLikelihoodsNoPMD[A]);
//            } else {
//                EXPECT_FLOAT_EQ(_baseLikelihoods[0][trueBase], errorOneThird);
            }
        }
    }
}
