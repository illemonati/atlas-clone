/*
 * TTestAlignmentMerger.cpp
 *
 *  Created on: Nov 11, 2022
 *      Author: Raphael
 */

#include "TAlignmentMerger.h"
#include "TAlignment.h"
#include "TCigar.h"
#include "gtest/gtest.h"

#include <algorithm>

using namespace GenomeTasks;

using namespace AlignmentMerger;

using namespace BAM;

//TODO:add unit test for highestQuality, check what happens when quality changes at border of overlap, anywhere in the overlap, first base etc.
//also check for all cases

//TODO: also test with soft clips and inserts etc. (especially the middle merge) 
//and put high quality inserts in and merge with highestquality --> should then merge the other read

TEST(TAlignmentMergerTest, forwardFirst_reverseSecond_mergeFirst){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 100);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar);

    TAlignmentMerger_firstMate mergeFirstMate;
    mergeFirstMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 10);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 100);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 90);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 100);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 100);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 4);
    cigar2.add('M', 3);
    cigar2.add('D', 5);
    cigar2.add('M', 32);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 40);
    cigar2.add('I', 10);
    cigar2.add('D', 2);
    cigar2.add('M', 3);
    cigar2.add('S', 6);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar2);

    mergeFirstMate.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 24);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 100);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 90);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 100);
}

TEST(TAlignmentMergerTest, forwardFirst_reverseSecond_mergeSecond){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 100);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar);

    TAlignmentMerger_secondMate mergeSecondMate;
    mergeSecondMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 10);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 110);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 90);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 100);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 4);
    cigar2.add('M', 3);
    cigar2.add('D', 5);
    cigar2.add('M', 32);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 40);
    cigar2.add('I', 10);
    cigar2.add('D', 2);
    cigar2.add('M', 3);
    cigar2.add('S', 6);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar2);

    mergeSecondMate.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 17);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 110);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 90);
}

TEST(TAlignmentMergerTest, forwardFirst_reverseSecond_mergeMiddle){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 100);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 5);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 5);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 105);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 95);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 95);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 100);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 4);
    cigar2.add('M', 3);
    cigar2.add('D', 5);
    cigar2.add('M', 32);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 40);
    cigar2.add('I', 10);
    cigar2.add('D', 2);
    cigar2.add('M', 3);
    cigar2.add('S', 6);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar2);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 19);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 17);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 110);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 95);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 90);
}

TEST(TAlignmentMergerTest, forwardFirst_reverseSecond_mergeMiddleOddOverlapLength){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 99);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(100);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(100);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(100);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar, vect, lowerQuality);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 5);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 6);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 105);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 95);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 94);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 99);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    lowerQuality[5] = genometools::PhredIntProbability(coretools::Probability (0.9));

    firstRead2.setSequenceQualities(cigar, vect, higherQuality);
    secondRead2.setSequenceQualities(cigar, vect, lowerQuality);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 5);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 104);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 94);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 95);


    TAlignment firstRead3(1, 10);
    TAlignment secondRead3(1, 99);

    firstRead3.setIsReverseStrand(false);
    secondRead3.setIsReverseStrand(true);

    TCigar newCigar;
    newCigar.add('S', 8);
    newCigar.add('M', 2);
    newCigar.add('I', 1);
    newCigar.add('M', 3);
    newCigar.add('D', 2);
    newCigar.add('M', 35);
    newCigar.add('I', 5);
    newCigar.add('M', 8);
    newCigar.add('D', 5);
    newCigar.add('M', 36);
    newCigar.add('I', 1);
    newCigar.add('M', 3);
    newCigar.add('D', 2);
    newCigar.add('M', 4);
    newCigar.add('S', 6);
 
    std::vector<genometools::PhredIntProbability> newHigherQuality;
    newHigherQuality.resize(112);
    std::fill(newHigherQuality.begin(),newHigherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> newLowerQuality;
    newLowerQuality.resize(112);
    std::fill(newLowerQuality.begin(),newLowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));
    vect.resize(112);

    firstRead3.setSequenceQualities(newCigar, vect, newHigherQuality);
    secondRead3.setSequenceQualities(newCigar, vect, newLowerQuality);

    mergeMiddle.merge(firstRead3, secondRead3);
    EXPECT_EQ(firstRead3.cigar().lengthSoftClippedRight(), 10);
    EXPECT_EQ(firstRead3.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(secondRead3.cigar().lengthSoftClippedLeft(), 14);
    EXPECT_EQ(secondRead3.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead3.position(), 10);  
    EXPECT_EQ(secondRead3.position(), 106);
    EXPECT_EQ(firstRead3.cigar().lengthMapped(), 94);
    EXPECT_EQ(secondRead3.cigar().lengthMapped(), 93);
}

TEST(TAlignmentMergerTest, forwardFirst_reverseSecond_mergeHighestQuality){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 100);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(100);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(100);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(100);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar, vect, lowerQuality);

    TAlignmentMerger_highestQuality mergeQual;
    mergeQual.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 10);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 110);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 90);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 100);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    firstRead2.setSequenceQualities(cigar, vect, lowerQuality);
    secondRead2.setSequenceQualities(cigar, vect, higherQuality);

    mergeQual.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 10);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 100);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 90);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 100);

    TAlignment firstRead3(1, 10);
    TAlignment secondRead3(1, 100);

    firstRead3.setIsReverseStrand(false);
    secondRead3.setIsReverseStrand(true);

    higherQuality[5] = genometools::PhredIntProbability(coretools::Probability (0.01));

    firstRead3.setSequenceQualities(cigar, vect, lowerQuality);
    secondRead3.setSequenceQualities(cigar, vect, higherQuality);

    mergeQual.merge(firstRead3, secondRead3);
    EXPECT_EQ(firstRead3.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(secondRead3.cigar().lengthSoftClippedLeft(), 10);
    EXPECT_EQ(firstRead3.position(), 10);  
    EXPECT_EQ(secondRead3.position(), 110);
    EXPECT_EQ(firstRead3.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead3.cigar().lengthMapped(), 90);


    TAlignment firstRead4(1, 10);
    TAlignment secondRead4(1, 100);

    firstRead4.setIsReverseStrand(false);
    secondRead4.setIsReverseStrand(true);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    std::vector<genometools::Base> vect2;
    vect2.push_back(genometools::char2base('G'));
    vect2.resize(112);

    std::vector<genometools::PhredIntProbability> higherQuality2;
    higherQuality2.resize(112);
    std::fill(higherQuality2.begin(),higherQuality2.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality2;
    lowerQuality2.resize(112);
    std::fill(lowerQuality2.begin(),lowerQuality2.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));

    lowerQuality2[12] = genometools::PhredIntProbability(coretools::Probability (0.9));
    firstRead4.setSequenceQualities(cigar2, vect2, higherQuality2);
    secondRead4.setSequenceQualities(cigar2, vect2, lowerQuality2);

    mergeQual.merge(firstRead4, secondRead4);
    EXPECT_EQ(firstRead4.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead4.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(secondRead4.cigar().lengthSoftClippedLeft(), 17);
    EXPECT_EQ(secondRead4.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead4.position(), 10);  
    EXPECT_EQ(secondRead4.position(), 110);
    EXPECT_EQ(firstRead4.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead4.cigar().lengthMapped(), 90);
}

TEST(TAlignmentMergerTest, sameStartAndEndPos_mergeFirst){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar);

    TAlignmentMerger_firstMate mergeFirstMate;
    mergeFirstMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 100);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 10);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 0);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 100);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 10);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar2);

    mergeFirstMate.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 112);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 10);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 0);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 100);
}

TEST(TAlignmentMergerTest, sameStartAndEndPos_mergeSecond){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar);

    TAlignmentMerger_secondMate mergeSecondMate;
    mergeSecondMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 100);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 110);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 0);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 10);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar2);

    mergeSecondMate.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 112);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 110);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 0);
}

TEST(TAlignmentMergerTest, sameStartAndEndPos_mergeMiddle){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 50);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 50);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 60);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 50);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 50);


    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 10);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar2);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 50);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 62);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 65);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 50);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 45);
}

TEST(TAlignmentMergerTest, sameStartAndEndPos_mergeMiddleOddOverlapLength){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 101);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(101);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(101);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(101);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar, vect, lowerQuality);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 50);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 51);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 61);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 51);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 50);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 10);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    lowerQuality[50] = genometools::PhredIntProbability(coretools::Probability (0.9));

    firstRead2.setSequenceQualities(cigar, vect, higherQuality);
    secondRead2.setSequenceQualities(cigar, vect, lowerQuality);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 51);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 50);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 60);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 50);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 51);


    TAlignment firstRead3(1, 10);
    TAlignment secondRead3(1, 10);

    firstRead3.setIsReverseStrand(false);
    secondRead3.setIsReverseStrand(true);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 5);
    cigar2.add('S', 6);

    std::vector<genometools::Base> vect2;
    vect2.push_back(genometools::char2base('G'));
    vect2.resize(113);

    std::vector<genometools::PhredIntProbability> higherQuality2;
    higherQuality2.resize(113);
    std::fill(higherQuality2.begin(),higherQuality2.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality2;
    lowerQuality2.resize(113);
    std::fill(lowerQuality2.begin(),lowerQuality2.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));

    firstRead3.setSequenceQualities(cigar2, vect2, higherQuality2);
    secondRead3.setSequenceQualities(cigar2, vect2, lowerQuality2);

    mergeMiddle.merge(firstRead3, secondRead3);
    EXPECT_EQ(firstRead3.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead3.cigar().lengthSoftClippedRight(), 51);
    EXPECT_EQ(secondRead3.cigar().lengthSoftClippedLeft(), 62);
    EXPECT_EQ(secondRead3.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead3.position(), 10);  
    EXPECT_EQ(secondRead3.position(), 65);
    EXPECT_EQ(firstRead3.cigar().lengthMapped(), 50);
    EXPECT_EQ(secondRead3.cigar().lengthMapped(), 46);
}

TEST(TAlignmentMergerTest, forwardStrandStartsFirstEndsLast_mergeFirst){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 80);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_firstMate mergeFirstMate;
    mergeFirstMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 90);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 20);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 10);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 80);


    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TCigar cigar3;
    cigar3.add('S', 8);
    cigar3.add('M', 2);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 20);
    cigar3.add('I', 5);
    cigar3.add('M', 8);
    cigar3.add('D', 5);
    cigar3.add('M', 20);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 4);
    cigar3.add('S', 6);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 20);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar3);

    mergeFirstMate.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 95);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 20);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 10);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 69);
}

TEST(TAlignmentMergerTest, forwardStrandStartsFirstEndsLast_mergeSecond){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 80);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_secondMate mergeSecondMate;
    mergeSecondMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 80);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 100);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 0);


    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TCigar cigar3;
    cigar3.add('S', 8);
    cigar3.add('M', 2);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 20);
    cigar3.add('I', 5);
    cigar3.add('M', 8);
    cigar3.add('D', 5);
    cigar3.add('M', 20);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 4);
    cigar3.add('S', 6);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 20);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar3);

    mergeSecondMate.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 81);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 89);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 0);
}

TEST(TAlignmentMergerTest, forwardStrandStartsFirstEndsLast_mergeMiddle){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 80);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 50);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 40);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 60);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 50);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 40);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TCigar cigar3;
    cigar3.add('S', 8);
    cigar3.add('M', 2);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 20);
    cigar3.add('I', 5);
    cigar3.add('M', 8);
    cigar3.add('D', 5);
    cigar3.add('M', 21);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 4);
    cigar3.add('S', 6);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 20);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar3);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 55);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 47);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 60);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 45);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 30);
}

TEST(TAlignmentMergerTest, forwardStrandStartsFirstEndsLast_mergeMiddleOddOverlapLength){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 19);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 81);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(100);

    std::vector<genometools::Base> vect2;
    vect2.push_back(genometools::char2base('G'));
    vect2.resize(81);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(100);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(81);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar1, vect2, lowerQuality);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 50);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 41);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 60);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 50);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 40);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 19);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    lowerQuality[40] = genometools::PhredIntProbability(coretools::Probability (0.9));

    firstRead2.setSequenceQualities(cigar, vect, higherQuality);
    secondRead2.setSequenceQualities(cigar1, vect2, lowerQuality);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 51);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 40);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 59);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 49);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 41);


    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TCigar cigar3;
    cigar3.add('S', 8);
    cigar3.add('M', 2);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 20);
    cigar3.add('I', 5);
    cigar3.add('M', 8);
    cigar3.add('D', 5);
    cigar3.add('M', 20);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 4);
    cigar3.add('S', 6);

    TAlignment firstRead3(1, 10);
    TAlignment secondRead3(1, 20);

    firstRead3.setIsReverseStrand(false);
    secondRead3.setIsReverseStrand(true);

    std::vector<genometools::Base> vect3;
    vect3.push_back(genometools::char2base('G'));
    vect3.resize(112);

    std::vector<genometools::Base> vect4;
    vect4.push_back(genometools::char2base('G'));
    vect4.resize(81);

    std::vector<genometools::PhredIntProbability> higherQuality2;
    higherQuality2.resize(81);
    std::fill(higherQuality2.begin(),higherQuality2.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality2;
    lowerQuality2.resize(112);
    std::fill(lowerQuality2.begin(),lowerQuality2.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));

    firstRead3.setSequenceQualities(cigar2, vect3, lowerQuality2);
    secondRead3.setSequenceQualities(cigar3, vect4, higherQuality2);

    mergeMiddle.merge(firstRead3, secondRead3);
    EXPECT_EQ(firstRead3.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead3.cigar().lengthSoftClippedRight(), 56);
    EXPECT_EQ(secondRead3.cigar().lengthSoftClippedLeft(), 46);
    EXPECT_EQ(secondRead3.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead3.position(), 10);  
    EXPECT_EQ(secondRead3.position(), 54);
    EXPECT_EQ(firstRead3.cigar().lengthMapped(), 44);
    EXPECT_EQ(secondRead3.cigar().lengthMapped(), 35);
}


TEST(TAlignmentMergerTest, reverseStrandStartsFirstEndsLast_mergeFirst){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 80);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_firstMate mergeFirstMate;
    mergeFirstMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 90);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 100);  
    EXPECT_EQ(secondRead.position(), 20);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 10);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 80);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TCigar cigar3;
    cigar3.add('S', 8);
    cigar3.add('M', 2);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 20);
    cigar3.add('I', 5);
    cigar3.add('M', 8);
    cigar3.add('D', 5);
    cigar3.add('M', 21);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 4);
    cigar3.add('S', 6);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 20);

    firstRead2.setIsReverseStrand(true);
    secondRead2.setIsReverseStrand(false);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar3);

    mergeFirstMate.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 87);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(firstRead2.position(), 90);  
    EXPECT_EQ(secondRead2.position(), 20);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 20);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 70);
}

TEST(TAlignmentMergerTest, reverseStrandStartsFirstEndsLast_mergeSecond){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 80);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_secondMate mergeSecondMate;
    mergeSecondMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 80);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 20);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 0);

    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TCigar cigar3;
    cigar3.add('S', 8);
    cigar3.add('M', 2);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 20);
    cigar3.add('I', 5);
    cigar3.add('M', 8);
    cigar3.add('D', 5);
    cigar3.add('M', 21);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 4);
    cigar3.add('S', 6);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 20);

    firstRead2.setIsReverseStrand(true);
    secondRead2.setIsReverseStrand(false);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar3);

    mergeSecondMate.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 82);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 20);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 100);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 0);
}

TEST(TAlignmentMergerTest, reverseStrandStartsFirstEndsLast_mergeMiddle){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 80);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 50);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 40);
    EXPECT_EQ(firstRead.position(), 60);  
    EXPECT_EQ(secondRead.position(), 20);
    EXPECT_EQ(firstRead.cigar().lengthMapped(), 50);
    EXPECT_EQ(secondRead.cigar().lengthMapped(), 40);


    TCigar cigar2;
    cigar2.add('S', 8);
    cigar2.add('M', 2);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 35);
    cigar2.add('I', 5);
    cigar2.add('M', 8);//-3
    cigar2.add('D', 5);
    cigar2.add('M', 36);
    cigar2.add('I', 1);
    cigar2.add('M', 3);
    cigar2.add('D', 2);
    cigar2.add('M', 4);
    cigar2.add('S', 6);

    TCigar cigar3;
    cigar3.add('S', 8);
    cigar3.add('M', 2);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 20);
    cigar3.add('I', 5);
    cigar3.add('M', 8);
    cigar3.add('D', 5);
    cigar3.add('M', 21);
    cigar3.add('I', 1);
    cigar3.add('M', 3);
    cigar3.add('D', 2);
    cigar3.add('M', 4);
    cigar3.add('S', 6);

    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 20);

    firstRead2.setIsReverseStrand(true);
    secondRead2.setIsReverseStrand(false);

    firstRead2.setCigarForUnitTest(cigar2);
    secondRead2.setCigarForUnitTest(cigar3);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 57);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 6);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 8);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 35);
    EXPECT_EQ(firstRead2.position(), 55);  
    EXPECT_EQ(secondRead2.position(), 20);
    EXPECT_EQ(firstRead2.cigar().lengthMapped(), 55);
    EXPECT_EQ(secondRead2.cigar().lengthMapped(), 35);
}

TEST(TAlignmentMergerTest, reverseStrandStartsFirstEndsLast_mergeMiddleOddOverlapLength){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 19);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 81);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(100);

    std::vector<genometools::Base> vect2;
    vect2.push_back(genometools::char2base('G'));
    vect2.resize(81);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(100);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(81);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar1, vect2, lowerQuality);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 49);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 41);
    EXPECT_EQ(firstRead.position(), 59);  
    EXPECT_EQ(secondRead.position(), 19);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 19);

    firstRead2.setIsReverseStrand(true);
    secondRead2.setIsReverseStrand(false);

    lowerQuality[40] = genometools::PhredIntProbability(coretools::Probability (0.9));

    firstRead2.setSequenceQualities(cigar, vect, higherQuality);
    secondRead2.setSequenceQualities(cigar1, vect2, lowerQuality);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 50);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 40);
    EXPECT_EQ(firstRead2.position(), 60);  
    EXPECT_EQ(secondRead2.position(), 19);
}


TEST(TAlignmentMergerTest, sameStartForwardEndsLast_mergeFirst){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 80);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_firstMate mergeFirstMate;
    mergeFirstMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 100);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 10);
}

TEST(TAlignmentMergerTest, sameStartForwardEndsLast_mergeSecond){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 80);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_secondMate mergeSecondMate;
    mergeSecondMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 80);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 90);
}

TEST(TAlignmentMergerTest, sameStartForwardEndsLast_mergeMiddle){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 80);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 60);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 40);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 50);
}

TEST(TAlignmentMergerTest, sameStartForwardEndsLast_mergeMiddleOddOverlapLength){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 81);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(100);

    std::vector<genometools::Base> vect2;
    vect2.push_back(genometools::char2base('G'));
    vect2.resize(81);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(100);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(81);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar1, vect2, lowerQuality);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 59);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 41);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 51);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 10);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    lowerQuality[40] = genometools::PhredIntProbability(coretools::Probability (0.9));

    firstRead2.setSequenceQualities(cigar, vect, higherQuality);
    secondRead2.setSequenceQualities(cigar1, vect2, lowerQuality);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 60);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 40);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 50);
}


TEST(TAlignmentMergerTest, sameStartReverseEndsLast_mergeFirst){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 80);
    TCigar cigar1;
    cigar1.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_firstMate mergeFirstMate;
    mergeFirstMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 80);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 10);
}

TEST(TAlignmentMergerTest, sameStartReverseEndsLast_mergeSecond){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 80);
    TCigar cigar1;
    cigar1.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_secondMate mergeSecondMate;
    mergeSecondMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 80);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 90);
}

TEST(TAlignmentMergerTest, sameStartReverseEndsLast_mergeMiddle){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 80);
    TCigar cigar1;
    cigar1.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 40);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 40);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 50);
}

TEST(TAlignmentMergerTest, sameStartReverseEndsLast_mergeMiddleOddOverlapLength){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 10);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 81);
    TCigar cigar1;
    cigar1.add('M', 100);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(81);

    std::vector<genometools::Base> vect2;
    vect2.push_back(genometools::char2base('G'));
    vect2.resize(100);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(81);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(100);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar1, vect2, lowerQuality);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 40);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 41);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 51);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 10);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    lowerQuality[40] = genometools::PhredIntProbability(coretools::Probability (0.9));

    firstRead2.setSequenceQualities(cigar, vect, higherQuality);
    secondRead2.setSequenceQualities(cigar1, vect2, lowerQuality);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 41);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 40);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 50);
}


TEST(TAlignmentMergerTest, forwardStartsFirstSameEnd_mergeFirst){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 90);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_firstMate mergeFirstMate;
    mergeFirstMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 90);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 20);
}

TEST(TAlignmentMergerTest, forwardStartsFirstSameEnd_mergeSecond){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 90);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_secondMate mergeSecondMate;
    mergeSecondMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 90);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 110);
}

TEST(TAlignmentMergerTest, forwardStartsFirstSameEnd_mergeMiddle){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 90);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 45);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 45);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 65);
}

TEST(TAlignmentMergerTest, forwardStartsFirstSameEnd_mergeMiddleOddOverlapLength){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 21);

    firstRead.setIsReverseStrand(false);
    secondRead.setIsReverseStrand(true);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 89);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(100);

    std::vector<genometools::Base> vect2;
    vect2.push_back(genometools::char2base('G'));
    vect2.resize(89);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(100);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(89);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar1, vect2, lowerQuality);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 44);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 45);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 66);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 21);

    firstRead2.setIsReverseStrand(false);
    secondRead2.setIsReverseStrand(true);

    lowerQuality[44] = genometools::PhredIntProbability(coretools::Probability (0.9));

    firstRead2.setSequenceQualities(cigar, vect, higherQuality);
    secondRead2.setSequenceQualities(cigar1, vect2, lowerQuality);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedRight(), 45);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedLeft(), 44);
    EXPECT_EQ(firstRead2.position(), 10);  
    EXPECT_EQ(secondRead2.position(), 65);
}


TEST(TAlignmentMergerTest,reverseStartsFirstSameEnd_mergeFirst){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 90);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_firstMate mergeFirstMate;
    mergeFirstMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 100);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 110);  
    EXPECT_EQ(secondRead.position(), 20);
}

TEST(TAlignmentMergerTest, reverseStartsFirstSameEnd_mergeSecond){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 90);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_secondMate mergeSecondMate;
    mergeSecondMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 90);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 20);
}

TEST(TAlignmentMergerTest, reverseStartsFirstSameEnd_mergeMiddle){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 90);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar1);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 55);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 45);
    EXPECT_EQ(firstRead.position(), 65);  
    EXPECT_EQ(secondRead.position(), 20);
}

TEST(TAlignmentMergerTest, reverseStartsFirstSameEnd_mergeMiddleOddOverlapLength){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 21);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);
    TCigar cigar1;
    cigar1.add('M', 89);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(100);

    std::vector<genometools::Base> vect2;
    vect2.push_back(genometools::char2base('G'));
    vect2.resize(89);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(100);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(89);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar1, vect2, lowerQuality);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 55);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 45);
    EXPECT_EQ(firstRead.position(), 65);  
    EXPECT_EQ(secondRead.position(), 21);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 21);

    firstRead2.setIsReverseStrand(true);
    secondRead2.setIsReverseStrand(false);

    lowerQuality[44] = genometools::PhredIntProbability(coretools::Probability (0.9));

    firstRead2.setSequenceQualities(cigar, vect, higherQuality);
    secondRead2.setSequenceQualities(cigar1, vect2, lowerQuality);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 56);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 44);
    EXPECT_EQ(firstRead2.position(), 66);  
    EXPECT_EQ(secondRead2.position(), 21);
}


TEST(TAlignmentMergerTest, reverseFirst_forwardSecond_mergeFirst){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar);

    TAlignmentMerger_firstMate mergeFirstMate;
    mergeFirstMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 100);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 110);  
    EXPECT_EQ(secondRead.position(), 20);
}

TEST(TAlignmentMergerTest, reverseFirst_forwardSecond_mergeSecond){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar);

    TAlignmentMerger_secondMate mergeSecondMate;
    mergeSecondMate.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 100);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 20);
}

TEST(TAlignmentMergerTest, reverseFirst_forwardSecond_mergeMiddle){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 20);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);

    firstRead.setCigarForUnitTest(cigar);
    secondRead.setCigarForUnitTest(cigar);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 55);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 55);
    EXPECT_EQ(firstRead.position(), 65);  
    EXPECT_EQ(secondRead.position(), 20);
}

TEST(TAlignmentMergerTest, reverseFirst_forwardSecond_mergeMiddleOddOverlapLength){
    TAlignment firstRead(1, 10);
    TAlignment secondRead(1, 21);

    firstRead.setIsReverseStrand(true);
    secondRead.setIsReverseStrand(false);

    TCigar cigar;
    cigar.add('M', 100);

    std::vector<genometools::Base> vect;
    vect.push_back(genometools::char2base('G'));
    vect.resize(100);

    std::vector<genometools::Base> vect2;
    vect2.push_back(genometools::char2base('G'));
    vect2.resize(100);

    std::vector<genometools::PhredIntProbability> higherQuality;
    higherQuality.resize(100);
    std::fill(higherQuality.begin(),higherQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.5)));

    std::vector<genometools::PhredIntProbability> lowerQuality;
    lowerQuality.resize(100);
    std::fill(lowerQuality.begin(),lowerQuality.end(),genometools::PhredIntProbability(coretools::Probability (0.1)));


    firstRead.setSequenceQualities(cigar, vect, higherQuality);
    secondRead.setSequenceQualities(cigar, vect2, lowerQuality);

    TAlignmentMerger_middle mergeMiddle;
    mergeMiddle.merge(firstRead, secondRead);
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedLeft(), 55);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 56);
    EXPECT_EQ(firstRead.position(), 65);  
    EXPECT_EQ(secondRead.position(), 21);


    TAlignment firstRead2(1, 10);
    TAlignment secondRead2(1, 21);

    firstRead2.setIsReverseStrand(true);
    secondRead2.setIsReverseStrand(false);

    lowerQuality[44] = genometools::PhredIntProbability(coretools::Probability (0.9));

    firstRead2.setSequenceQualities(cigar, vect, higherQuality);
    secondRead2.setSequenceQualities(cigar, vect2, lowerQuality);

    mergeMiddle.merge(firstRead2, secondRead2);
    EXPECT_EQ(firstRead2.cigar().lengthSoftClippedLeft(), 56);
    EXPECT_EQ(secondRead2.cigar().lengthSoftClippedRight(), 55);
    EXPECT_EQ(firstRead2.position(), 66);  
    EXPECT_EQ(secondRead2.position(), 21);
}
