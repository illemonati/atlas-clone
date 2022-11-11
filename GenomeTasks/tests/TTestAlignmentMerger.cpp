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
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 10);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 100);
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
    EXPECT_EQ(firstRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 10);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 110);
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
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedLeft(), 0);
    EXPECT_EQ(secondRead.cigar().lengthSoftClippedRight(), 0);
    EXPECT_EQ(firstRead.position(), 10);  
    EXPECT_EQ(secondRead.position(), 10);
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
}