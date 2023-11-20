#include "TBamFile.h"
#include "gtest/gtest.h"

#include <stdio.h>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "TBamFilter.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "TCigar.h"
#include "coretools/Files/TFile.h"
#include "TGenome_OLD.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TReadGroups.h"
#include "TSamFlags.h"
#include "TSamHeader.h"
#include "TSequencedBase.h"
#include "TSite.h"
#include "TTestBamFile.h"
#include "TWindow.h"
#include "coretools/Math/counters.h"
#include "coretools/Strings/stringFunctions.h"

using coretools::instances::parameters;

using genometools::Base;
using genometools::PhredIntProbability;

//-------------------------------------------------------------
// TBamFile - simple writing and reading
//-------------------------------------------------------------

class TBamFile_Test_ReadWrite : public ::testing::Test {
protected:
    std::string _filename = "testBAM.bam";

public:
    std::unique_ptr<BAM::TTestBamFile> outputBam;
    std::unique_ptr<BAM::TBamFile> inputBam;

    void write(){
        //settings
        std::vector<size_t> chrLength = {1000, 2000, 3000};
        size_t numReadGroups = 2;

        //open BAM file for writing
        outputBam = std::make_unique<BAM::TTestBamFile>(_filename, chrLength, numReadGroups);

        //write alignments
        outputBam->writeDummyAlignments(100);
        outputBam->closeOutput();
    }

    void read(){
        //open BAM file for reading
        inputBam = std::make_unique<BAM::TBamFile>(_filename);
    }

    void SetUp() override{
        write();
        read();
    }

    void TearDown() override {};
};

TEST_F(TBamFile_Test_ReadWrite, chromosomes){
    // get chromosomes
    genometools::TChromosomes chromosomesWritten = outputBam->chromosomes();
    genometools::TChromosomes chromosomesRead = inputBam->chromosomes();

    // compare size
    EXPECT_EQ(chromosomesWritten.size(), chromosomesRead.size());

    // compare values of iterators
    auto itRead = chromosomesRead.begin();
    for (auto itWritten = chromosomesWritten.begin(); itWritten != chromosomesWritten.end(); itWritten++, itRead++){
        EXPECT_EQ(itWritten->name(), itRead->name());
        EXPECT_EQ(itWritten->length(), itRead->length());
        EXPECT_EQ(itWritten->refID(), itRead->refID());
        EXPECT_EQ(itWritten->start().position(), itRead->start().position());
        EXPECT_EQ(itWritten->end().position(), itRead->end().position());
    }
}

TEST_F(TBamFile_Test_ReadWrite, samHeader){
    // get header
    BAM::TSamHeader headerWritten = outputBam->header();
    BAM::TSamHeader headerRead = inputBam->samHeader();

    // compare values
    EXPECT_EQ(headerWritten.version(), headerRead.version());
    EXPECT_EQ(headerWritten.grouping(), headerRead.grouping());
    EXPECT_EQ(headerWritten.sortOrder(), headerRead.sortOrder());
    EXPECT_EQ(headerWritten.subSorting(), headerRead.subSorting());
}

TEST_F(TBamFile_Test_ReadWrite, readGroups){
    // get read groups
    BAM::TReadGroups readGroupsWritten = outputBam->readGroups();
    BAM::TReadGroups readGroupsRead = inputBam->readGroups();

    // compare size
    EXPECT_EQ(readGroupsWritten.size(), readGroupsRead.size());

    // compare values of iterators
    auto itRead = readGroupsRead.begin();
    for (auto itWritten = readGroupsWritten.begin(); itWritten != readGroupsWritten.end(); itWritten++, itRead++){
        EXPECT_EQ(itWritten->name_ID, itRead->name_ID);
        EXPECT_EQ(itWritten->sequencingCenter_CN, itRead->sequencingCenter_CN);
        EXPECT_EQ(itWritten->description_DS, itRead->description_DS);
        EXPECT_EQ(itWritten->sample_SM, itRead->sample_SM);
        EXPECT_EQ(itWritten->sequencingTechnology_PL, itRead->sequencingTechnology_PL);
    }
}

TEST_F(TBamFile_Test_ReadWrite, alignments){
    // read
    BAM::TAlignment alignmentRead;
    auto alignmentWritten = outputBam->beginWrittenAlignments();
    while (inputBam->readNextAlignment(alignmentRead)){
    	alignmentRead.parse();

        // basic attributes of TAlignment
        EXPECT_EQ(alignmentWritten->name(), alignmentRead.name());
        EXPECT_EQ(alignmentWritten->readGroupId(), alignmentRead.readGroupId());
        EXPECT_EQ(alignmentWritten->parsedLength(), alignmentRead.parsedLength());
        EXPECT_EQ(alignmentWritten->insertSize(), alignmentRead.insertSize());
        EXPECT_EQ(alignmentWritten->mappingQuality(), alignmentRead.mappingQuality());
        EXPECT_EQ(alignmentWritten->flags(), alignmentRead.flags());

        // position attributes of TAlignment
        EXPECT_EQ(alignmentWritten->lastAlingedInternalPos(), alignmentRead.lastAlingedInternalPos());
        EXPECT_EQ(alignmentWritten->lastAlignedPositionWithRespectToRef(), alignmentRead.lastAlignedPositionWithRespectToRef());
        for (size_t i = 0; i < alignmentWritten->parsedLength(); i++){
            EXPECT_EQ(alignmentWritten->isAlignedAtInternalPos(i), alignmentRead.isAlignedAtInternalPos(i));
			if (alignmentWritten->isAlignedAtInternalPos(i)) {
				EXPECT_EQ(alignmentWritten->positionInRef(i), alignmentRead.positionInRef(i));
			}
		}
        EXPECT_EQ(alignmentWritten->mateGenomicPosition(), alignmentRead.mateGenomicPosition());
        EXPECT_EQ(alignmentWritten->matePosition(), alignmentRead.matePosition());
        EXPECT_EQ(alignmentWritten->mateRefID(), alignmentRead.mateRefID());

        // attributes of TCigar
        EXPECT_EQ(alignmentWritten->cigar().lengthAligned(), alignmentRead.cigar().lengthAligned());
        EXPECT_EQ(alignmentWritten->cigar().lengthInserted(), alignmentRead.cigar().lengthInserted());
        EXPECT_EQ(alignmentWritten->cigar().lengthDeleted(), alignmentRead.cigar().lengthDeleted());
        EXPECT_EQ(alignmentWritten->cigar().lengthSoftClippedLeft(), alignmentRead.cigar().lengthSoftClippedLeft());
        EXPECT_EQ(alignmentWritten->cigar().lengthSoftClippedRight(), alignmentRead.cigar().lengthSoftClippedRight());
        EXPECT_EQ(alignmentWritten->cigar().lengthSoftClipped(), alignmentRead.cigar().lengthSoftClipped());
        EXPECT_EQ(alignmentWritten->cigar().lengthSequenced(), alignmentRead.cigar().lengthSequenced());
        EXPECT_EQ(alignmentWritten->cigar().lengthRead(), alignmentRead.cigar().lengthRead());

        // sequence attributes of TAlignment
        EXPECT_EQ(alignmentWritten->sequence(), alignmentRead.sequence());
        EXPECT_EQ(alignmentWritten->qualities(), alignmentRead.qualities());
        EXPECT_EQ(alignmentWritten->isReverseStrand(), alignmentRead.isReverseStrand());
        EXPECT_EQ(alignmentWritten->isPaired(), alignmentRead.isPaired());
        EXPECT_EQ(alignmentWritten->isProperPair(), alignmentRead.isProperPair());
        EXPECT_EQ(alignmentWritten->isParsed(), alignmentRead.isParsed());

        // iterate over bases of TAlignment
        auto baseRead = alignmentRead.begin();
        for (auto baseWritten = alignmentWritten->begin(); baseWritten != alignmentWritten->end(); baseWritten++, baseRead++){
            // all attributes of TBase
            EXPECT_EQ(baseWritten->originalQuality_phredInt, baseRead->originalQuality_phredInt);
            /*
            EXPECT_EQ(baseWritten->recalibratedQualityAsPhredInt, baseRead->recalibratedQualityAsPhredInt);
            EXPECT_EQ(baseWritten->distFrom3Prime, baseRead->distFrom3Prime);
            EXPECT_EQ(baseWritten->distFrom5Prime, baseRead->distFrom5Prime);
            EXPECT_EQ(baseWritten->readGroupID, baseRead->readGroupID);
            EXPECT_EQ(baseWritten->fragmentLength, baseRead->fragmentLength);
            EXPECT_EQ(baseWritten->mappingQuality, baseRead->mappingQuality);
            EXPECT_EQ(baseWritten->isReverseStrand(), baseRead->isReverseStrand());
            EXPECT_EQ(baseWritten->isAligned(), baseRead->isAligned());
            EXPECT_EQ(baseWritten->isSecondMate(), baseRead->isSecondMate());
            EXPECT_TRUE(baseWritten->base == baseRead->base);
            EXPECT_TRUE(baseWritten->context == baseRead->context);
            */
        }
        EXPECT_EQ(baseRead, alignmentRead.end());

        alignmentWritten++;
    }

    EXPECT_EQ(alignmentWritten, outputBam->endWrittenAlignments());
}


//-------------------------------------------------------------
// TBamFile - windows
//-------------------------------------------------------------

class TGenomeWindow_Test : public GenomeTasks::old::TGenome_windows {
protected:
    std::vector<genometools::TGenomeWindow> _windows_visited;

    void _handleWindow() override{
        // store sites
        std::vector<GenotypeLikelihoods::TSite> tmp;
        for (auto & it : _window)
            tmp.emplace_back(it);
        sites.emplace_back(tmp);
        // store genometools::TGenomeWindow
        _windows_visited.emplace_back(this->_window);
        // store GenotypeLikelihoods::TWindow attributes
        depth.emplace_back(_window.depth());
        fractionSitesNoData.emplace_back(_window.fractionSitesNoData());
        fractionDepthAtLeastTwo.emplace_back(_window.fractionDepthAtLeastTwo());
        numSitesWithData.emplace_back(_window.numSitesWithData());
        numReadsInWindow.emplace_back(_window.numReadsInWindow());
    };
	virtual void _handleAlignment() override {}

public:
    // storage
    std::vector<double> depth;
    std::vector<double> fractionSitesNoData;
    std::vector<double> fractionDepthAtLeastTwo;
    std::vector<size_t> numSitesWithData;
    std::vector<size_t> numReadsInWindow;
    std::vector<std::vector<GenotypeLikelihoods::TSite>> sites;

    TGenomeWindow_Test() : GenomeTasks::old::TGenome_windows() {};

    void traverse(){
        _traverseBAMWindows();
    }

    // size
    size_t numWindows(){ return _windows_visited.size(); };
    // loop
    std::vector<genometools::TGenomeWindow >::iterator begin(){ return _windows_visited.begin(); };
    std::vector<genometools::TGenomeWindow >::iterator end(){ return _windows_visited.end(); };
    // access
    genometools::TGenomeWindow& operator[](size_t pos){ return _windows_visited[pos]; };
};

class TBamFile_Test_Windows : public ::testing::Test {
public:
    std::unique_ptr<BAM::TTestBamFile> outputBam;
    std::unique_ptr<TGenomeWindow_Test> genomeWindow;
    std::string filename = "testBAM.bam";

    void write(){
        //settings
        std::vector<size_t> chrLength = {250, 50, 199, 80, 177};
        size_t numReadGroups = 2;

        //open BAM file for writing
        outputBam = std::make_unique<BAM::TTestBamFile>(filename, chrLength, numReadGroups);

        //write alignments

        // 1) overlapping alignments inside one window, cigar is all M
        BAM::TCigar onlyMCigar;
        onlyMCigar.add('M', 20);

        outputBam->writeDummyAlignment(Base::A, PhredIntProbability(1), genometools::TGenomePosition(0, 0), onlyMCigar);
        outputBam->writeDummyAlignment(Base::C, PhredIntProbability(2), genometools::TGenomePosition(0, 10), onlyMCigar);
        outputBam->writeDummyAlignment(Base::G, PhredIntProbability(3), genometools::TGenomePosition(0, 20), onlyMCigar);

        // 2) alignments overlap 2 windows
        outputBam->writeDummyAlignment(Base::T, PhredIntProbability(4), genometools::TGenomePosition(0, 80), onlyMCigar);
        outputBam->writeDummyAlignment(Base::A, PhredIntProbability(5), genometools::TGenomePosition(0, 90), onlyMCigar);
        outputBam->writeDummyAlignment(Base::C, PhredIntProbability(6), genometools::TGenomePosition(0, 95), onlyMCigar);
        outputBam->writeDummyAlignment(Base::G, PhredIntProbability(7), genometools::TGenomePosition(0, 100), onlyMCigar);

        // 3) one alignment inside 1 window
        outputBam->writeDummyAlignment(Base::T, PhredIntProbability(8), genometools::TGenomePosition(0, 220), onlyMCigar);

        // 4) only 1 window per chromosome
        outputBam->writeDummyAlignment(Base::A, PhredIntProbability(9), genometools::TGenomePosition(1, 10), onlyMCigar);

        // 5) empty window
        outputBam->writeDummyAlignment(Base::C, PhredIntProbability(10), genometools::TGenomePosition(2, 10), onlyMCigar);

        // 6) empty chromosome

        // 7) cigar strings more complicated (with soft-clips, insertions and deletions)
        BAM::TCigar mixedCigar;
        mixedCigar.add('S', 3);
        mixedCigar.add('M', 5);
        mixedCigar.add('D', 5);
        mixedCigar.add('I', 5);
        mixedCigar.add('S', 2);
        outputBam->writeDummyAlignment(Base::A, PhredIntProbability(1), genometools::TGenomePosition(4, 0), mixedCigar);
        outputBam->writeDummyAlignment(Base::C, PhredIntProbability(2), genometools::TGenomePosition(4, 4), mixedCigar);
        outputBam->writeDummyAlignment(Base::G, PhredIntProbability(3), genometools::TGenomePosition(4, 8), mixedCigar);

        // 8) last window is empty


        outputBam->closeOutput();
    }

    void read(){
        // first set window size in parameters
        parameters().add("window", "100");
        parameters().add("bam", filename);
        parameters().add("filterReadLength", "[0,20]");

        // create instance of TGenomeWindow
        genomeWindow = std::make_unique<TGenomeWindow_Test>();

        // now traverse bam file
        genomeWindow->traverse();
    }

    void SetUp() override{
        write();
        read();
    }

    void TearDown() override {};
};

TEST_F(TBamFile_Test_Windows, numWindows){
    EXPECT_EQ(genomeWindow->numWindows(),9);
}

TEST_F(TBamFile_Test_Windows, refIDWindows){
    // do windows store correct refID?

    // 1. chr (cases 1), 2) and 3))
    EXPECT_EQ((*genomeWindow)[0].refID(), 0);
    EXPECT_EQ((*genomeWindow)[1].refID(), 0);
    EXPECT_EQ((*genomeWindow)[2].refID(), 0);

    // 2. chr (case 4))
    EXPECT_EQ((*genomeWindow)[3].refID(), 1);

    // 3. chr (case 5))
    EXPECT_EQ((*genomeWindow)[4].refID(), 2);
    EXPECT_EQ((*genomeWindow)[5].refID(), 2);

    // 4. chr (case6))
    EXPECT_EQ((*genomeWindow)[6].refID(), 3);

    // 5. chr (case 7) and 8))
    EXPECT_EQ((*genomeWindow)[7].refID(), 4);
    EXPECT_EQ((*genomeWindow)[8].refID(), 4);
}

TEST_F(TBamFile_Test_Windows, positionsWindows){
    // do windows store correct positions? I.e., if chromosome is shorter than window, is window correctly adjusted?

    // 1. chr (cases 1), 2) and 3))
    EXPECT_EQ((*genomeWindow)[0].fromOnChr(), 0); EXPECT_EQ((*genomeWindow)[0].toOnChr(), 100);
    EXPECT_EQ((*genomeWindow)[1].fromOnChr(), 100); EXPECT_EQ((*genomeWindow)[1].toOnChr(), 200);
    EXPECT_EQ((*genomeWindow)[2].fromOnChr(), 200); EXPECT_EQ((*genomeWindow)[2].toOnChr(), 250);

    // 2. chr (case 4))
    EXPECT_EQ((*genomeWindow)[3].fromOnChr(), 0); EXPECT_EQ((*genomeWindow)[3].toOnChr(), 50);

    // 3. chr (case 5))
    EXPECT_EQ((*genomeWindow)[4].fromOnChr(), 0); EXPECT_EQ((*genomeWindow)[4].toOnChr(), 100);
    EXPECT_EQ((*genomeWindow)[5].fromOnChr(), 100); EXPECT_EQ((*genomeWindow)[5].toOnChr(), 199);

    // 4. chr (case6))
    EXPECT_EQ((*genomeWindow)[6].fromOnChr(), 0); EXPECT_EQ((*genomeWindow)[6].toOnChr(), 80);

    // 5. chr (case 7) and 8))
    EXPECT_EQ((*genomeWindow)[7].fromOnChr(), 0); EXPECT_EQ((*genomeWindow)[7].toOnChr(), 100);
    EXPECT_EQ((*genomeWindow)[8].fromOnChr(), 100); EXPECT_EQ((*genomeWindow)[8].toOnChr(), 177);
}

TEST_F(TBamFile_Test_Windows, depthPerWindow){
    // is depth per window correctly calculated?

    // 1. chr (cases 1), 2) and 3))
    EXPECT_EQ(genomeWindow->depth[0], 0.95);
    EXPECT_EQ(genomeWindow->depth[1], 0.45);
    EXPECT_EQ(genomeWindow->depth[2], 0.4);

    // 2. chr (case 4))
    EXPECT_EQ(genomeWindow->depth[3], 0.4);

    // 3. chr (case 5))
    EXPECT_EQ(genomeWindow->depth[4], 0.2);
    EXPECT_EQ(genomeWindow->depth[5], 0.);

    // 4. chr (case6))
    EXPECT_EQ(genomeWindow->depth[6], 0.);

    // 5. chr (case 7) and 8))
    EXPECT_EQ(genomeWindow->depth[7], 0.15); // only M count for depth; S, I and D of cigar are ignored
    EXPECT_EQ(genomeWindow->depth[8], 0.);
}

TEST_F(TBamFile_Test_Windows, fractionSitesNoData_PerWindow){
    // 1. chr (cases 1), 2) and 3))
    EXPECT_EQ(genomeWindow->fractionSitesNoData[0], 0.4);
    EXPECT_EQ(genomeWindow->fractionSitesNoData[1], 0.8);
    EXPECT_EQ(genomeWindow->fractionSitesNoData[2], 0.6);

    // 2. chr (case 4))
    EXPECT_EQ(genomeWindow->fractionSitesNoData[3], 0.6);

    // 3. chr (case 5))
    EXPECT_EQ(genomeWindow->fractionSitesNoData[4], 0.8);
    EXPECT_EQ(genomeWindow->fractionSitesNoData[5], 1.);

    // 4. chr (case6))
    EXPECT_EQ(genomeWindow->fractionSitesNoData[6], 1.);

    // 5. chr (case 7) and 8))
    EXPECT_EQ(genomeWindow->fractionSitesNoData[7], 0.87);
    EXPECT_EQ(genomeWindow->fractionSitesNoData[8], 1.);
}

TEST_F(TBamFile_Test_Windows, fractionDepthAtLeastTwo_PerWindow){
    // 1. chr (cases 1), 2) and 3))
    EXPECT_EQ(genomeWindow->fractionDepthAtLeastTwo[0], 0.3);
    EXPECT_EQ(genomeWindow->fractionDepthAtLeastTwo[1], 0.15);
    EXPECT_EQ(genomeWindow->fractionDepthAtLeastTwo[2], 0.);

    // 2. chr (case 4))
    EXPECT_EQ(genomeWindow->fractionDepthAtLeastTwo[3], 0.);

    // 3. chr (case 5))
    EXPECT_EQ(genomeWindow->fractionDepthAtLeastTwo[4], 0.);
    EXPECT_EQ(genomeWindow->fractionDepthAtLeastTwo[5], 0.);

    // 4. chr (case6))
    EXPECT_EQ(genomeWindow->fractionDepthAtLeastTwo[6], 0.);

    // 5. chr (case 7) and 8))
    EXPECT_EQ(genomeWindow->fractionDepthAtLeastTwo[7], 0.02);
    EXPECT_EQ(genomeWindow->fractionDepthAtLeastTwo[8], 0.);
}

TEST_F(TBamFile_Test_Windows, numSitesWithData_PerWindow){
    // 1. chr (cases 1), 2) and 3))
    EXPECT_EQ(genomeWindow->numSitesWithData[0], 60);
    EXPECT_EQ(genomeWindow->numSitesWithData[1], 20);
    EXPECT_EQ(genomeWindow->numSitesWithData[2], 20);

    // 2. chr (case 4))
    EXPECT_EQ(genomeWindow->numSitesWithData[3], 20);

    // 3. chr (case 5))
    EXPECT_EQ(genomeWindow->numSitesWithData[4], 20);
    EXPECT_EQ(genomeWindow->numSitesWithData[5], 0);

    // 4. chr (case6))
    EXPECT_EQ(genomeWindow->numSitesWithData[6], 0);

    // 5. chr (case 7) and 8))
    EXPECT_EQ(genomeWindow->numSitesWithData[7], 13);
    EXPECT_EQ(genomeWindow->numSitesWithData[8], 0);
}

TEST_F(TBamFile_Test_Windows, numReadsInWindow_PerWindow){
    // 1. chr (cases 1), 2) and 3))
    EXPECT_EQ(genomeWindow->numReadsInWindow[0], 6);
    EXPECT_EQ(genomeWindow->numReadsInWindow[1], 3);
    EXPECT_EQ(genomeWindow->numReadsInWindow[2], 1);

    // 2. chr (case 4))
    EXPECT_EQ(genomeWindow->numReadsInWindow[3], 1);

    // 3. chr (case 5))
    EXPECT_EQ(genomeWindow->numReadsInWindow[4], 1);
    EXPECT_EQ(genomeWindow->numReadsInWindow[5], 0);

    // 4. chr (case6))
    EXPECT_EQ(genomeWindow->numReadsInWindow[6], 0);

    // 5. chr (case 7) and 8))
    EXPECT_EQ(genomeWindow->numReadsInWindow[7], 3);
    EXPECT_EQ(genomeWindow->numReadsInWindow[8], 0);
}

TEST_F(TBamFile_Test_Windows, sites_getBases){
    // 1. chr (cases 1), 2) and 3))
    int c = 0;
    for (const auto& site : genomeWindow->sites[0]){ // first window
        if (c < 10)
            EXPECT_EQ(site.getBases(), "A");
        else if (c >= 10 && c < 20)
            EXPECT_EQ(site.getBases(), "AC");
        else if (c >= 20 && c < 30)
            EXPECT_EQ(site.getBases(), "CG");
        else if (c >= 30 && c < 40)
            EXPECT_EQ(site.getBases(), "G");
        else if (c >= 80 && c < 90)
            EXPECT_EQ(site.getBases(), "T");
        else if (c >= 90 && c < 95)
            EXPECT_EQ(site.getBases(), "TA");
        else if (c >= 95 && c < 100)
            EXPECT_EQ(site.getBases(), "TAC");
        else
            EXPECT_EQ(site.getBases(), "-");
        c++;
    }
    c = 0;
    for (const auto& site : genomeWindow->sites[1]){ // second window
        if (c < 10)
            EXPECT_EQ(site.getBases(), "ACG");
        else if (c >= 10 && c < 15)
            EXPECT_EQ(site.getBases(), "CG");
        else if (c >= 15 && c < 20)
            EXPECT_EQ(site.getBases(), "G");
        else
            EXPECT_EQ(site.getBases(), "-");
        c++;
    }
    c = 0;
    for (const auto& site : genomeWindow->sites[2]){ // third window
        if (c >= 20 && c < 40)
            EXPECT_EQ(site.getBases(), "T");
        else
            EXPECT_EQ(site.getBases(), "-");
        c++;
    }

    // 2. chr (case 4))
    c = 0;
    for (const auto& site : genomeWindow->sites[3]){ // 4. window
        if (c >= 10 && c < 30)
            EXPECT_EQ(site.getBases(), "A");
        else
            EXPECT_EQ(site.getBases(), "-");
        c++;
    }

    // 3. chr (case 5))
    c = 0;
    for (const auto& site : genomeWindow->sites[4]){ // 5. window
        if (c >= 10 && c < 30)
            EXPECT_EQ(site.getBases(), "C");
        else
            EXPECT_EQ(site.getBases(), "-");
        c++;
    }
    c = 0;
    for (const auto& site : genomeWindow->sites[5]){ // 6. window
        EXPECT_EQ(site.getBases(), "-");
        c++;
    }

    // 4. chr (case6))
    c = 0;
    for (const auto& site : genomeWindow->sites[6]){ // 7. window
        EXPECT_EQ(site.getBases(), "-");
        c++;
    }

    // 5. chr (case 7) and 8))
    c = 0;
    for (const auto& site : genomeWindow->sites[7]){ // 8. window
        if (c >= 0 && c < 4)
            EXPECT_EQ(site.getBases(), "A");
        else if (c == 4)
            EXPECT_EQ(site.getBases(), "AC");
        else if (c >= 5 && c < 8)
            EXPECT_EQ(site.getBases(), "C");
        else if (c == 8)
            EXPECT_EQ(site.getBases(), "CG");
        else if (c >= 9 && c < 13)
            EXPECT_EQ(site.getBases(), "G");
        else
            EXPECT_EQ(site.getBases(), "-");
        c++;
    }
    c = 0;
    for (const auto& site : genomeWindow->sites[8]){ // 9. window
        EXPECT_EQ(site.getBases(), "-");
    }
}

TEST_F(TBamFile_Test_Windows, sites_getQualities){
    // 1. chr (cases 1), 2) and 3))
    int c = 0;
    for (const auto& site : genomeWindow->sites[0]){ // first window
        if (c < 10)
            EXPECT_EQ(site.getQualities(), "\"");
        else if (c >= 10 && c < 20)
            EXPECT_EQ(site.getQualities(), "\"#");
        else if (c >= 20 && c < 30)
            EXPECT_EQ(site.getQualities(), "#$");
        else if (c >= 30 && c < 40)
            EXPECT_EQ(site.getQualities(), "$");
        else if (c >= 80 && c < 90)
            EXPECT_EQ(site.getQualities(), "%");
        else if (c >= 90 && c < 95)
            EXPECT_EQ(site.getQualities(), "%&");
        else if (c >= 95 && c < 100)
            EXPECT_EQ(site.getQualities(), "%&'");
        else
            EXPECT_EQ(site.getQualities(), "-");
        c++;
    }
    c = 0;
    for (const auto& site : genomeWindow->sites[1]){ // second window
        if (c < 10)
            EXPECT_EQ(site.getQualities(), "&'(");
        else if (c >= 10 && c < 15)
            EXPECT_EQ(site.getQualities(), "'(");
        else if (c >= 15 && c < 20)
            EXPECT_EQ(site.getQualities(), "(");
        else
            EXPECT_EQ(site.getQualities(), "-");
        c++;
    }
    c = 0;
    for (const auto& site : genomeWindow->sites[2]){ // third window
        if (c >= 20 && c < 40)
            EXPECT_EQ(site.getQualities(), ")");
        else
            EXPECT_EQ(site.getQualities(), "-");
        c++;
    }

    // 2. chr (case 4))
    c = 0;
    for (const auto& site : genomeWindow->sites[3]){ // 4. window
        if (c >= 10 && c < 30)
            EXPECT_EQ(site.getQualities(), "*");
        else
            EXPECT_EQ(site.getQualities(), "-");
        c++;
    }

    // 3. chr (case 5))
    c = 0;
    for (const auto& site : genomeWindow->sites[4]){ // 5. window
        if (c >= 10 && c < 30)
            EXPECT_EQ(site.getQualities(), "+");
        else
            EXPECT_EQ(site.getQualities(), "-");
        c++;
    }
    c = 0;
    for (const auto& site : genomeWindow->sites[5]){ // 6. window
        EXPECT_EQ(site.getQualities(), "-");
        c++;
    }

    // 4. chr (case6))
    c = 0;
    for (const auto& site : genomeWindow->sites[6]){ // 7. window
        EXPECT_EQ(site.getQualities(), "-");
        c++;
    }

    // 5. chr (case 7) and 8))
    c = 0;
    for (const auto& site : genomeWindow->sites[7]){ // 8. window
        if (c >= 0 && c < 4)
            EXPECT_EQ(site.getQualities(), "\"");
        else if (c == 4)
            EXPECT_EQ(site.getQualities(), "\"#");
        else if (c >= 5 && c < 8)
            EXPECT_EQ(site.getQualities(), "#");
        else if (c == 8)
            EXPECT_EQ(site.getQualities(), "#$");
        else if (c >= 9 && c < 13)
            EXPECT_EQ(site.getQualities(), "$");
        else
            EXPECT_EQ(site.getQualities(), "-");
        c++;
    }
    c = 0;
    for (const auto& site : genomeWindow->sites[8]){ // 9. window
        EXPECT_EQ(site.getQualities(), "-");
    }
}

//-------------------------------------------------------------
// TBamFile - filters
//-------------------------------------------------------------

using coretools::TCountDistribution;
using coretools::TCountDistributionVector;

class TBamFilter : public GenomeTasks::old::TGenome_filtered {
    // class very similar to TBamDiagnoser, but inherits from TGenome_filtered -> can apply all filters
public:
    // distributions
    TCountDistribution<> totalReads;

    TCountDistribution<> duplicates;
    TCountDistribution<> improperPairs;
    TCountDistribution<> failedQC;
    TCountDistribution<> unmapped;
    TCountDistribution<> secondary;
    TCountDistribution<> supplementary;
    TCountDistribution<> longerThanFragmentLength;
    TCountDistribution<> fwdStrand;
    TCountDistribution<> revStrand;
    TCountDistribution<> firstMate;
    TCountDistribution<> secondMate;
    std::vector<std::string> names; // for blacklist

    TCountDistributionVector<> readLength;
    TCountDistributionVector<> usableLength;
    TCountDistributionVector<> softClippedLength;
    TCountDistributionVector<> mappingQuality;
    TCountDistributionVector<> fragmentLength;

    TCountDistribution<> readGroup;
    TCountDistribution<> refIDs;

    GenomeTasks::TQualityFilter qualFilter;
    std::vector<std::string> readGroupNames;

    void _handleAlignment() override {
        //get read group
        uint16_t curReadGroup = _bamFile.curReadGroupID();

        //add to counters
        totalReads.add(curReadGroup);

        duplicates.add(curReadGroup, _bamFile.curIsDuplicate());
        improperPairs.add(curReadGroup, _bamFile.curIsPaired() && !_bamFile.curIsProperPair());
        failedQC.add(curReadGroup, _bamFile.curIsFailedQC());
        unmapped.add(curReadGroup, !_bamFile.curIsMapped());
        secondary.add(curReadGroup, _bamFile.curIsSecondary());
        supplementary.add(curReadGroup, _bamFile.curIsSupplementary());
        longerThanFragmentLength.add(curReadGroup, _bamFile.curIsLongerThanFragment());
        fwdStrand.add(curReadGroup, !_bamFile.curIsReverseStrand());
        revStrand.add(curReadGroup, _bamFile.curIsReverseStrand());
        firstMate.add(curReadGroup, _bamFile.curIsFirstMate());
        secondMate.add(curReadGroup, _bamFile.curIsSecondMate());
        names.push_back(_bamFile.curName());

        readLength.add(curReadGroup, _bamFile.curCIGAR().lengthRead());
        softClippedLength.add(curReadGroup, _bamFile.curCIGAR().lengthSoftClipped());
        mappingQuality.add(curReadGroup, _bamFile.curMappingQuality());

        readGroup.add(curReadGroup);
        refIDs.add(_bamFile.curPosition().refID());

        //fragment length: only once
        if(!_bamFile.curIsReverseStrand()){
            fragmentLength.add(curReadGroup, _bamFile.curFragmentLength());
        }
    }

    TBamFilter() {
        //settings
        _bamFile.readGroups().fillVectorWithNames(readGroupNames);
    };

    void filter(){
        //initialize counters
        size_t numRG = _bamFile.readGroups().size();

        // resize distributions
        totalReads.resize(numRG);

        duplicates.resize(numRG);
        improperPairs.resize(numRG);
        failedQC.resize(numRG);
        unmapped.resize(numRG);
        secondary.resize(numRG);
        supplementary.resize(numRG);
        longerThanFragmentLength.resize(numRG);
        fwdStrand.resize(numRG);
        revStrand.resize(numRG);
        firstMate.resize(numRG);
        secondMate.resize(numRG);

        readLength.resize(numRG);
        usableLength.resize(numRG);
        softClippedLength.resize(numRG);
        mappingQuality.resize(numRG);
        fragmentLength.resize(numRG);

        readGroup.resize(numRG);

        //now parse through bam file
        _traverseBAMPassedQC();
    };
};

class TBamFilter_Test : public ::testing::Test {
protected:
    void _disableAllFilters(){
        // keep all reads
		parameters().clear();
        parameters().add("keepDuplicates");
        parameters().add("keepImproperPairs");
        parameters().add("keepUnmappedReads");
        parameters().add("keepFailedQC");
        parameters().add("keepSecondaryReads");
        parameters().add("keepSupplementaryReads");
        parameters().add("keepReadsLongerThanFragment");
    }

public:
    std::unique_ptr<BAM::TTestBamFile> outputBam;
    std::string filename = "testBAM.bam";
    uint16_t numReads = 2000;

    std::unique_ptr<TBamFilter> bamFilter;

    void write(bool paired){
        //settings
        std::vector<size_t> chrLength = {1000, 2000, 3000};
        size_t numReadGroups = 5;

        //open BAM file for writing
        if (paired)
            outputBam = std::make_unique<BAM::TTestBamFilePairedEnd>(filename, chrLength, numReadGroups);
        else
            outputBam = std::make_unique<BAM::TTestBamFile>(filename, chrLength, numReadGroups);

        //write alignments
        outputBam->writeDummyAlignments(numReads, true);

        outputBam->closeOutput();
    }

    void read(){ // called by every test, as parameters (= filters) are different
        // set bam parameters
        parameters().add("bam", filename);

        // create instance of TBamFilter
        bamFilter = std::make_unique<TBamFilter>();
        bamFilter->filter();
    }
};

TEST_F(TBamFilter_Test, maxReadLength){
    write(false);

    // 1) filter: maxReadLength
    _disableAllFilters();
    parameters().add("filterReadLength", "[20,30]");

    read();

    //count number of simulated alignments outside this range
	size_t numAligmentsOutsideRange = 0;
	for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
		if(a->cigar().lengthRead() < 20 || a->cigar().lengthRead() > 30){
			++numAligmentsOutsideRange;
		}
	}

	std::cout << "numAligmentsOutsideRange = " << numAligmentsOutsideRange << std::endl;
	std::cout << "read length = " << bamFilter->readLength.counts() << std::endl;

	coretools::TOutputFile out("RL.out");
	bamFilter->readLength.write(out);

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numAligmentsOutsideRange);
};

TEST_F(TBamFilter_Test, keepDuplicates){
    write(false);
    _disableAllFilters();
    // 2) do not filter: 'keepDuplicates'
    parameters().add("keepDuplicates");
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads);
    EXPECT_TRUE(bamFilter->duplicates.counts() > 0);
};

TEST_F(TBamFilter_Test, doNotKeepDuplicates){
    write(false);
	parameters().clear();
    // 2) filter: do not specify 'keepDuplicates'
    parameters().add("keepImproperPairs");
    parameters().add("keepUnmappedReads");
    parameters().add("keepFailedQC");
    parameters().add("keepSecondaryReads");
    parameters().add("keepSupplementaryReads");
    parameters().add("keepReadsLongerThanFragment");
    read();

    // count number of duplicates in simulated alignments
    size_t numDupWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        // get sam flag of alignment as int and convert
        BAM::TSamFlags samFlags(a->flags());
        if (samFlags.isDuplicate())
            numDupWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numDupWritten);
    EXPECT_TRUE(bamFilter->duplicates.counts() == 0);
}

TEST_F(TBamFilter_Test, filterSoftClips){
    write(false);
    _disableAllFilters();
    // 3) filter: 'filterSoftClips'
    parameters().add("filterSoftClips");
    read();

    // count number of soft clips in simulated alignments
    size_t numSoftClipsWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        if (a->cigar().lengthSoftClipped() > 0)
            numSoftClipsWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numSoftClipsWritten);
    EXPECT_TRUE(bamFilter->softClippedLength.counts() == 0);
}

TEST_F(TBamFilter_Test, doNotfilterSoftClips){
    write(false);
    _disableAllFilters();
    // 3) do not filter: 'filterSoftClips'
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads);
    EXPECT_TRUE(bamFilter->softClippedLength.counts() > 0);
}

TEST_F(TBamFilter_Test, keepImproperPairs){
    write(true);
    _disableAllFilters();
    // 2) do not filter: 'keepImproperPairs'
    parameters().add("keepImproperPairs");
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads);
    EXPECT_TRUE(bamFilter->improperPairs.counts() > 0);
}

TEST_F(TBamFilter_Test, doNotKeepImproperPairs){
    write(true);
	parameters().clear();
    // 2) filter: do not specify 'keepImproperPairs'
    parameters().add("keepDuplicates");
    parameters().add("keepUnmappedReads");
    parameters().add("keepFailedQC");
    parameters().add("keepSecondaryReads");
    parameters().add("keepSupplementaryReads");
    parameters().add("keepReadsLongerThanFragment");
    read();

    // count number of improper pairs in simulated alignments
    size_t numImproperWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        // get sam flag of alignment as int and convert
        BAM::TSamFlags samFlags(a->flags());
        if (!samFlags.isProperPair())
            numImproperWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numImproperWritten);
    EXPECT_TRUE(bamFilter->improperPairs.counts() == 0);
}

TEST_F(TBamFilter_Test, keepUnmappedReads){
    write(false);
    _disableAllFilters();
    // 3) do not filter: 'keepUnmappedReads'
    parameters().add("keepUnmappedReads");
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads);
    EXPECT_TRUE(bamFilter->unmapped.counts() > 0);
}

TEST_F(TBamFilter_Test, doNotKeepUnmappedReads){
    write(false);
	parameters().clear();
    // 3) filter: do not specify 'keepUnmappedReads'
    parameters().add("keepDuplicates");
    parameters().add("keepImproperPairs");
    parameters().add("keepFailedQC");
    parameters().add("keepSecondaryReads");
    parameters().add("keepSupplementaryReads");
    parameters().add("keepReadsLongerThanFragment");
    read();

    // count number of unmapped reads in simulated alignments
    size_t numUnmappedWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        // get sam flag of alignment as int and convert
        BAM::TSamFlags samFlags(a->flags());
        if (samFlags.isUnmapped())
            numUnmappedWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numUnmappedWritten);
    EXPECT_TRUE(bamFilter->unmapped.counts() == 0);
}

TEST_F(TBamFilter_Test, keepFailedQC){
    write(false);
    _disableAllFilters();
    // 3) do not filter: 'keepFailedQC'
    parameters().add("keepFailedQC");
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads);
    EXPECT_TRUE(bamFilter->failedQC.counts() > 0);
}

TEST_F(TBamFilter_Test, doNotKeepFailedQC){
    write(false);
	parameters().clear();
    // 3) filter: do not specify 'keepFailedQC'
    parameters().add("keepDuplicates");
    parameters().add("keepImproperPairs");
    parameters().add("keepUnmappedReads");
    parameters().add("keepSecondaryReads");
    parameters().add("keepSupplementaryReads");
    parameters().add("keepReadsLongerThanFragment");
    read();

    // count number of failed QC reads in simulated alignments
    size_t numFailedQCWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        // get sam flag of alignment as int and convert
        BAM::TSamFlags samFlags(a->flags());
        if (samFlags.isQCFailed())
            numFailedQCWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numFailedQCWritten);
    EXPECT_TRUE(bamFilter->failedQC.counts() == 0);
}

TEST_F(TBamFilter_Test, keepSecondaryReads){
    write(false);
    _disableAllFilters();
    // 4) do not filter: 'keepSecondaryReads'
    parameters().add("keepSecondaryReads");
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads);
    EXPECT_TRUE(bamFilter->secondary.counts() > 0);
}

TEST_F(TBamFilter_Test, doNotKeepSecondaryReads){
    write(false);
	parameters().clear();
    // 4) filter: do not specify 'keepSecondaryReads'
    parameters().add("keepDuplicates");
    parameters().add("keepImproperPairs");
    parameters().add("keepUnmappedReads");
    parameters().add("keepFailedQC");
    parameters().add("keepSupplementaryReads");
    parameters().add("keepReadsLongerThanFragment");
    read();

    // count number of secondary reads in simulated alignments
    size_t numSecReadsWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        // get sam flag of alignment as int and convert
        BAM::TSamFlags samFlags(a->flags());
        if (samFlags.isSecondary())
            numSecReadsWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numSecReadsWritten);
    EXPECT_TRUE(bamFilter->secondary.counts() == 0);
}

TEST_F(TBamFilter_Test, keepSupplementaryReads){
    write(false);
    _disableAllFilters();
    // 5) do not filter: 'keepSupplementaryReads'
    parameters().add("keepSupplementaryReads");
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads);
    EXPECT_TRUE(bamFilter->supplementary.counts() > 0);
}

TEST_F(TBamFilter_Test, doNotKeepSupplementaryReads){
    write(false);
	parameters().clear();
    // 5) filter: do not specify 'keepSupplementaryReads'
    parameters().add("keepDuplicates");
    parameters().add("keepImproperPairs");
    parameters().add("keepUnmappedReads");
    parameters().add("keepFailedQC");
    parameters().add("keepSecondaryReads");
    parameters().add("keepReadsLongerThanFragment");
    read();

    // count number of supplementary reads in simulated alignments
    size_t numSupplementaryWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        // get sam flag of alignment as int and convert
        BAM::TSamFlags samFlags(a->flags());
        if (samFlags.isSupplementary())
            numSupplementaryWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numSupplementaryWritten);
    EXPECT_TRUE(bamFilter->supplementary.counts() == 0);
}

TEST_F(TBamFilter_Test, keepReadsLongerThanFragment){
    write(true);
    _disableAllFilters();
    // 6) do not filter: 'keepReadsLongerThanFragment'
    parameters().add("keepReadsLongerThanFragment");
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads);
    EXPECT_TRUE(bamFilter->longerThanFragmentLength.counts() > 0);
}

TEST_F(TBamFilter_Test, filterReadsLongerThanFragment){
    write(true);
	parameters().clear();
    // 6) filter: do not specify 'keepReadsLongerThanFragment'
    parameters().add("keepDuplicates");
    parameters().add("keepImproperPairs");
    parameters().add("keepUnmappedReads");
    parameters().add("keepFailedQC");
    parameters().add("keepSecondaryReads");
    parameters().add("keepSupplementaryReads");
    parameters().add("filterReadsLongerThanFragment");
    read();

    // count number of reads that are longer than fragments in simulated alignments
    size_t numLongerWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        size_t insertSize;
        if (a->matePosition() > a->position())
            insertSize = a->matePosition() - a->position();
        else insertSize = a->position() - a->matePosition();
        if (a->isProperPair() && insertSize < a->cigar().lengthAligned())
            numLongerWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numLongerWritten);
    EXPECT_TRUE(bamFilter->longerThanFragmentLength.counts() == 0);
}

TEST_F(TBamFilter_Test, keepOnlyFwd){
    write(false);
    _disableAllFilters();
    // 7) filter: 'keepOnlyFwd', remove all reverse reads
    parameters().add("keepOnlyFwd");
    read();

    // count number of reverse reads in simulated alignments
    size_t numRevWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        if (a->isReverseStrand())
            numRevWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numRevWritten);
    EXPECT_TRUE(bamFilter->fwdStrand.counts() > 0);
    EXPECT_TRUE(bamFilter->revStrand.counts() == 0);
}

TEST_F(TBamFilter_Test, keepOnlyRev){
    write(false);
    _disableAllFilters();
    // 7) filter: 'keepOnlyRev', remove all forward reads
    parameters().add("keepOnlyRev");
    read();

    // count number of forward reads in simulated alignments
    size_t numFwdWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        if (!a->isReverseStrand())
            numFwdWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numFwdWritten);
    EXPECT_TRUE(bamFilter->revStrand.counts() > 0);
    EXPECT_TRUE(bamFilter->fwdStrand.counts() == 0);
}

TEST_F(TBamFilter_Test, keepOnlyFirst){
    write(true);
    _disableAllFilters();
    // 8) filter: 'keepOnlyFirst', remove all second mates
    parameters().add("keepOnlyFirst");
    read();

    // count number of second mates in simulated alignments
    size_t numSecondWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        // get sam flag of alignment as int and convert
        BAM::TSamFlags samFlags(a->flags());
        if (samFlags.isSecondMate())
            numSecondWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numSecondWritten);
    EXPECT_TRUE(bamFilter->firstMate.counts() > 0);
    EXPECT_TRUE(bamFilter->secondMate.counts() == 0);
}

TEST_F(TBamFilter_Test, keepOnlySecond){
    write(true);
    _disableAllFilters();
    // 8) filter: 'keepOnlySecond', remove all first mates
    parameters().add("keepOnlySecond");
    read();

    // count number of first mates in simulated alignments
    size_t numFirstWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        // get sam flag of alignment as int and convert
        BAM::TSamFlags samFlags(a->flags());
        if (samFlags.isFirstMate())
            numFirstWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numFirstWritten);
    EXPECT_TRUE(bamFilter->secondMate.counts() > 0);
    EXPECT_TRUE(bamFilter->firstMate.counts() == 0);
}

TEST_F(TBamFilter_Test, blacklist){
    write(false);
    _disableAllFilters();
    // 9) filter: 'blacklist'

    // generate blacklist
    coretools::TOutputFile blackList("blacklist.txt", "");
    for (int i = 0; i < 200; i++){ // add first 200 reads to blacklist
        blackList.writeln("alignment_", i);
    }
    blackList.close();

    parameters().add("blacklist", "blacklist.txt");
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - 200);
    for (auto & it : bamFilter->names){
        // all alignments 0-199 should have been removed
       int alignmentNum = coretools::str::fromString<int>(coretools::str::extractAfter(it, "_"));
       EXPECT_TRUE(alignmentNum >= 200);
    }

    // remove blacklist file
    remove("blacklist.txt");
}

TEST_F(TBamFilter_Test, MQ){
    write(false);
    _disableAllFilters();
    // 10) filter: 'minMQ' and 'maxMQ'
    parameters().add("filterMQ", "[80,100]");
    read();

    // count number of reads with lower and higher MQ in simulated alignments
    size_t numReadsWithOutsideMQWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        if (a->mappingQuality().get() < 80 || a->mappingQuality().get() > 100)
            numReadsWithOutsideMQWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numReadsWithOutsideMQWritten);

    for (size_t id = 0; id < bamFilter->mappingQuality.size(); id++){
        EXPECT_TRUE(bamFilter->mappingQuality[id].min() >= 80);
        EXPECT_TRUE(bamFilter->mappingQuality[id].max() <= 100);
    }
}

TEST_F(TBamFilter_Test, fragmentLength){
    write(true);
    _disableAllFilters();
    // 11) filter: 'minFragmentLength' and 'maxFragmentLength'
    parameters().add("filterFragmentLength", "[80,100]");
    read();


    // count number of reads with lower fragment length in simulated alignments
    size_t numReadsWithOutsideFragLengthWritten = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        // calculate insert size
        size_t insertSize;
        if (a->matePosition() > a->position())
            insertSize = a->matePosition() - a->position();
        else insertSize = a->position() - a->matePosition();
        if (a->isProperPair()) {
            if (insertSize + a->cigar().lengthInserted() - a->cigar().lengthDeleted() < 80 || insertSize + a->cigar().lengthInserted() - a->cigar().lengthDeleted() > 100)
                numReadsWithOutsideFragLengthWritten++;
        } else if (a->cigar().lengthSequenced() < 80 || a->cigar().lengthSequenced() > 100)
            numReadsWithOutsideFragLengthWritten++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numReadsWithOutsideFragLengthWritten);

    for (size_t id = 0; id < bamFilter->mappingQuality.size(); id++){
        if (!bamFilter->fragmentLength[id].empty()) {
            EXPECT_TRUE(bamFilter->fragmentLength[id].min() >= 80);
            EXPECT_TRUE(bamFilter->fragmentLength[id].max() <= 100);
        }
    }
}

TEST_F(TBamFilter_Test, readGroups){
    write(false);
    _disableAllFilters();
    // 12) filter: 'readGroup'
    parameters().add("readGroup", "ReadGroup0,ReadGroup1,ReadGroup2");
    read();

    // count number of reads with read group ID 3 or 4 in simulated alignments
    size_t numReadGroups34 = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        if (a->readGroupId() == 3 || a->readGroupId() == 4)
            numReadGroups34++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numReadGroups34);

    EXPECT_TRUE(bamFilter->readGroup[0] > 0);
    EXPECT_TRUE(bamFilter->readGroup[1] > 0);
    EXPECT_TRUE(bamFilter->readGroup[2] > 0);
    EXPECT_TRUE(bamFilter->readGroup[3] == 0); // filtered out
    EXPECT_TRUE(bamFilter->readGroup[4] == 0); // filtered out
}

TEST_F(TBamFilter_Test, limitReads){
    write(false);
    _disableAllFilters();
    // 13) filter: 'limitReads'
    parameters().add("limitReads", "50");
    read();

    EXPECT_TRUE(bamFilter->totalReads.counts() == 50);
}

TEST_F(TBamFilter_Test, chr){
    write(false);
    _disableAllFilters();
    // 14) filter: 'chr'
    parameters().add("chr", "Chr1,Chr2");
    read();

    TCountDistribution writtenRefIDs;
    for (auto it = outputBam->beginWrittenAlignments(); it != outputBam->endWrittenAlignments(); it++){
        writtenRefIDs.add(it->refID());
    }

    // count number of reads with chromosome 3 in simulated alignments
    size_t numReadsChr3Written = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        if (a->refID() == 2)
            numReadsChr3Written++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numReadsChr3Written);

    EXPECT_EQ(writtenRefIDs[0], bamFilter->refIDs[0]);
    EXPECT_EQ(writtenRefIDs[1], bamFilter->refIDs[1]);
    EXPECT_NE(writtenRefIDs[2], bamFilter->refIDs[2]); // filtered out while reading -> expect not equal to written
}

TEST_F(TBamFilter_Test, limitChr){
    write(false);
    _disableAllFilters();
    // 14) filter: 'chr'
    parameters().add("limitChr", "Chr2");
    read();

    TCountDistribution writtenRefIDs;
    for (auto it = outputBam->beginWrittenAlignments(); it != outputBam->endWrittenAlignments(); it++){
        writtenRefIDs.add(it->refID());
    }

    // count number of reads with chromosome 3 in simulated alignments
    size_t numReadsChr3Written = 0;
    for (auto a = outputBam->beginWrittenAlignments(); a != outputBam->endWrittenAlignments(); a++){
        if (a->refID() == 2)
            numReadsChr3Written++;
    }

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads - numReadsChr3Written);

    EXPECT_EQ(writtenRefIDs[0], bamFilter->refIDs[0]);
    EXPECT_EQ(writtenRefIDs[1], bamFilter->refIDs[1]);
    EXPECT_NE(writtenRefIDs[2], bamFilter->refIDs[2]); // filtered out while reading -> expect not equal to written
}

TEST_F(TBamFilter_Test, keepAllReads){
    write(false);
    _disableAllFilters();
    // 14) filter: 'keepAllReads'
    parameters().add("keepAllReads");
    read();

    EXPECT_EQ(bamFilter->totalReads.counts(), numReads);
}
