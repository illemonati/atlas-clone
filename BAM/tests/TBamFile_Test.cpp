
#include "TBamFile.h"
#include "TTestBamFile.h"
#include "TAlignment.h"
#include "TGenome.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

//-------------------------------------------------------------
// TBamFile - simple writing and reading
//-------------------------------------------------------------

class TBamFile_Test_ReadWrite : public ::testing::Test {
protected:
    std::string _filename = "testBAM.bam";
    TLog _logfile;
    GenotypeLikelihoods::TGenotypeMap genoMap;
    BAM::TQualityMap qualMap;

public:
    std::unique_ptr<TestUtilities::TTestBamFile> outputBam;
    std::unique_ptr<BAM::TBamFile> inputBam;

    void write(){
        //settings
        std::vector<uint32_t> chrLength = {1000, 2000, 3000};
        uint32_t numReadGroups = 2;

        //open BAM file for writing
        outputBam = std::make_unique<TestUtilities::TTestBamFile>(_filename, chrLength, numReadGroups);

        //write alignments
        outputBam->writeDummyAlignments(100);
        outputBam->closeOutput();
    }

    void read(){
        //open BAM file for reading
        inputBam = std::make_unique<BAM::TBamFile>();
        inputBam->open(_filename, true, &_logfile);
    }

    void SetUp() override{
        write();
        read();
    }

    void TearDown() override {};
};

TEST_F(TBamFile_Test_ReadWrite, chromosomes){
    // get chromosomes
    BAM::TChromosomes chromosomesWritten = outputBam->chromosomes();
    BAM::TChromosomes chromosomesRead = inputBam->chromosomes();

    // compare size
    EXPECT_EQ(chromosomesWritten.size(), chromosomesRead.size());

    // compare values of iterators
    auto itRead = chromosomesRead.begin();
    for (auto itWritten = chromosomesWritten.begin(); itWritten != chromosomesWritten.end(); itWritten++, itRead++){
        EXPECT_EQ(itWritten->name, itRead->name);
        EXPECT_EQ(itWritten->length, itRead->length);
        EXPECT_EQ(itWritten->refID(), itRead->refID());
        EXPECT_EQ(itWritten->chrStart.position(), itRead->chrStart.position());
        EXPECT_EQ(itWritten->chrEnd.position(), itRead->chrEnd.position());
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
        for (uint32_t i = 0; i < alignmentWritten->parsedLength(); i++){
            EXPECT_EQ(alignmentWritten->isAlignedAtInternalPos(i), alignmentRead.isAlignedAtInternalPos(i));
            EXPECT_EQ(alignmentWritten->referenceAtInternalPos(i), alignmentRead.referenceAtInternalPos(i));
            EXPECT_EQ(alignmentWritten->positionInRef(i), alignmentRead.positionInRef(i));
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
        EXPECT_EQ(alignmentWritten->sequence(genoMap, qualMap), alignmentRead.sequence(genoMap, qualMap));
        EXPECT_EQ(alignmentWritten->qualities(genoMap, qualMap), alignmentRead.qualities(genoMap, qualMap));
        EXPECT_EQ(alignmentWritten->isReverseStrand(), alignmentRead.isReverseStrand());
        EXPECT_EQ(alignmentWritten->isPaired(), alignmentRead.isPaired());
        EXPECT_EQ(alignmentWritten->isProperPair(), alignmentRead.isProperPair());
        EXPECT_EQ(alignmentWritten->isParsed(), alignmentRead.isParsed());

        // iterate over bases of TAlignment
        auto baseRead = alignmentRead.begin();
        for (auto baseWritten = alignmentWritten->begin(); baseWritten != alignmentWritten->end(); baseWritten++, baseRead++){
            // all attributes of TBase
            EXPECT_EQ(baseWritten->originalQuality_phredInt, baseRead->originalQuality_phredInt);
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
        }
        EXPECT_EQ(baseRead, alignmentRead.end());

        alignmentWritten++;
    }

    EXPECT_EQ(alignmentWritten, outputBam->endWrittenAlignments());
}


//-------------------------------------------------------------
// TBamFile - windows
//-------------------------------------------------------------

class TGenomeWindow_Test : public GenomeTasks::TGenome_windows {
protected:
    std::vector<BAM::TGenomeWindow> _windows_visited;

    void _handleWindow() override{
        _windows_visited.emplace_back(this->_window);
        depth.emplace_back(_window.depth());
    };

public:
    std::vector<double> depth;

    TGenomeWindow_Test(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator) : GenomeTasks::TGenome_windows(Params, Logfile, RandomGenerator) {};

    void traverse(){
        _traverseBAMWindows();
    }

    // size
    uint32_t numWindows(){ return _windows_visited.size(); };
    // loop
    std::vector<BAM::TGenomeWindow >::iterator begin(){ return _windows_visited.begin(); };
    std::vector<BAM::TGenomeWindow >::iterator end(){ return _windows_visited.end(); };
    // access
    BAM::TGenomeWindow& operator[](uint32_t pos){ return _windows_visited[pos]; };
};

class TBamFile_Test_Windows : public ::testing::Test {
protected:
    TLog _logfile;
    TRandomGenerator _randomGenerator;
    TParameters _parameters;

public:
    std::unique_ptr<TestUtilities::TTestBamFile> outputBam;
    std::unique_ptr<TGenomeWindow_Test> genomeWindow;
    std::string filename = "testBAM.bam";

    void write(){
        //settings
        std::vector<uint32_t> chrLength = {250, 50, 199, 80, 77};
        uint32_t numReadGroups = 2;

        //open BAM file for writing
        outputBam = std::make_unique<TestUtilities::TTestBamFile>(filename, chrLength, numReadGroups);

        //write alignments

        // 1) overlapping alignments inside one window, cigar is all M
        BAM::TCigar onlyMCigar;
        onlyMCigar.add('M', 20);
        outputBam->writeDummyAlignment('A', '1', BAM::TGenomePosition(0, 0), onlyMCigar);
        outputBam->writeDummyAlignment('C', '2', BAM::TGenomePosition(0, 10), onlyMCigar);
        outputBam->writeDummyAlignment('G', '3', BAM::TGenomePosition(0, 20), onlyMCigar);

        // 2) alignments overlap 2 windows
        outputBam->writeDummyAlignment('T', '4', BAM::TGenomePosition(0, 80), onlyMCigar);
        outputBam->writeDummyAlignment('A', '5', BAM::TGenomePosition(0, 90), onlyMCigar);
        outputBam->writeDummyAlignment('C', '6', BAM::TGenomePosition(0, 95), onlyMCigar);
        outputBam->writeDummyAlignment('G', '7', BAM::TGenomePosition(0, 100), onlyMCigar);

        // 3) one alignment inside 1 window
        outputBam->writeDummyAlignment('T', '8', BAM::TGenomePosition(0, 220), onlyMCigar);

        // 4) only 1 window per chromosome
        outputBam->writeDummyAlignment('A', '9', BAM::TGenomePosition(1, 10), onlyMCigar);

        // 5) empty window
        outputBam->writeDummyAlignment('C', '0', BAM::TGenomePosition(2, 10), onlyMCigar);

        // 6) empty chromosome

        // 7) cigar strings more complicated (with soft-clips, insertions and deletions)
        BAM::TCigar mixedCigar;
        mixedCigar.add('S', 3);
        mixedCigar.add('M', 5);
        mixedCigar.add('D', 5);
        mixedCigar.add('I', 5);
        mixedCigar.add('S', 2);
        outputBam->writeDummyAlignment('A', '1', BAM::TGenomePosition(4, 0), mixedCigar);
        outputBam->writeDummyAlignment('C', '2', BAM::TGenomePosition(4, 4), mixedCigar);
        outputBam->writeDummyAlignment('G', '3', BAM::TGenomePosition(4, 8), mixedCigar);


        outputBam->closeOutput();
    }

    void read(){
        // first set window size in parameters
        _parameters.addParameter("window", "100");
        _parameters.addParameter("bam", filename);
        _parameters.addParameter("maxReadLength", "20");

        // create instance of TGenomeWindow
        genomeWindow = std::make_unique<TGenomeWindow_Test>(_parameters, &_logfile, &_randomGenerator);

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
    EXPECT_EQ(genomeWindow->numWindows(),6); // 6 of 8 windows contain data and are stored
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

    // 5. chr (case 7))
    EXPECT_EQ((*genomeWindow)[5].refID(), 2);
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

    // 5. chr (case 7))
    EXPECT_EQ((*genomeWindow)[5].fromOnChr(), 0); EXPECT_EQ((*genomeWindow)[5].toOnChr(), 77);
}

TEST_F(TBamFile_Test_Windows, depthPerWindow){
    // do windows store correct positions? I.e., if chromosome is shorter than window, is window correctly adjusted?

    // 1. chr (cases 1), 2) and 3))
    EXPECT_EQ(genomeWindow->depth[0], 0.95);
    EXPECT_EQ(genomeWindow->depth[1], 0.45);
    EXPECT_EQ(genomeWindow->depth[2], 0.4);

    // 2. chr (case 4))
    EXPECT_EQ(genomeWindow->depth[3], 0.4);

    // 3. chr (case 5))
    EXPECT_EQ(genomeWindow->depth[4], 0.2);

    // 5. chr (case 7))
    EXPECT_EQ(genomeWindow->depth[5], 0.1948052);
}