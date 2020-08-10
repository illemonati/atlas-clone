
#include "TBamFile.h"
#include "TTestBamFile.h"
#include "TAlignment.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

//-------------------------------------------------------------
// TBamFile
//-------------------------------------------------------------

class TBamFile_Test_Easy : public ::testing::Test {
private:
    std::string _filename = "testBAM.bam";
    TLog _logfile;

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
        outputBam->writeDummyAlignments(10);
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

TEST_F(TBamFile_Test_Easy, chromosomes){
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

TEST_F(TBamFile_Test_Easy, samHeader){
    // get header
    BAM::TSamHeader headerWritten = outputBam->header();
    BAM::TSamHeader headerRead = inputBam->samHeader();

    // compare values
    EXPECT_EQ(headerWritten.version(), headerRead.version());
    EXPECT_EQ(headerWritten.grouping(), headerRead.grouping());
    EXPECT_EQ(headerWritten.sortOrder(), headerRead.sortOrder());
    EXPECT_EQ(headerWritten.subSorting(), headerRead.subSorting());
}

TEST_F(TBamFile_Test_Easy, readGroups){
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

TEST_F(TBamFile_Test_Easy, alignments){
    // read
    auto it = outputBam->beginWrittenAlignments();
    while (inputBam->readNextAlignment()){
        EXPECT_EQ(it->name(), inputBam->curName());
        EXPECT_EQ(it->readGroupId(), inputBam->curReadGroupID());
        EXPECT_EQ(it->mappingQuality(), inputBam->curMappingQuality());

        it++;
    }

    EXPECT_EQ(it, outputBam->endWrittenAlignments());
}