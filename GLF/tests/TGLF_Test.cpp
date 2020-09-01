//
// Created by caduffm on 9/1/20.
//

#include "TGLF.h"
#include "TTestGLFFile.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

//-------------------------------------------------------------
// TGLF_Test_ReadWrite - simple writing and reading
//-------------------------------------------------------------

class TGLF_Test_ReadWrite : public ::testing::Test {
protected:
    std::string _filename = "testGLF.glf";

public:
    std::unique_ptr<TestUtilities::TTestGLFFile> outputGLF;
    std::unique_ptr<TGlfReader> inputGLF;

    // store stuff from input GLF to later compare with output
    std::vector<uint32_t> positions;
    std::vector<uint8_t> depths;

    void write(){
        //settings
        std::vector<uint32_t> chrLength = {100, 200, 300};

        //open GLF file for writing
        outputGLF = std::make_unique<TestUtilities::TTestGLFFile>(_filename, chrLength);

        //write sites
        outputGLF->writeDummySites(300);
        outputGLF->closeOutput();
    }

    void read(){
        //open GLF for reading
        inputGLF = std::make_unique<TGlfReader>(_filename);

        // read!
        while(inputGLF->readNext()){
            positions.push_back(inputGLF->position());
            depths.push_back(inputGLF->depth());
        }
    }

    void SetUp() override{
        write();
        read();
    }

    void TearDown() override {};
};

TEST_F(TGLF_Test_ReadWrite, positions){
    // check if written and read positions are equal
    int c = 0;
    for (auto writtenPosition = outputGLF->beginPositions(); writtenPosition != outputGLF->endPositions(); writtenPosition++, c++){
        EXPECT_EQ(*writtenPosition, positions[c]);
    }
}

TEST_F(TGLF_Test_ReadWrite, depth){
    // check if written and read depth are equal
    int c = 0;
    for (auto writtenDepth = outputGLF->beginDepths(); writtenDepth != outputGLF->endDepths(); writtenDepth++, c++){
        EXPECT_EQ(*writtenDepth, (int) depths[c]);
    }
}

TEST_F(TGLF_Test_ReadWrite, chromosomes){
    // check if written and read chromosomes are equal

    // first chromosome
    TGlfChromosome * chr;
    inputGLF->fillPointerToChr(0, chr);

    EXPECT_EQ(chr->refId, 0);
    EXPECT_EQ(chr->name, "Chr1");
    EXPECT_EQ(chr->length, 100);
    EXPECT_EQ(chr->ploidy, 2);
    EXPECT_EQ(chr->numLikelihoodValues, 10); // diploid
    EXPECT_EQ(chr->maxNumLikelihoodValues, 10); // diploid

    // second chromosome
    inputGLF->fillPointerToChr(1, chr);

    EXPECT_EQ(chr->refId, 1);
    EXPECT_EQ(chr->name, "Chr2");
    EXPECT_EQ(chr->length, 200);
    EXPECT_EQ(chr->ploidy, 2);
    EXPECT_EQ(chr->numLikelihoodValues, 10); // diploid
    EXPECT_EQ(chr->maxNumLikelihoodValues, 10); // diploid

    // third chromosome
    inputGLF->fillPointerToChr(2, chr);

    EXPECT_EQ(chr->refId, 2);
    EXPECT_EQ(chr->name, "Chr3");
    EXPECT_EQ(chr->length, 300);
    EXPECT_EQ(chr->ploidy, 2);
    EXPECT_EQ(chr->numLikelihoodValues, 10); // diploid
    EXPECT_EQ(chr->maxNumLikelihoodValues, 10); // diploid
}





