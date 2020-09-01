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

    void write(){
        //settings
        std::vector<uint32_t> chrLength = {1000, 2000, 3000};

        //open GLF file for writing
        outputGLF = std::make_unique<TestUtilities::TTestGLFFile>(_filename, chrLength);

        //write sites
        outputGLF->writeDummySites(100);
        outputGLF->closeOutput();
    }

    void read(){
        //open GLF for reading
        inputGLF = std::make_unique<TGlfReader>(_filename);
    }

    void SetUp() override{
        write();
        read();
    }

    void TearDown() override {};
};

TEST_F(TGLF_Test_ReadWrite, init){
    EXPECT_EQ(2, 2);
}





