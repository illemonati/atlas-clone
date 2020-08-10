
#include "TBamFile.h"
#include "TTestBamFile.h"
#include "TAlignment.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

//-------------------------------------------------------------
// TBamFile
//-------------------------------------------------------------

TEST(TBamFile, writeRead){
	//settings
	std::string filename = "testBAM.bam";
	std::vector<uint32_t> chrLength = {1000, 2000, 3000};
	uint32_t numReadGroups = 2;

	//open BAM file for writing
    TestUtilities::TTestBamFile testBam(filename, chrLength, numReadGroups);

    //write alignments
    testBam.writeDummyAlignments(10);


    //EXPECT_EQ(chr.refId, 1);
    //EXPECT_EQ(chr.name, "1");
}

