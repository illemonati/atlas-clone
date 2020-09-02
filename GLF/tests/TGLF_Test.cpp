//
// Created by caduffm on 9/1/20.
//

#include "TGLF.h"
#include "TTestGLFFile.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

//-------------------------------------------------------------
// TGLF_Test_ReadWrite
// -- simple writing and reading
//-------------------------------------------------------------

class TGLF_Test_WriteRead : public ::testing::Test {
protected:
    std::string _filename = "testGLF.glf";
    std::vector<uint32_t> chrLength;
public:
    std::unique_ptr<TestUtilities::TTestGLFFile> outputGLF;
    std::unique_ptr<TGlfReader> inputGLF;

    // store stuff from input GLF to later compare with output
    std::vector<uint32_t> positions;
    std::vector<uint8_t> depths;
    std::vector<uint16_t *> genotypeLikelihoods;

    // convert
    TGlfConverter converter;

    virtual void write(){
        //settings
        chrLength = {100, 200, 300};

        //open GLF file for writing
        outputGLF = std::make_unique<TestUtilities::TTestGLFFile>(_filename, chrLength);

        //write sites
        outputGLF->writeDummySites(300);
        outputGLF->closeOutput();
    }

    virtual void read(){
        //open GLF for reading
        inputGLF = std::make_unique<TGlfReader>(_filename);

        // read!
        while(inputGLF->readNext()){
            // store
            positions.push_back(inputGLF->position());
            depths.push_back(inputGLF->depth());

            genotypeLikelihoods.push_back(new uint16_t [10]);
            inputGLF->fillGenotypeLikelihoodsGLF(genotypeLikelihoods.back());
        }
    }

    void SetUp() override{
        write();
        read();
    }

    void TearDown() override {
        for (auto & it : genotypeLikelihoods)
            delete [] it;
    };
};

TEST_F(TGLF_Test_WriteRead, name){
    EXPECT_EQ(inputGLF->name(), _filename);
}

TEST_F(TGLF_Test_WriteRead, positions){
    // check if written and read positions are equal
    int c = 0;
    for (auto writtenPosition = outputGLF->beginPositions(); writtenPosition != outputGLF->endPositions(); writtenPosition++, c++){
        EXPECT_EQ(writtenPosition->position(), positions[c]);
    }
}

TEST_F(TGLF_Test_WriteRead, depth){
    // check if written and read depth are equal
    int c = 0;
    for (auto writtenDepth = outputGLF->beginDepths(); writtenDepth != outputGLF->endDepths(); writtenDepth++, c++){
        EXPECT_EQ(*writtenDepth, (int) depths[c]);
    }
}

TEST_F(TGLF_Test_WriteRead, chromosomes){
    // check if written and read chromosomes are equal
    TGlfChromosome * chr;
    for (int i = 0; i < 3; i++){
        inputGLF->fillPointerToChr(i, chr);

        EXPECT_EQ(chr->refId, i);
        EXPECT_EQ(chr->name, "Chr" + toString(i + 1));
        EXPECT_EQ(chr->length, chrLength[i]);
        EXPECT_EQ(chr->ploidy, 2);
        EXPECT_EQ(chr->numLikelihoodValues, 10); // diploid
        EXPECT_EQ(chr->maxNumLikelihoodValues, 10); // diploid
    }
}

void normalizeByMax_Diploid(GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods){
    // find maxLL
    double maxLL = genotypeLikelihoods[0];
    for(int i=1; i < 10; ++i){
        if(genotypeLikelihoods[i] > maxLL)
            maxLL = genotypeLikelihoods[i];
    }

    // normalize
    for(int i=0; i < 10; ++i){
        genotypeLikelihoods[i] /= maxLL;
    }
}

void normalizeByMax_Haploid(GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods){
    // find maxLL
    double maxLL = genotypeLikelihoods[AA];
    if(genotypeLikelihoods[CC] > maxLL) maxLL = genotypeLikelihoods[CC];
    if(genotypeLikelihoods[GG] > maxLL) maxLL = genotypeLikelihoods[GG];
    if(genotypeLikelihoods[TT] > maxLL) maxLL = genotypeLikelihoods[TT];

    // normalize
    genotypeLikelihoods[AA] /= maxLL;
    genotypeLikelihoods[CC] /= maxLL;
    genotypeLikelihoods[GG] /= maxLL;
    genotypeLikelihoods[TT] /= maxLL;

}

TEST_F(TGLF_Test_WriteRead, genotypeLikelihoods){
    // check if written and read genotype likelihoods are equal
    int c = 0;
    for (auto writtenGTL = outputGLF->beginGenotypeLikelihoods(); writtenGTL != outputGLF->endGenotypeLikelihoods(); writtenGTL++, c++){
        // need to normalize the written likelihoods by maximal LL in order to compare
        normalizeByMax_Diploid(*writtenGTL);
        for (int g = 0; g < 10; g++){ // go over all 10 possible genotypes
            // compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood -> likelihood)
            EXPECT_EQ(converter.toGlfFormat((*writtenGTL)[g]), genotypeLikelihoods[c][g]);
        }
    }
}

TEST_F(TGLF_Test_WriteRead, rewind){
    positions.clear();
    inputGLF->rewind();
    read();

    // check if written and read positions are equal
    int c = 0;
    for (auto writtenPosition = outputGLF->beginPositions(); writtenPosition != outputGLF->endPositions(); writtenPosition++, c++){
        EXPECT_EQ(writtenPosition->position(), positions[c]);
    }
}

TEST_F(TGLF_Test_WriteRead, jumpToNextChr){
    inputGLF->rewind();

    // we start on chromosome 1
    EXPECT_EQ(0, inputGLF->refId()); EXPECT_EQ(100, inputGLF->chrLength());

    // jump to first site on chromosome 2
    EXPECT_TRUE(inputGLF->jumpToNextChr());
    EXPECT_EQ(1, inputGLF->refId()); EXPECT_EQ(200, inputGLF->chrLength());
    uint32_t whichSite = 0;
    for (auto pos = outputGLF->beginPositions(); pos != outputGLF->endPositions(); pos++, whichSite++){
        if (pos->refID() == 1) break;
    }
    EXPECT_EQ(outputGLF->position(whichSite).position(), inputGLF->position()); EXPECT_EQ(outputGLF->depth(whichSite), inputGLF->depth());

    // jump to next chromosome (=3) -> skip all other sites of chromosome 2
    EXPECT_TRUE(inputGLF->jumpToNextChr());
    EXPECT_EQ(2, inputGLF->refId()); EXPECT_EQ(300, inputGLF->chrLength());
    whichSite = 0;
    for (auto pos = outputGLF->beginPositions(); pos != outputGLF->endPositions(); pos++, whichSite++){
        if (pos->refID() == 2) break;
    }
    EXPECT_EQ(outputGLF->position(whichSite).position(), inputGLF->position()); EXPECT_EQ(outputGLF->depth(whichSite), inputGLF->depth());

    // chromosome 3 was the last one -> jumping to the next one returns false
    EXPECT_FALSE(inputGLF->jumpToNextChr());
}

//-------------------------------------------------------------
// TGLF_Test_WriteRead_Missing
// -- writing and reading with missing data
//-------------------------------------------------------------

class TGLF_Test_WriteRead_Missing : public TGLF_Test_WriteRead {
public:
    void write() override{
        // write some more complicated stuff: missing sites, depth = 0, mapping qual = 0, missing chromosomes etc.

        //settings
        chrLength = {50, 100, 150, 200, 250};

        //open GLF file for writing
        outputGLF = std::make_unique<TestUtilities::TTestGLFFile>(_filename, chrLength);

        //write sites
        // 1) first chromosome is empty
        outputGLF->writeNewChromosome();

        // second chromosome has sites, but some sites have depth 0 or base N
        outputGLF->writeNewChromosome();
        // 2) ok site
        outputGLF->writeDummySite(10, 10);
        // 3) depth = 0
        outputGLF->writeDummySite(20, 0);
        // 4) depth = 10, but all bases are N
        std::vector<GenotypeLikelihoods::TBaseData> bases;
        bases.reserve(10);
        Base base = N;
        double error = 0.001;
        for (uint32_t d = 0; d < 10; d++){
            TBaseData baseData(base, error);
            bases.emplace_back(baseData);
        }
        TGenotypeLikelihoods gtL;
        gtL.fill(bases);
        outputGLF->writeDummySite(30, 0, gtL);
        // 5) depth = 10, all bases are A, but mapping quality is zero
        bases.clear();
        base = A;
        for (uint32_t d = 0; d < 10; d++){
            TBaseData baseData(base, error);
            bases.emplace_back(baseData);
        }
        gtL.fill(bases);
        outputGLF->writeSite(40, 0, gtL, 0);

        // 6) third chromosome is empty
        outputGLF->writeNewChromosome();

        // 7) fourth chromosome has sites
        outputGLF->writeNewChromosome();
        // only last site of chromosome is not missing
        outputGLF->writeDummySite(200, 10);

        // 8) last chromosome is empty
        outputGLF->writeNewChromosome();

        outputGLF->closeOutput();
    }
};

TEST_F(TGLF_Test_WriteRead_Missing, positions){
    // check if written and read positions are equal
    int c = 0;
    for (auto writtenPosition = outputGLF->beginPositions(); writtenPosition != outputGLF->endPositions(); writtenPosition++, c++){
        EXPECT_EQ(writtenPosition->position(), positions[c]);
    }
}

TEST_F(TGLF_Test_WriteRead_Missing, depth){
    // check if written and read depth are equal
    int c = 0;
    for (auto writtenDepth = outputGLF->beginDepths(); writtenDepth != outputGLF->endDepths(); writtenDepth++, c++){
        EXPECT_EQ(*writtenDepth, (int) depths[c]);
    }
}

TEST_F(TGLF_Test_WriteRead_Missing, chromosomes){
    // check if written and read chromosomes are equal
    // first chromosome
    TGlfChromosome * chr;
    for (int i = 0; i < 5; i++){
        inputGLF->fillPointerToChr(i, chr);

        EXPECT_EQ(chr->refId, i);
        EXPECT_EQ(chr->name, "Chr" + toString(i + 1));
        EXPECT_EQ(chr->length, chrLength[i]);
        EXPECT_EQ(chr->ploidy, 2);
        EXPECT_EQ(chr->numLikelihoodValues, 10); // diploid
        EXPECT_EQ(chr->maxNumLikelihoodValues, 10); // diploid
    }
}

TEST_F(TGLF_Test_WriteRead_Missing, genotypeLikelihoods){
    // check if written and read genotype likelihoods are equal
    int c = 0;
    for (auto writtenGTL = outputGLF->beginGenotypeLikelihoods(); writtenGTL != outputGLF->endGenotypeLikelihoods(); writtenGTL++, c++){
        // need to normalize the written likelihoods by maximal LL in order to compare
        normalizeByMax_Diploid(*writtenGTL);
        for (int g = 0; g < 10; g++){ // go over all 10 possible genotypes
            // compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood -> likelihood)
            EXPECT_EQ(converter.toGlfFormat((*writtenGTL)[g]), genotypeLikelihoods[c][g]);
        }
    }
}

TEST_F(TGLF_Test_WriteRead_Missing, jumpToNextChr){
    positions.clear(); depths.clear();
    inputGLF->rewind();

    // we start on chromosome 1
    EXPECT_EQ(0, inputGLF->refId()); EXPECT_EQ(50, inputGLF->chrLength());

    // jump to first site on chromosome 2
    inputGLF->jumpToNextChr();
    EXPECT_EQ(10, inputGLF->position()); EXPECT_EQ(10, inputGLF->depth());
    EXPECT_EQ(1, inputGLF->refId()); EXPECT_EQ(100, inputGLF->chrLength());

    // jump to next chromosome (=3), but this is empty, so jump to chromosome 4 and reads first site on chr4 -> skip all other sites of chromosome 2
    inputGLF->jumpToNextChr();
    EXPECT_EQ(3, inputGLF->refId()); EXPECT_EQ(200, inputGLF->chrLength());
    EXPECT_EQ(200, inputGLF->position()); EXPECT_EQ(10, inputGLF->depth());

    // jump to chromosome 5 -> this is the last one -> returns false
    EXPECT_FALSE(inputGLF->jumpToNextChr());
    EXPECT_EQ(4, inputGLF->refId()); EXPECT_EQ(250, inputGLF->chrLength());
}



//-------------------------------------------------------------
// TGLF_Test_WriteRead_Ploidies
// -- writing and reading with chromosomes of different ploidies
//-------------------------------------------------------------

class TGLF_Test_WriteRead_Ploidies : public TGLF_Test_WriteRead {
public:
    void write() override{
        // write chromosomes with different ploidies
        chrLength = {50, 100, 150};
        std::vector<uint8_t> ploidies = {2, 1, 2}; // second chromosome is haploid

        //open GLF file for writing
        outputGLF = std::make_unique<TestUtilities::TTestGLFFile>(_filename, chrLength, ploidies);

        //write sites
        // 1) first chromosome is diploid
        outputGLF->writeNewChromosome();
        outputGLF->writeDummySite(10);
        outputGLF->writeDummySite(20);

        // 2) second chromosome is haploid
        outputGLF->writeNewChromosome();
        outputGLF->writeDummySite(10);
        outputGLF->writeDummySite(20);

        // 3) third chromosome is diploid
        outputGLF->writeNewChromosome();
        outputGLF->writeDummySite(10);
        outputGLF->writeDummySite(20);

        outputGLF->closeOutput();
    }

    void read() override{
        //open GLF for reading
        inputGLF = std::make_unique<TGlfReader>(_filename);

        // read!
        while(inputGLF->readNext()){
            // store
            positions.push_back(inputGLF->position());
            depths.push_back(inputGLF->depth());

            if (inputGLF->chrIsHaploid()) genotypeLikelihoods.push_back(new uint16_t [4]);
            else genotypeLikelihoods.push_back(new uint16_t [10]);

            inputGLF->fillGenotypeLikelihoodsGLF(genotypeLikelihoods.back());
        }
    }
};

TEST_F(TGLF_Test_WriteRead_Ploidies, positions){
    // check if written and read positions are equal
    int c = 0;
    for (auto writtenPosition = outputGLF->beginPositions(); writtenPosition != outputGLF->endPositions(); writtenPosition++, c++){
        EXPECT_EQ(writtenPosition->position(), positions[c]);
    }
}

TEST_F(TGLF_Test_WriteRead_Ploidies, depth){
    // check if written and read depth are equal
    int c = 0;
    for (auto writtenDepth = outputGLF->beginDepths(); writtenDepth != outputGLF->endDepths(); writtenDepth++, c++){
        EXPECT_EQ(*writtenDepth, (int) depths[c]);
    }
}

TEST_F(TGLF_Test_WriteRead_Ploidies, chromosomes){
    // check if written and read chromosomes are equal
    // first chromosome
    TGlfChromosome * chr;
    for (int i = 0; i < 3; i++){
        inputGLF->fillPointerToChr(i, chr);

        EXPECT_EQ(chr->refId, i);
        EXPECT_EQ(chr->name, "Chr" + toString(i + 1));
        EXPECT_EQ(chr->length, chrLength[i]);
        EXPECT_EQ(chr->maxNumLikelihoodValues, 10);
        if (i == 1){ // haploid
            EXPECT_EQ(chr->ploidy, 1);
            EXPECT_EQ(chr->numLikelihoodValues, 4);
        }
        else { // diploid
            EXPECT_EQ(chr->ploidy, 2);
            EXPECT_EQ(chr->numLikelihoodValues, 10);
        }
    }
}

TEST_F(TGLF_Test_WriteRead_Ploidies, genotypeLikelihoods){
    // check if written and read genotype likelihoods are equal
    int c = 0;
    for (auto writtenGTL = outputGLF->beginGenotypeLikelihoods(); writtenGTL != outputGLF->endGenotypeLikelihoods(); writtenGTL++, c++){
        if (c == 2 || c == 3){
            normalizeByMax_Haploid(*writtenGTL);
            // haploid -> compare homozygous genotypes! ATTENTION: if haploid, only 4 values are returned -> access with ACGT, not AA CC GG TT!!
            EXPECT_EQ(converter.toGlfFormat((*writtenGTL)[AA]), genotypeLikelihoods[c][A]);
            EXPECT_EQ(converter.toGlfFormat((*writtenGTL)[GG]), genotypeLikelihoods[c][G]);
            EXPECT_EQ(converter.toGlfFormat((*writtenGTL)[TT]), genotypeLikelihoods[c][T]);
            EXPECT_EQ(converter.toGlfFormat((*writtenGTL)[CC]), genotypeLikelihoods[c][C]);
        } else {
            // need to normalize the written likelihoods by maximal LL in order to compare
            normalizeByMax_Diploid(*writtenGTL);
            for (int g = 0; g < 10; g++) { // go over all possible genotypes
                // compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood -> likelihood)
                EXPECT_EQ(converter.toGlfFormat((*writtenGTL)[g]), genotypeLikelihoods[c][g]);
            }
        }
    }
}

//-------------------------------------------------------------
// TGLF_Test_WriteRead_Windows
// -- read in windows
//-------------------------------------------------------------

class TGLF_Test_WriteRead_Windows : public TGLF_Test_WriteRead {
protected:
    uint32_t windowLen = 20;
    std::vector< std::vector<uint16_t *> > genotypeLikelihoods_perWindow;

public:
    void read() override{
        //open GLF for reading
        inputGLF = std::make_unique<TGlfReader>(_filename);

        //initialize variables
        bool isGood = true;

        //prepare variables
        uint32_t curRefId;
        long curChrLen;
        long windowStart;
        long windowEnd;
        //parse GLFs in windows
        while(!inputGLF->eof()){
            //move to new chromosome
            curRefId = inputGLF->refId();
            curChrLen = inputGLF->chrLength();
            windowStart = 1;
            windowEnd = windowLen + 1;

            std::cout << "NEW CHROMOSOME!!" << std::endl;
            std::cout << "\t -> curRefId = " <<  curRefId << std::endl;
            std::cout << "\t -> curChrLen = " <<  curChrLen << std::endl;

            //parse all windows of chromosome
            while(windowStart < curChrLen) {
                std::cout << "\t NEW WINDOW!!" << std::endl;
                std::cout << "\t\t -> windowStart = " <<  windowStart << std::endl;
                std::cout << "\t\t -> windowEnd = " <<  windowEnd << std::endl;
                // resize storage
                std::vector<uint16_t*> genoLikelihoods_oneWindow(windowLen);
                for(int i=0; i<windowLen; ++i){
                    genoLikelihoods_oneWindow[i] = new uint16_t[10];
                }

                //read data
                isGood = inputGLF->readNextWindow(genoLikelihoods_oneWindow, curRefId, windowStart, windowEnd);
                if (isGood || inputGLF->eof()) {
                    // store
                    genotypeLikelihoods_perWindow.push_back(genoLikelihoods_oneWindow);
                    for (auto & it : genoLikelihoods_oneWindow){
                        for (int i = 0; i < 10; i++)
                            std::cout << (int) it[i] << "  ";
                        std::cout << std::endl;
                    }
                }
                //move window
                windowStart = windowEnd;
                windowEnd = windowStart + windowLen;
            }
        }
    }

    void TearDown() override {
        //clean up memory
        for (auto & it : genotypeLikelihoods_perWindow) {
            for(int i=0; i<windowLen; ++i) {
                delete[] it[i];
            }
        }
    }
};

/*
TEST_F(TGLF_Test_WriteRead_Windows, genotypeLikelihoods){
    // check if written and read genotype likelihoods are equal
    auto writtenGTL = outputGLF->beginGenotypeLikelihoodsWithMissingSites();

    for (; writtenGTL != outputGLF->endGenotypeLikelihoodsWithMissingSites(); writtenGTL++){
        normalizeByMax_Diploid(*writtenGTL);
        for (int g = 0; g < 10; g++) { // go over all 10 possible genotypes
            std::cout <<  converter.toGlfFormat((*writtenGTL)[g]) << "  ";
        }
        std::cout << std::endl;
    }
}
*/

TEST_F(TGLF_Test_WriteRead_Windows, genotypeLikelihoods){
    // check if written and read genotype likelihoods are equal
    auto writtenGTL = outputGLF->beginGenotypeLikelihoodsWithMissingSites();

    for (auto & window : genotypeLikelihoods_perWindow){
        for (auto genotypeLikelihood_read : window){

            // need to normalize the written likelihoods by maximal LL in order to compare
            normalizeByMax_Diploid(*writtenGTL);

            for (int g = 0; g < 10; g++){ // go over all 10 possible genotypes
                // compare in GLF format (quite a large imprecision when going from likelihood -> GLF likelihood -> likelihood)
                EXPECT_EQ(converter.toGlfFormat((*writtenGTL)[g]), genotypeLikelihood_read[g]);
                std::cout << converter.toGlfFormat((*writtenGTL)[g]) << " <--> " << genotypeLikelihood_read[g] << std::endl;

            }
            writtenGTL++;
        }
    }
    EXPECT_EQ(writtenGTL, outputGLF->endGenotypeLikelihoodsWithMissingSites());
}