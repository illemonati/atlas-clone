//
// Created by caduffm on 10/27/21.
//

#include "TestCase.h"
#include "../TVcfConverter.h"
#include "TVcfFile.h"
#include "TDataFile.h"
#include "TObservationTypedDefault.h"

using namespace testing;

class TVCFConverterTest : public Test {
protected:
    std::vector<uint16_t> phred_g1 = {0,1,0,0,3,0,1,2,0,0,2,0,1,0,0,1,0,0,1,1,1,0,0,3,1,0,0,1,0,1,1,0,0,1,1,1,2,3,2,0,3,0,0,1,0,0,0,1,2,0};
    std::vector<uint16_t> phred_g2 = {1,2,0,0,0,1,1,0,1,1,1,0,0,2,1,2,2,0,1,1,0,2,0,0,0,0,0,0,2,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,1,0,0,1};
    std::vector<uint16_t> phred_g3 = {1,0,1,3,1,2,0,1,0,1,0,0,1,0,2,0,0,0,0,0,0,0,1,0,0,1,3,0,0,1,0,1,0,0,1,1,0,1,0,1,2,3,0,0,1,0,2,1,0,1};
    std::vector<uint16_t> depth = {14,11,9,9,6,8,12,9,10,11,5,14,7,13,13,11,11,5,9,8,9,8,6,7,10,14,10,11,13,10,7,9,11,13,8,7,9,6,10,6,15,15,10,9,11,9,6,14,9,14};
    std::vector<std::string> sampleNames = {"Indiv1", "Indiv2", "Indiv3", "Indiv4", "Indiv5"};
    std::vector<std::string> samplesToKeep = {"Indiv4", "Indiv1", "Indiv5"}; // omit Indiv3 and Indiv4; and shuffle
    std::vector<size_t> indexInSampleNames = {3, 0, 4};

    void _writeHeader(gz::ogzstream & VCF, const std::vector<std::string> & SampleNames) {
        //write info
        VCF << "##fileformat=VCFv4.2\n";
        VCF << "##source=Simulation\n";
        VCF << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        VCF << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
        VCF << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n";
        VCF << "##FORMAT=<ID=DP,Number=G,Type=Integer,Description=\"Depth at site\">\n";

        // write header
        VCF << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (const auto & SampleName : SampleNames){
            VCF << "\t" << SampleName;
        }
    }

    void _writeCell(gz::ogzstream & VCF, const genometools::TSampleLikelihoods<genometools::PhredIntProbability> & SampleLikelihoods, const size_t & Depth){
        // genotype with lowest phred-score is the observed genotype
        genometools::BiallelicGenotype observedGenotype = SampleLikelihoods.mostLikelyGenotype();
        genometools::BiallelicGenotype secondBestGenotype = SampleLikelihoods.secondMostLikelyGenotype();

        // write to vcf
        VCF << "\t" << (std::string) observedGenotype << ":" << SampleLikelihoods[secondBestGenotype] << ":";
        VCF << coretools::str::toString(SampleLikelihoods[genometools::BiallelicGenotype(genometools::homoFirst)])
                + "," + coretools::str::toString(SampleLikelihoods[genometools::BiallelicGenotype(genometools::het)])
                + "," + coretools::str::toString(SampleLikelihoods[genometools::BiallelicGenotype(genometools::homoSecond)]);
        VCF << ":" << Depth;
    };

public:
    coretools::TLog logfile;
    coretools::TParameters parameters;
    coretools::TRandomGenerator randomGenerator;
    std::string filename;

    uint32_t numLoci;
    uint32_t numIndiv;

    void SetUp(){
        numLoci = 10;
        numIndiv = 5;
        filename = "test.vcf.gz";
        parameters.addParameter("vcf", filename);
        parameters.addParameter("minMAF", "0");
    }

    void writeVcfFile() {
        gz::ogzstream vcf;
        vcf.open(filename.c_str());
        _writeHeader(vcf, sampleNames);

        std::string chr = "junk1";
        size_t linearIndex = 0;
        for (size_t l = 0; l < numLoci; ++l) { // loop over loci
            vcf << '\n' << chr << '\t' << l+1 << "\t.\tA\tC\t50\t.\t.\tGT:GQ:PL:DP";
            for (size_t i = 0; i < numIndiv; ++i, ++linearIndex) { // loop over individuals
                // write to vcf
                auto GTL0 = genometools::PhredIntProbability(phred_g1[linearIndex]);
                auto GTL1 = genometools::PhredIntProbability(phred_g2[linearIndex]);
                auto GTL2 = genometools::PhredIntProbability(phred_g3[linearIndex]);
                _writeCell(vcf, {GTL0, GTL1, GTL2}, depth[linearIndex]);
            }
        }
        vcf << "\n";
        vcf.close();
    }

    void writeSampleFile(){
        coretools::TOutputFile sampleFile("test.samples");
        for (auto & it : samplesToKeep){
            sampleFile << it << std::endl;
        }
        // add to parameters
        parameters.addParameter("samples", "test.samples");
    }
};

TEST_F(TVCFConverterTest, beagle){
    writeVcfFile();

    // convert to beagle
    VCF::TVcfToBeagle converter(&logfile);
    converter.vcfToBeagle(parameters);

    // read beagle
    coretools::TInputFile beagle("test.beagle.gz", coretools::TFile_Filetype::header);
    EXPECT_EQ(beagle.headerAt(0), "marker");
    EXPECT_EQ(beagle.headerAt(1), "allele1");
    EXPECT_EQ(beagle.headerAt(2), "allele2");
    for (size_t i = 0; i < numIndiv; ++i){
        EXPECT_EQ(beagle.headerAt(3+i*3), sampleNames[i]);
        EXPECT_EQ(beagle.headerAt(3+i*3+1), sampleNames[i]);
        EXPECT_EQ(beagle.headerAt(3+i*3+2), sampleNames[i]);
    }

    // check if genotype likelihoods are as expected
    std::vector<std::string> line;
    size_t linearIndex = 0;
    for (size_t l = 0; l < numLoci; ++l) { // loop over loci
        beagle.read(line);
        EXPECT_EQ(line[0], "junk1_"+coretools::str::toString(l+1));
        EXPECT_EQ(line[1], "A");
        EXPECT_EQ(line[2], "C");
        for (size_t i = 0; i < numIndiv; ++i, ++linearIndex) { // loop over individuals
            // genotype likelihoods (we loose some precision when reading/writing, so only expect near)
            EXPECT_NEAR(coretools::str::convertString<coretools::Probability>(line[3+i*3]), (coretools::Probability) genometools::PhredIntProbability(phred_g1[linearIndex]), 0.00001);
            EXPECT_NEAR(coretools::str::convertString<coretools::Probability>(line[3+i*3+1]), (coretools::Probability) genometools::PhredIntProbability(phred_g2[linearIndex]), 0.00001);
            EXPECT_NEAR(coretools::str::convertString<coretools::Probability>(line[3+i*3+2]), (coretools::Probability) genometools::PhredIntProbability(phred_g3[linearIndex]), 0.00001);
        }
    }
}

TEST_F(TVCFConverterTest, beagle_withSamples){
    writeVcfFile();
    writeSampleFile();

    // convert to beagle
    VCF::TVcfToBeagle converter(&logfile);
    converter.vcfToBeagle(parameters);

    // read beagle
    coretools::TInputFile beagle("test.beagle.gz", coretools::TFile_Filetype::header);
    EXPECT_EQ(beagle.headerAt(0), "marker");
    EXPECT_EQ(beagle.headerAt(1), "allele1");
    EXPECT_EQ(beagle.headerAt(2), "allele2");
    for (size_t i = 0; i < samplesToKeep.size(); ++i){
        EXPECT_EQ(beagle.headerAt(3+i*3), samplesToKeep[i]);
        EXPECT_EQ(beagle.headerAt(3+i*3+1), samplesToKeep[i]);
        EXPECT_EQ(beagle.headerAt(3+i*3+2), samplesToKeep[i]);
    }

    // check if genotype likelihoods are as expected
    std::vector<std::string> line;
    for (size_t l = 0; l < numLoci; ++l) { // loop over loci
        beagle.read(line);
        EXPECT_EQ(line[0], "junk1_"+coretools::str::toString(l+1));
        EXPECT_EQ(line[1], "A");
        EXPECT_EQ(line[2], "C");
        for (size_t i = 0; i < samplesToKeep.size(); ++i) { // loop over individuals
            size_t relevantIndex = l*numIndiv + indexInSampleNames[i];
            // genotype likelihoods (we loose some precision when reading/writing, so only expect near)
            EXPECT_NEAR(coretools::str::convertString<coretools::Probability>(line[3+i*3]), (coretools::Probability) genometools::PhredIntProbability(phred_g1[relevantIndex]), 0.00001);
            EXPECT_NEAR(coretools::str::convertString<coretools::Probability>(line[3+i*3+1]), (coretools::Probability) genometools::PhredIntProbability(phred_g2[relevantIndex]), 0.00001);
            EXPECT_NEAR(coretools::str::convertString<coretools::Probability>(line[3+i*3+2]), (coretools::Probability) genometools::PhredIntProbability(phred_g3[relevantIndex]), 0.00001);
        }
    }
}

TEST_F(TVCFConverterTest, geno){
    writeVcfFile();

    // convert to geno
    VCF::TVcfToGeno converter(&logfile);
    converter.vcfToGeno(parameters);

    // read geno
    coretools::TInputFile geno("test.geno", coretools::TFile_Filetype::fixed);

    // check if genotypes are as expected
    std::vector<std::string> line;
    for (size_t l = 0; l < numLoci; ++l) { // loop over loci
        geno.read(line);
        std::vector<char> genotypes(line[0].begin(), line[0].end());
        for (size_t i = 0; i < numIndiv; ++i) { // loop over individuals
            // genotypes
            size_t relevantIndex = l*numIndiv + i;
            auto GTL0 = genometools::PhredIntProbability(phred_g1[relevantIndex]);
            auto GTL1 = genometools::PhredIntProbability(phred_g2[relevantIndex]);
            auto GTL2 = genometools::PhredIntProbability(phred_g3[relevantIndex]);
            genometools::TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
            genometools::BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

            EXPECT_EQ(coretools::str::convertString<uint8_t>(std::string(1, genotypes[i])), 2 - (uint8_t) observedGenotype); // genotype defined based on reference allele
        }
    }
}

TEST_F(TVCFConverterTest, geno_withSamples){
    writeVcfFile();
    writeSampleFile();

    // convert to geno
    VCF::TVcfToGeno converter(&logfile);
    converter.vcfToGeno(parameters);

    // read geno
    coretools::TInputFile geno("test.geno", coretools::TFile_Filetype::fixed);

    // check if genotypes are as expected
    std::vector<std::string> line;
    for (size_t l = 0; l < numLoci; ++l) { // loop over loci
        geno.read(line);
        std::vector<char> genotypes(line[0].begin(), line[0].end());
        for (size_t i = 0; i < samplesToKeep.size(); ++i) { // loop over individuals
            // genotypes
            size_t relevantIndex = l*numIndiv + indexInSampleNames[i];
            auto GTL0 = genometools::PhredIntProbability(phred_g1[relevantIndex]);
            auto GTL1 = genometools::PhredIntProbability(phred_g2[relevantIndex]);
            auto GTL2 = genometools::PhredIntProbability(phred_g3[relevantIndex]);
            genometools::TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
            genometools::BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

            EXPECT_EQ(coretools::str::convertString<uint8_t>(std::string(1, genotypes[i])), 2 - (uint8_t) observedGenotype); // genotype defined based on reference allele
        }
    }
}

TEST_F(TVCFConverterTest, lfmmCalledGeno){
    writeVcfFile();

    // convert to lfmm
    VCF::TVcfToLFMMCalledGeno converter(&logfile);
    converter.vcfToLFMM(parameters);

    // read lfmm
    coretools::TInputFile lfmm("test.lfmm", coretools::TFile_Filetype::fixed);

    // check if genotypes are as expected
    std::vector<std::string> line;
    for (size_t i = 0; i < numIndiv; ++i) { // loop over individuals
        lfmm.read(line);
        for (size_t l = 0; l < numLoci; ++l) { // loop over loci
            // genotypes
            size_t relevantIndex = l*numIndiv + i;
            auto GTL0 = genometools::PhredIntProbability(phred_g1[relevantIndex]);
            auto GTL1 = genometools::PhredIntProbability(phred_g2[relevantIndex]);
            auto GTL2 = genometools::PhredIntProbability(phred_g3[relevantIndex]);
            genometools::TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
            genometools::BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

            EXPECT_EQ(coretools::str::convertString<uint8_t>(line[l]), (uint8_t) observedGenotype);
        }
    }
}

TEST_F(TVCFConverterTest, lfmmCalledGeno_withSamples){
    writeVcfFile();
    writeSampleFile();

    // convert to lfmm
    VCF::TVcfToLFMMCalledGeno converter(&logfile);
    converter.vcfToLFMM(parameters);

    // read lfmm
    coretools::TInputFile lfmm("test.lfmm", coretools::TFile_Filetype::fixed);

    // check if genotypes are as expected
    std::vector<std::string> line;
    for (size_t i = 0; i < samplesToKeep.size(); ++i) { // loop over individuals
        lfmm.read(line);
        for (size_t l = 0; l < numLoci; ++l) { // loop over loci
            // genotypes
            size_t relevantIndex = l*numIndiv + indexInSampleNames[i];
            auto GTL0 = genometools::PhredIntProbability(phred_g1[relevantIndex]);
            auto GTL1 = genometools::PhredIntProbability(phred_g2[relevantIndex]);
            auto GTL2 = genometools::PhredIntProbability(phred_g3[relevantIndex]);
            genometools::TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
            genometools::BiallelicGenotype observedGenotype = sampleLikelihoods.mostLikelyGenotype();

            EXPECT_EQ(coretools::str::convertString<uint8_t>(line[l]), (uint8_t) observedGenotype);
        }
    }
}

TEST_F(TVCFConverterTest, lfmmMeanPosteriorGeno){
    writeVcfFile();

    // convert to lfmm
    VCF::TVcfToLFMMPostGeno converter(&logfile);
    converter.vcfToLFMM(parameters);

    // read lfmm
    coretools::TInputFile lfmm("test.lfmm", coretools::TFile_Filetype::fixed);

    // check if genotype likelihoods are as expected
    std::vector<std::string> line;
    for (size_t i = 0; i < numIndiv; ++i) { // loop over individuals
        lfmm.read(line);
        for (size_t l = 0; l < numLoci; ++l) { // loop over loci
            // genotypes
            size_t relevantIndex = l*numIndiv + i;
            auto GTL0 = genometools::PhredIntProbability(phred_g1[relevantIndex]);
            auto GTL1 = genometools::PhredIntProbability(phred_g2[relevantIndex]);
            auto GTL2 = genometools::PhredIntProbability(phred_g3[relevantIndex]);
            genometools::TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
            double posteriorGenotype = sampleLikelihoods.meanPosteriorGenotype();

            EXPECT_NEAR(coretools::str::convertString<double>(line[l]), posteriorGenotype, 0.00001);
        }
    }
}

TEST_F(TVCFConverterTest, lfmmMeanPosteriorGeno_withSamples){
    writeVcfFile();
    writeSampleFile();

    // convert to lfmm
    VCF::TVcfToLFMMPostGeno converter(&logfile);
    converter.vcfToLFMM(parameters);

    // read lfmm
    coretools::TInputFile lfmm("test.lfmm", coretools::TFile_Filetype::fixed);

    // check if genotypes are as expected
    std::vector<std::string> line;
    for (size_t i = 0; i < samplesToKeep.size(); ++i) { // loop over individuals
        lfmm.read(line);
        for (size_t l = 0; l < numLoci; ++l) { // loop over loci
            // genotypes
            size_t relevantIndex = l*numIndiv + indexInSampleNames[i];
            auto GTL0 = genometools::PhredIntProbability(phred_g1[relevantIndex]);
            auto GTL1 = genometools::PhredIntProbability(phred_g2[relevantIndex]);
            auto GTL2 = genometools::PhredIntProbability(phred_g3[relevantIndex]);
            genometools::TSampleLikelihoods sampleLikelihoods(GTL0, GTL1, GTL2);
            double posteriorGenotype = sampleLikelihoods.meanPosteriorGenotype();

            EXPECT_NEAR(coretools::str::convertString<double>(line[l]), posteriorGenotype, 0.00001);
        }
    }
}

