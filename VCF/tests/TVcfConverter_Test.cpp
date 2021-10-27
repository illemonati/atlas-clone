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

    void writeVcfFile(const std::vector<std::string> & SampleNames) {
        gz::ogzstream vcf;
        vcf.open(filename.c_str());
        _writeHeader(vcf, SampleNames);

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
};

TEST_F(TVCFConverterTest, dummy){
    std::vector<std::string> sampleNames = {"Indiv1", "Indiv2", "Indiv3", "Indiv4", "Indiv5"};
    writeVcfFile(sampleNames);

    // convert to beagle
    VCF::TVcfToBeagle converter(&logfile);
    converter.vcfToBeagle(parameters);

    // read beagle using stattools reader
    stattools::TDataReaderBase<coretools::WeakType<double>> reader;
    auto observation = std::make_shared<stattools::TObservationTypedDefault<coretools::WeakType<double>>>("observation", &randomGenerator);
    std::shared_ptr<stattools::TNamesEmpty> individualNames = std::make_shared<stattools::TNamesStrings>();
    std::shared_ptr<stattools::TNamesEmpty> lociNames = std::make_shared<stattools::TNamesStrings>();
    reader.read("test.beagle.gz", 2, observation, {lociNames, individualNames}, {'\n', '\t'}, {true, true}, {numLoci,numIndiv}, {true, true}, numLoci, stattools::colNames_multiline);

    // check if genotype likelihoods are as expected
    size_t linearIndex = 0;
    for (size_t l = 0; l < numLoci; ++l) { // loop over loci
        for (size_t i = 0; i < numIndiv; ++i, ++linearIndex) { // loop over individuals
            // genotype likelihoods
            EXPECT_DOUBLE_EQ(observation->value(linearIndex), (coretools::Probability) genometools::PhredIntProbability(phred_g1[linearIndex]));
            EXPECT_DOUBLE_EQ(observation->value(linearIndex+1), (coretools::Probability) genometools::PhredIntProbability(phred_g2[linearIndex]));
            EXPECT_DOUBLE_EQ(observation->value(linearIndex+2), (coretools::Probability) genometools::PhredIntProbability(phred_g3[linearIndex]));
        }
    }
}

