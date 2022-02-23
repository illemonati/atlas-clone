//
// Created by reynac on 2/22/22.
//

#include "TF2Estimator.h"
namespace PopulationTools {

    using coretools::TParameters;
    using coretools::TLog;
    using coretools::TRandomGenerator;
    using coretools::TOutputFile;
    using coretools::TFactorial;

//------------------------------------------------
//TF2Estimator
//------------------------------------------------
    TF2Estimator::TF2Estimator(TParameters &Parameters, TLog *Logfile, TRandomGenerator *RandonGenerator) {
        _logfile = Logfile;
        _randonGenerator = RandonGenerator;

        //read list of provided samples
        if (Parameters.parameterExists("samples")) {
            _samples.readSamples(Parameters.getParameter<std::string>("samples"), _logfile);
        }

        //open VCF
        _vcfFilename = Parameters.getParameter<std::string>("vcf");
        _logfile->list("Reading VCF from file '" + _vcfFilename + "'.");
        _vcfFile.openStream(_vcfFilename);

        //enable parsers
        _vcfFile.enablePositionParsing(); //CHROM, POS
        _vcfFile.enableVariantParsing(); //REF, ALT
        _vcfFile.enableVariantQualityParsing(); //QUAL
        //_vcfFile.enableInfoParsing(); //INFO
        _vcfFile.enableFormatParsing(); //FORMAT
        _vcfFile.enableSampleParsing(); //Samples


        //Match samples
        if (_samples.hasSamples()) {
            _samples.fillVCFOrder(_vcfFile.parser.samples);
        } else {
            _samples.readSamplesFromVCFNames(_vcfFile.parser.samples);
        }

        //resize populations
        //_populations.resize(_samples.numPopulations());

        //get output name
        std::string tmp = coretools::str::extractBeforeLast(_vcfFilename, ".vcf");
        _outname = Parameters.getParameterWithDefault<std::string>("out", tmp);

        //limit lines - testing
        _limitLines = Parameters.parameterExists("limitLines");
        if (_limitLines) {
            _maxNumLines = Parameters.getParameter<long>("limitLines");
        } else {
            _maxNumLines = 0;
        }
    };


    void TF2Estimator::calculateF2() {
        //open output file
        std::string filename = _outname + "_F2.txt";
        _logfile->list("Writing F2 results to file '" + filename + "'.");
        std::vector<std::string> header; //header.emplace_back("Sample");
        for (uint32_t s = 0; s < _samples.numSamples(); ++s) {
            header.emplace_back(_samples.sampleName(s));
        }
        TOutputFile out(filename, header);

        //progress
        coretools::TTimer timer;
        uint64_t numFiltered = 0; //# multiallelic sites filtered
        uint64_t lineCounter = 0; //# lines

        //matrix total comparison and diff sites
        std::vector<uint64_t> vec_F2 (_samples.numSamples() * _samples.numSamples());

        //traverse VCF
        _logfile->startIndent("Traversing VCF file:");
        while (_vcfFile.next()) {
            ++lineCounter;
            //exclude multiallelic
            if (_vcfFile.getNumAlleles() == 2) {
                //compare genotypes one ind vs all others at current line
                for (uint32_t s1 = 0; s1 < _samples.numSamples()-1; ++s1) {
                    uint32_t vcfIndex = _samples.sampleIndexInVCF(s1);
                    if (!_vcfFile.sampleIsMissing(vcfIndex)) {
                        auto genotype_s1 = _vcfFile.sampleBiallelicGenotype(vcfIndex);

                        for (uint32_t s2 = s1+1; s2 < _samples.numSamples(); ++s2) {
                            vcfIndex = _samples.sampleIndexInVCF(s2);
                            if (!_vcfFile.sampleIsMissing(vcfIndex)) {
                                auto genotype_s2 = _vcfFile.sampleBiallelicGenotype(vcfIndex);
                                //store total # comparison per combination lower triangle
                                ++vec_F2[(s2*_samples.numSamples())+s1];
                                if (genotype_s1 != genotype_s2) { //genometools::BiallelicGenotype::haploidFirst
                                    //store # diff Sites per combination upper triangle
                                    ++vec_F2[(s1*_samples.numSamples())+s2];
                                }
                            }
                        }
                    }
                }

                //perform test
                //_populations.runTest(out);

                //progress / limit lines
                if (lineCounter % 10000 == 0) {
                    _logfile->list("Parsed ", lineCounter, " lines in " + timer.formattedTime());
                }
                if (_limitLines && lineCounter == _maxNumLines) {
                    break;
                }
            } else {
                ++numFiltered;
            }
        }
        //out.writeLine(std::vector<>);

        _logfile->list("Reached end of VCf file.");
        _logfile->conclude("Parsed ", lineCounter, " lines in ", timer.formattedTime());
        _logfile->conclude("Ignored ", numFiltered, " sites that were not bi-allelic SNPs.");

        std::vector<uint64_t> subvector;
        uint64_t tmp = 0;
        for (uint32_t s = 0; s < _samples.numSamples(); ++s) {
            tmp = s*_samples.numSamples();
            subvector = { vec_F2.begin() + tmp, vec_F2.begin() + tmp + (_samples.numSamples()) };
            out.writeLine(subvector);
        }

        _logfile->endIndent();

        //close file
        out.close();
    };
}