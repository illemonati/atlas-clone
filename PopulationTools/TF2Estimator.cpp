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
        //progress
        coretools::TTimer timer;
        uint64_t numFiltered = 0; // # multiallelic sites filtered
        uint64_t lineCounter = 0; // # lines


        //caclulate total comparisons and diff sites
        std::vector<uint64_t> countsDiff (_samples.numSamples() * _samples.numSamples());
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
                                ++countsDiff[(s2 * _samples.numSamples()) + s1];
                                if (genotype_s1 != genotype_s2) {
                                    //store # diff Sites per combination upper triangle
                                    ++countsDiff[(s1 * _samples.numSamples()) + s2];
                                }
                            }
                        }
                    }
                }

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

        _logfile->list("Reached end of VCf file.");
        _logfile->conclude("Parsed ", lineCounter, " lines in ", timer.formattedTime());
        _logfile->conclude("Ignored ", numFiltered, " sites that were not bi-allelic SNPs.");

        //open output file for counts
        std::string filename = _outname + "_counts.txt";
        _logfile->listFlush("Writing counts of different/compared sites in the upper/lower triangle to file '" + filename + "' ...");
        std::vector<std::string> header = {"Sample"};
        for (uint32_t s = 0; s < _samples.numSamples(); ++s) {
            header.emplace_back(_samples.sampleName(s));
        }
        TOutputFile out(filename, header);

        for (uint32_t s = 0; s < _samples.numSamples(); ++s) {
            uint64_t tmp = s * _samples.numSamples();
            std::vector<uint64_t> subvector = {countsDiff.begin() + tmp, countsDiff.begin() + tmp + (_samples.numSamples()) };
            out << _samples.sampleName(s) << subvector << std::endl;
        }
        _logfile->done();

        //calculate sample F2
        std::vector<double> sampleF2 (_samples.numSamples() * _samples.numSamples());
        for (uint32_t s1 = 0; s1 < _samples.numSamples()-1; ++s1) {
                for (uint32_t s2 = s1+1; s2 < _samples.numSamples(); ++s2) {
                    // diff Sites / total compared sites
                    if ( countsDiff[(s2 * _samples.numSamples()) + s1] != 0){
                        sampleF2[(s1 * _samples.numSamples()) + s2] = static_cast<float>(countsDiff[(s1 * _samples.numSamples()) + s2]) / countsDiff[(s2 * _samples.numSamples()) + s1];
                        sampleF2[(s2 * _samples.numSamples()) + s1] = sampleF2[(s1 * _samples.numSamples()) + s2];
                    }
                }
        }

        filename = _outname + "_sampleF2.txt";
        _logfile->listFlush("Writing sample F2 results to file '" + filename + "' ...");
        //open output file for sample F2
        TOutputFile outF2(filename, header);
        for (uint32_t s = 0; s < _samples.numSamples(); ++s) {
            uint64_t tmp = s * _samples.numSamples();
            std::vector<double> subvector = {sampleF2.begin() + tmp, sampleF2.begin() + tmp + (_samples.numSamples()) };
            outF2 << _samples.sampleName(s) << subvector << std::endl;
        }
        _logfile->done();


            //calculate within and between population average F2
            std::vector<double> popF2 (_samples.numPopulations() * _samples.numPopulations());
            for(uint32_t p1 = 0; p1 < _samples.numPopulations(); ++p1){
                for(uint32_t p2 = p1; p2 < _samples.numPopulations(); ++p2){

                    if ( p1 == p2 ){
                        uint32_t counts = 0;
                        for(uint32_t s1 = _samples.startIndex(p1); s1 < _samples.startIndex(p1) + _samples.numSamplesInPop(p1) - 1; ++s1){
                            for(uint32_t s2 = s1 + 1; s2 < _samples.startIndex(p2) + _samples.numSamplesInPop(p2); ++s2){
                                popF2[p1 * _samples.numPopulations() + p2] += sampleF2[s1 * _samples.numSamples() + s2];
                                ++counts;
                            }
                        }
                        std::cout << counts << std::endl;
                        if(counts > 0){ popF2[p1 * _samples.numPopulations() + p2] /= counts; }
                    } else{
                        uint32_t  counts = 0;
                        for(uint32_t s1 = _samples.startIndex(p1); s1 < _samples.startIndex(p1) + _samples.numSamplesInPop(p1); ++s1){
                            for(uint32_t s2 = _samples.startIndex(p2); s2 < _samples.startIndex(p2) + _samples.numSamplesInPop(p2); ++s2){
                                popF2[p1 * _samples.numPopulations() + p2] += sampleF2[s1 * _samples.numSamples() + s2];
                                ++counts;
                            }
                        }
                        std::cout << counts << std::endl;
                        popF2[p1 * _samples.numPopulations() + p2] /= counts;
                        popF2[p2 * _samples.numPopulations() + p1] = popF2[p1 * _samples.numPopulations() + p2];
                    }

                }
            }

            filename = _outname + "_popF2.txt";
            _logfile->listFlush("Writing population F2 results to file '" + filename + "' ...");
            //open output file for population F2
            std::vector<std::string> Pops = {"Population"};
            for (uint32_t p = 0; p < _samples.numPopulations(); ++p) {
                Pops.emplace_back(_samples.getPopulationName(p));
            }
            TOutputFile outPopF2(filename, Pops);
            for (uint32_t p = 0; p < _samples.numPopulations(); ++p) {
                uint32_t tmp = p * _samples.numPopulations();
                std::vector<double> subvector = {popF2.begin() + tmp, popF2.begin() + tmp + _samples.numPopulations() };
                outPopF2 << _samples.getPopulationName(p) << subvector << std::endl;
            }
            _logfile->done();

        _logfile->endIndent();
    };
}