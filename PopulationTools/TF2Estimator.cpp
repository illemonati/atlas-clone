//
// Created by reynac on 2/22/22.
//

#include "TF2Estimator.h"

auto &parameters = coretools::instances::parameters();
auto &logfile = coretools::instances::logfile();
auto &randomgenerator = coretools::instances::randomGenerator();

namespace PopulationTools {

    using coretools::TParameters;
    using coretools::TLog;
    using coretools::TRandomGenerator;
    using coretools::TOutputFile;
    using coretools::TFactorial;

//------------------------------------------------
//TF2Estimator
//------------------------------------------------
    TF2Estimator::TF2Estimator() {
        _init();
        _initialize();
    };

    void TF2Estimator::_init(){
        _initialized = false;
        //settings
        _limitLines = false;
        _maxNumLines = 0;
        _minDepth = 1;
        _minNumSamplesWithData = 0;
        _minVariantQuality = 0;
    }

    void TF2Estimator::_initialize(){
        // limit lines for testing purpose
        _limitLines = parameters.parameterExists("limitLines");
        if (_limitLines) {
            _maxNumLines = parameters.getParameter<long>("limitLines");
            logfile.list("Will limit analysis to the first " + coretools::str::toString(_maxNumLines) + " lines of the VCF file. (parameter 'limitLines')");
        }

        // set a depth filter
        _minDepth = parameters.getParameterWithDefault<int>("minDepth", 1);
        if(_minDepth < 1)
            throw "minDepth must be >= 1!";
        if(_minDepth > 1)
            logfile.list("Will filter samples to a minimum depth of " + coretools::str::toString(_minDepth) + ". (parameter 'minDepth')");


        // set a missingness filter
        _minNumSamplesWithData = parameters.getParameterWithDefault<uint32_t>("minSamplesWithData", 1);
        if(_minNumSamplesWithData > 1)
            logfile.list("Will remove loci where less than " + coretools::str::toString(_minNumSamplesWithData) + " samples have data. (parameter 'minSamplesWithData')");

        // filter on variant quality
        _minVariantQuality = parameters.getParameterWithDefault<int>("minVariantQuality", 0);
        if(_minVariantQuality < 0) throw "minVariantQuality must be >= 0!";
        if(_minVariantQuality > 0){
            logfile.list("Will only keep sites with variant quality >= " + coretools::str::toString(_minVariantQuality) + ". (parameter 'minVariantQuality')");
        }

        //read list of provided samples
        if (parameters.parameterExists("samples")) {
            _samples.readSamples(parameters.getParameter<std::string>("samples"), &logfile);
        }
        _initialized = true;
    }


    void TF2Estimator::run(){
        _openVCF();
        _readVCF();
    }

    void TF2Estimator::_openVCF(){
        //open VCF
        if(!_initialized){
            throw "Can not open VCF: TPopulationLikelihoodReader was never initialized!";
        }

        _vcfFilename = parameters.getParameter<std::string>("vcf");
        logfile.list("Reading VCF from file '" + _vcfFilename + "'.");
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
        _outname = parameters.getParameterWithDefault<std::string>("out", tmp);
    }

    int TF2Estimator::_filterOnDepth(){
        int numIndividualsWithData = 0;
        for(uint32_t s = 0; s < _samples.numSamples(); ++s){
            int vcfIndex = _samples.sampleIndexInVCF(s);

            // depth filter: if a locus has < minDepth reads, flag locus as missing (set all genotype likelihoods = 1)
            if (_vcfFile.sampleDepth(vcfIndex) < _minDepth){
                _vcfFile.setSampleMissing(vcfIndex);
            } else {
                numIndividualsWithData++;
            }
        }
        return numIndividualsWithData;
    };

    bool TF2Estimator::_filterSite(){
        //skip sites with != 2 alleles
        if(_vcfFile.getNumAlleles() != 2){
            return false;
        }

        //skip sites with too low variant quality
        if(_minVariantQuality > 0 && (_vcfFile.variantQualityIsMissing() || _vcfFile.variantQuality() < _minVariantQuality)){
            return false;
        }

        // missingness filter: if > percentMissingPerLocus of individuals per locus have are missing, remove locus
        uint32_t numIndividualsWithData = _filterOnDepth();
        if (numIndividualsWithData < _minNumSamplesWithData){
            return false;
        }

        //all fine!
        return true;
    };


    void TF2Estimator::_readVCF(){
        //progress
        coretools::TTimer timer;
        uint64_t numFiltered = 0; // # filtered sites
        uint64_t lineCounter = 0; // # lines

        //caclulate total comparisons and diff sites
        std::vector<uint64_t> countsDiff (_samples.numSamples() * _samples.numSamples());
        //traverse VCF
        logfile.startIndent("Traversing VCF file:");
        while (_vcfFile.next()) {
            ++lineCounter;
            //exclude sites
            if (_filterSite()) {
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
                    logfile.list("Parsed ", lineCounter, " lines in " + timer.formattedTime());
                }
                if (_limitLines && lineCounter == _maxNumLines) {
                    break;
                }
            } else {
                ++numFiltered;
            }
        }

        logfile.list("Reached end of VCf file.");
        logfile.conclude("Parsed ", lineCounter, " lines in ", timer.formattedTime());
        logfile.conclude("Ignored ", numFiltered, " sites that were not bi-allelic SNPs.");

        calculateF2(countsDiff);
    }

    void TF2Estimator::calculateF2(const std::vector<uint64_t>& countsDiff) {
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
                    if(counts > 0){ popF2[p1 * _samples.numPopulations() + p2] /= counts; }
                } else{
                    uint32_t  counts = 0;
                    for(uint32_t s1 = _samples.startIndex(p1); s1 < _samples.startIndex(p1) + _samples.numSamplesInPop(p1); ++s1){
                        for(uint32_t s2 = _samples.startIndex(p2); s2 < _samples.startIndex(p2) + _samples.numSamplesInPop(p2); ++s2){
                            popF2[p1 * _samples.numPopulations() + p2] += sampleF2[s1 * _samples.numSamples() + s2];
                            ++counts;
                        }
                    }
                    popF2[p1 * _samples.numPopulations() + p2] /= counts;
                    popF2[p2 * _samples.numPopulations() + p1] = popF2[p1 * _samples.numPopulations() + p2];
                }

            }
        }
        writeF2(countsDiff, sampleF2, popF2);
    };

    void TF2Estimator::writeF2(const std::vector<uint64_t>& countsDiff, const std::vector<double>& sampleF2, const std::vector<double>& popF2){
        //open output file for counts
        std::string filename = _outname + "_counts.txt";
        logfile.listFlush("Writing counts of different/compared sites in the upper/lower triangle to file '" + filename + "' ...");
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
        logfile.done();

        //Sample F2
        filename = _outname + "_sampleF2.txt";
        logfile.listFlush("Writing sample F2 results to file '" + filename + "' ...");
        //open output file for sample F2
        TOutputFile outF2(filename, header);
        for (uint32_t s = 0; s < _samples.numSamples(); ++s) {
            uint64_t tmp = s * _samples.numSamples();
            std::vector<double> subvector = {sampleF2.begin() + tmp, sampleF2.begin() + tmp + (_samples.numSamples()) };
            outF2 << _samples.sampleName(s) << subvector << std::endl;
        }
        logfile.done();

        //Pop F2
        filename = _outname + "_popF2.txt";
        logfile.listFlush("Writing population F2 results to file '" + filename + "' ...");
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
        logfile.done();
        logfile.endIndent();
    }

}