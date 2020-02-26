//
// Created by caduffm on 2/24/20.
//

#include <TPopulationLikelihoods.h>
#include "TVcfConverter.h"

/***************************************
 * 									   *
 * 			 TVcfReader 			   *
 * 									   *
 ***************************************/

TVcfReader::TVcfReader(){
    logfile = nullptr;
    isZipped = false;
    numSamples = 0;
}

void TVcfReader::init(TParameters & Params, TLog* Logfile){
    logfile = Logfile;

    //open vcf file
    openVCF(Params.getParameterString("vcf"));
}

inline double TVcfReader::phredToError(double phred){
    return pow(10.0, -phred/10.0);
}

//open input stream
void TVcfReader::openVCF(std::string name) {
    //open vcf file
    if (name.find(".gz") == std::string::npos) {
        logfile->list("Reading vcf from file '" + name + "'.");
    } else {
        logfile->list("Reading vcf from gzipped file '" + name + "'.");
        isZipped = true;
    }

    vcfFile.openStream(name, isZipped);
    filename = name;

    //enable parsers
    vcfFile.enablePositionParsing();
    vcfFile.enableVariantParsing();
    vcfFile.enableVariantQualityParsing();
    vcfFile.enableFormatParsing();
    vcfFile.enableSampleParsing();
    vcfFile.enableInfoParsing();

    numSamples = vcfFile.numSamples();
}

bool TVcfReader::readOneLineVcf(TVcfFilters & vcfFilters, double * data){
    // keep chromosomes
    if (!vcfFilters.filterParameters.chromosomesToKeep.empty() && // don't keep this chromosome
        std::find(vcfFilters.filterParameters.chromosomesToKeep.begin(), vcfFilters.filterParameters.chromosomesToKeep.end(), vcfFile.chr()) == vcfFilters.filterParameters.chromosomesToKeep.end()){
        vcfFilters.filterCounters.notOnChrCounter ++;
        return false;
    }

    //skip sites with != 2 alleles
    if(vcfFile.getNumAlleles() != 2){
        vcfFilters.filterCounters.notBialleleicCounter++;
        return false;
    }

    //skip sites with too low variant quality
    if(vcfFilters.filterParameters.minVariantQuality > 0 && (vcfFile.variantQualityIsMissing() || vcfFile.variantQuality() < vcfFilters.filterParameters.minVariantQuality)){
        vcfFilters.filterCounters.lowVariantQualityCounter++;
        return false;
    }

    //check if GL or PL is given
    int numIndividualsWithData = 0;
    if(vcfFile.formatColExists("GL")){
        numIndividualsWithData = vcfFilters.filterOnDepth(numSamples, vcfFile);
        double tmp[3];
        unsigned int index = 0;
        for(int s = 0; s < numSamples; ++s, index+=3){
            if (vcfFile.sampleIsMissing(s)){
                data[index] = 0.333; data[index + 1] = 0.333; data[index + 2] = 0.333;
            }
            else{
                vcfFile.fillLog10GenotypeLikelihoods(s, tmp[0], tmp[1], tmp[2]); // log10 -> take 10^
                data[index] = pow(10, tmp[0]);
                data[index + 1] = pow(10, tmp[1]);
                data[index + 2] = pow(10, tmp[2]);
            }
        }
    } else if(vcfFile.formatColExists("PL")){
        numIndividualsWithData = vcfFilters.filterOnDepth(numSamples, vcfFile);
        uint8_t tmp[3];
        unsigned int index = 0;
        for(int s = 0; s < numSamples; ++s, index+=3){
            if (vcfFile.sampleIsMissing(s)){
                data[index] = 0.333; data[index + 1] = 0.333; data[index + 2] = 0.333;
            }
            else{
                vcfFile.fillPhredScore(s, tmp[0], tmp[1], tmp[2]);
                data[index] = phredToError(tmp[0]);
                data[index + 1] = phredToError(tmp[1]);
                data[index + 2] = phredToError(tmp[2]);
            }
        }
    } else {
        ++vcfFilters.filterCounters.noPLCounter;
        return false;
    }

    // missingness filter: if less than <minNumSamplesWithData> individuals per locus have are missing, remove locus
    if (numIndividualsWithData < vcfFilters.filterParameters.minNumSamplesWithData){
        vcfFilters.filterCounters.missingSNPCounter++;
        return false;
    }

    //filter in MAF
    /*if(vcfFilters.filterParameters.freqFilter > 0.0 || vcfFilters.filterParameters.estimateGenotypeFrequencies) {
        genoFrequencies.estimate(data, numSamples, glfConverter, epsilonF);

        if (genoFrequencies.MAF < freqFilter) {
            vcfFilters.filterCounters.lowFreqSNPCounter++;
            return false;
        }
    }*/

    //SNP is accepted!
    vcfFilters.filterCounters.numAcceptedLoci++;
    return true;
}

/***************************************
 * 									   *
 * 			 TVcfFiltering			   *
 * 									   *
 ***************************************/

TFilterCounters::TFilterCounters() {
    lineCounter = 0;
    notBialleleicCounter = 0;
    missingSNPCounter = 0;
    lowFreqSNPCounter = 0;
    lowVariantQualityCounter = 0;
    noPLCounter = 0;
    notOnChrCounter = 0;
    numAcceptedLoci = 0;
}

TFilterParameters::TFilterParameters() {
    limitLines = 0;
    minDepth = 1;
    minNumSamplesWithData = 0;
    minVariantQuality = 0;
    freqFilter = 0.0;
    epsilonF = 0; //F for EM algorithm to estimate allele frequencies
    estimateGenotypeFrequencies = false;
}

TVcfFilters::TVcfFilters(TParameters & Params, TLog * Logfile){
    progressFrequency = 0;
    logfile = Logfile;

    readFilteringOptions(Params);
}

void TVcfFilters::readFilteringOptions(TParameters & Parameters){
    // do we limit the lines to read?
    filterParameters.limitLines = Parameters.getParameterLongWithDefault("limitLines", -1);
    if(filterParameters.limitLines > 0)
        logfile->list("Will limit analysis to the first " + toString(filterParameters.limitLines) + " lines of the VCF file. (parameter 'limitLines')");

    // do we set a depth filter?
    filterParameters.minDepth = Parameters.getParameterIntWithDefault("minDepth", 1);
    if(filterParameters.minDepth < 1)
        throw std::runtime_error("minDepth must be >= 1!");
    if(filterParameters.minDepth > 1)
        logfile->list("Will filter samples to a minimum depth of " + toString(filterParameters.minDepth) + ". (parameter 'minDepth')");

    // do we set a missingness filter?
    filterParameters.minNumSamplesWithData = Parameters.getParameterIntWithDefault("minSamplesWithData", 1);
    if(filterParameters.minNumSamplesWithData < 1)
        throw std::runtime_error("minNumSamplesWithData must be >= 1!");
    if(filterParameters.minNumSamplesWithData > 1)
        logfile->list("Will remove loci where less than " + toString(filterParameters.minNumSamplesWithData) + " samples have data. (parameter 'minSamplesWithData')");

    // parameters to set a filter on the allele frequency?
    filterParameters.freqFilter = Parameters.getParameterDoubleWithDefault("minMAF", 0.0); // MAF = minor allele frequency
    if(filterParameters.freqFilter < 0.0 || filterParameters.freqFilter >= 0.5)
        throw std::runtime_error("MAF filter must be within (0.0,0.5)!");
    if(filterParameters.freqFilter > 0.0){
        filterParameters.estimateGenotypeFrequencies = true;
        filterParameters.epsilonF = Parameters.getParameterDoubleWithDefault("epsF", 0.0000001);
        logfile->list("Will filter on an allele frequency of " + toString(filterParameters.freqFilter) + ". (parameter 'epsF')");
    } else {
        filterParameters.estimateGenotypeFrequencies = false;
    }

    //filter on variant quality?
    filterParameters.minVariantQuality = Parameters.getParameterIntWithDefault("minVariantQuality", 0);
    if(filterParameters.minVariantQuality < 0) throw std::runtime_error("minVariantQuality must be >= 0!");
    if(filterParameters.minVariantQuality > 0){
        logfile->list("Will only keep sites with variant quality >= " + toString(filterParameters.minVariantQuality) + ". (parameter 'minVariantQuality')");
    }

    // filter for specific chromosomes?
    if(Parameters.parameterExists("keepChromosomes"))
        specifyChromosomesToKeep(Parameters);

    //set progress frequency
    progressFrequency = Parameters.getParameterIntWithDefault("reportFreq", 10000);
}

void TVcfFilters::printProgressFrequencyFiltering(const struct timeval & startTime){
    struct timeval end{};
    gettimeofday(&end, nullptr);
    double runtime = static_cast<double>(end.tv_sec  - startTime.tv_sec)/60.0;
    logfile->list("Parsing line " + toString(filterCounters.lineCounter) + ", retained " + toString(filterCounters.numAcceptedLoci) + " loci in " + toString(runtime) + " min");
}

void TVcfFilters::specifyChromosomesToKeep(TParameters & Parameters){
    std::string argument = Parameters.getParameterString("keepChromosomes");
    if(stringContains(argument, ".txt")){ // specified as a file name
        logfile->startIndent("Reading chromosomes that should be kept from '" + argument + "'");
        std::ifstream keepChromosomesFile(argument.c_str());
        if(!keepChromosomesFile)
            throw std::runtime_error("Failed to open file '" + argument + "'!");
        while(keepChromosomesFile.good() && !keepChromosomesFile.eof()){
            std::string line;
            std::getline(keepChromosomesFile, line);
            std::vector<std::string> vec;
            fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
            //skip empty lines
            if(!vec.empty())
                filterParameters.chromosomesToKeep.push_back(vec[0]);
        }
        keepChromosomesFile.close();
    }
    else { // specified as a vector on command line
        logfile->startIndent("Reading chromosomes from command line.");
        fillVectorFromString(Parameters.getParameterString("keepChromosomes"), filterParameters.chromosomesToKeep, ',');
    }

    // write to logfile
    logfile->startIndent("Will keep the following chromosomes in the output file:");
    for (auto it = filterParameters.chromosomesToKeep.begin(); it < filterParameters.chromosomesToKeep.end(); it ++)
        logfile->list(*it);

    logfile->endIndent();
    logfile->endIndent();
}

void TVcfFilters::concludeFilters(const struct timeval & startTime){
    printProgressFrequencyFiltering(startTime);
    if(filterCounters.notBialleleicCounter > 0)
        logfile->conclude(toString(filterCounters.notBialleleicCounter) + " loci were not bi-allelic.");
    if(filterCounters.lowVariantQualityCounter > 0)
        logfile->conclude(toString(filterCounters.lowVariantQualityCounter) + " loci had variant quality < " + toString(filterParameters.minVariantQuality) + ".");
    if(filterCounters.noPLCounter > 0)
        logfile->conclude(toString(filterCounters.noPLCounter) + " loci had no PL or GL field.");
    if(filterCounters.missingSNPCounter > 0)
        logfile->conclude(toString(filterCounters.missingSNPCounter) + " loci had < " + toString(filterParameters.minNumSamplesWithData) + " samples with data.");
    if(filterCounters.notOnChrCounter > 0)
        logfile->conclude(toString(filterCounters.notOnChrCounter) + " loci were on other chromosomes than specified.");
    if(filterCounters.lowFreqSNPCounter > 0)
        logfile->conclude(toString(filterCounters.lowFreqSNPCounter) + " loci had MAF < " + toString(filterParameters.freqFilter) + ".");
}

int TVcfFilters::filterOnDepth(unsigned int numSamples, TVcfFileSingleLine & vcfFile){
    int numIndividualsWithData = 0;
    unsigned int index = 0;
    for(int s = 0; s < numSamples; ++s, index+=3){
        // depth filter: if a locus has < minDepth reads, flag locus as missing (set all genotype likelihoods = 0.333)
        if (vcfFile.sampleDepth(s) < filterParameters.minDepth)
            vcfFile.setSampleMissing(s);
        else
            numIndividualsWithData++;
    }

    return numIndividualsWithData;
}

/***************************************
 * 									   *
 * 			 TVcfConverter			   *
 * 									   *
 ***************************************/

TVcfConverter::TVcfConverter(TLog * Logfile, TParameters & Params){
    logfile = Logfile;

    vcfReader.init(Params, logfile);
    readOutputName(Params);
    readFilters(Params);
}

TVcfConverter::~TVcfConverter()= default;

void TVcfConverter::readOutputName(TParameters & Params){
    //outputname
    outname = Params.getParameterStringWithDefault("out", "");
    if(outname.empty()){
        //guess from filename
        outname = vcfReader.filename;
        outname = extractBeforeLast(outname, ".");
        if(vcfReader.isZipped)
            //if zipped there is extra .gz
            outname = extractBeforeLast(outname, ".");

    }
    logfile->list("Writing output files with prefix '" + outname + "'.");
}

void TVcfConverter::readFilters(TParameters & Params) {
    // read vcf filters
    TVcfFilters vcfFilters(Params, logfile);
}

bool TVcfConverter::readVcf(TVcfFilters & vcfFilters){
    //set time at beginning
    logfile->startIndent("Start parsing VCF-file...");
    struct timeval startTime{};
    gettimeofday(&startTime, nullptr);

    // init data for one locus (will be overwritten every line)
    double data[3*vcfReader.numSamples];

    //read next
    while(vcfReader.vcfFile.next()){ // new locus
        ++vcfFilters.filterCounters.lineCounter;
        // limit lines
        if(vcfFilters.filterParameters.limitLines > 0 && vcfFilters.filterCounters.lineCounter > vcfFilters.filterParameters.limitLines){
            logfile->list("Reached limit of " + toString(vcfFilters.filterParameters.limitLines) + " lines.");
            break;
        }
        //print progress
        if(vcfFilters.filterCounters.lineCounter % vcfFilters.progressFrequency == 0)
            vcfFilters.printProgressFrequencyFiltering(startTime);

        // filter, keep or not
        if (vcfReader.readOneLineVcf(vcfFilters, data)){
            // filter passed
            writeData(data);
        }
    }

    //return false at end of file
    logfile->endIndent("Reached end of VCF file.");
    vcfFilters.concludeFilters(startTime);

    return false;
}

void TVcfConverter::writeData(double * data){}

int TVcfConverter::baseToNumber(char base, std::string & marker){
    if(base == 'A') return 0;
    else if(base == 'C') return 1;
    else if(base == 'G') return 2;
    else if (base == 'T') return 3;
    else throw std::runtime_error("unknown base " + toString(base) + " at marker " + marker);
}


/***************************************
 * 									   *
 * 			 Vcf to Beagle			   *
 * 									   *
 ***************************************/
TVcfToBeagle::TVcfToBeagle(TParameters &Params, TLog *Logfile) : TVcfConverter(Logfile, Params) {}

void TVcfToBeagle::writeBeagleHeader(){
    //header string
    std::vector <std::string> header {"marker", "allele1", "allele2"};
    for(int i=0; i<vcfReader.vcfFile.numSamples(); ++i){
        for(int r=0; r<3; ++r)
            header.push_back(vcfReader.vcfFile.sampleName(i));
    }
    beagleFile->writeHeader(header);
}

void TVcfToBeagle::writeData(double * data){
    //write line
    std::string marker = vcfReader.vcfFile.chr() + ":" + toString(vcfReader.vcfFile.position());
    (*beagleFile) << marker; // marker
    (*beagleFile) << baseToNumber(vcfReader.vcfFile.getRefAllele(), marker) << baseToNumber(vcfReader.vcfFile.getFirstAltAllele(), marker); // ref and alt allele

    for (int s = 0; s < vcfReader.vcfFile.numSamples() * 3; s+=3){
        (*beagleFile) << data[s] << data[s + 1] << data[s + 2];
    }

    beagleFile->endLine();
}

void TVcfToBeagle::vcfToBeagle(TParameters & Params){
    //open output files
    beagleFile = new TOutputFileZipped(outname + ".beagle.gz");
    writeBeagleHeader();

    // clean up
    beagleFile->close();
}
