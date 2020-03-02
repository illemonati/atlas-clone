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

bool TVcfReader::readOneLineVcf(TVcfFilters & vcfFilters, TSampleLikelihoods * data, std::string & locusName){
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
        numIndividualsWithData = vcfFilters.filterOnDepth(data, vcfFile, numSamples);
        double tmp[3];
        for(int s = 0; s < numSamples; ++s){
            vcfFile.fillLog10GenotypeLikelihoods(s, tmp[0], tmp[1], tmp[2]);
            data[s].glfLikelihood_0 = glfConverter.log10ToGlfFormat(tmp[0]);
            data[s].glfLikelihood_1 = glfConverter.log10ToGlfFormat(tmp[1]);
            data[s].glfLikelihood_2 = glfConverter.log10ToGlfFormat(tmp[2]);
        }
    } else if(vcfFile.formatColExists("PL")){
        numIndividualsWithData = vcfFilters.filterOnDepth(data, vcfFile, numSamples);
        uint8_t tmp[3];
        for(int s = 0; s < numSamples; ++s){
            vcfFile.fillPhredScore(s, tmp[0], tmp[1], tmp[2]);
            data[s].glfLikelihood_0 = glfConverter.phredToGlfFormat(tmp[0]);
            data[s].glfLikelihood_1 = glfConverter.phredToGlfFormat(tmp[1]);
            data[s].glfLikelihood_2 = glfConverter.phredToGlfFormat(tmp[2]);
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
    if(vcfFilters.filterParameters.freqFilter > 0.0 || vcfFilters.filterParameters.estimateGenotypeFrequencies) {
        genoFrequencies.estimate(data, static_cast<int>(numSamples), glfConverter, vcfFilters.filterParameters.epsilonF);

        if (genoFrequencies.MAF < vcfFilters.filterParameters.freqFilter) {
            vcfFilters.filterCounters.lowFreqSNPCounter++;
            return false;
        }
    }

    //SNP is accepted!
    locusName = vcfFile.chr() + ":" + toString(vcfFile.position());
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

TVcfFilters::TVcfFilters(){
    logfile = nullptr;
    progressFrequency = 0;
}

void TVcfFilters::init(TParameters & Params, TLog * Logfile){
    logfile = Logfile;
    readFilteringOptions(Params);
}

void TVcfFilters::readFilteringOptions(TParameters & Parameters){
    logfile->startIndent("Filter options:");
    int filterCounter = 0;

    // do we limit the lines to read?
    filterParameters.limitLines = Parameters.getParameterLongWithDefault("limitLines", -1);
    if(filterParameters.limitLines > 0) {
        logfile->list("Will limit analysis to the first " + toString(filterParameters.limitLines) +
                      " lines of the VCF file. (parameter 'limitLines')");
        filterCounter++;
    }

    // do we set a depth filter?
    filterParameters.minDepth = Parameters.getParameterIntWithDefault("minDepth", 1);
    if(filterParameters.minDepth < 1)
        throw std::runtime_error("minDepth must be >= 1!");
    if(filterParameters.minDepth > 1) {
        logfile->list("Will filter samples to a minimum depth of " + toString(filterParameters.minDepth) +
                      ". (parameter 'minDepth')");
        filterCounter++;
    }

    // do we set a missingness filter?
    filterParameters.minNumSamplesWithData = Parameters.getParameterIntWithDefault("minSamplesWithData", 1);
    if(filterParameters.minNumSamplesWithData < 1)
        throw std::runtime_error("minNumSamplesWithData must be >= 1!");
    if(filterParameters.minNumSamplesWithData > 1) {
        logfile->list("Will remove loci where less than " + toString(filterParameters.minNumSamplesWithData) +
                      " samples have data. (parameter 'minSamplesWithData')");
        filterCounter++;
    }

    // parameters to set a filter on the allele frequency?
    filterParameters.freqFilter = Parameters.getParameterDoubleWithDefault("minMAF", 0.0); // MAF = minor allele frequency
    if(filterParameters.freqFilter < 0.0 || filterParameters.freqFilter >= 0.5)
        throw std::runtime_error("MAF filter must be within (0.0,0.5)!");
    if(filterParameters.freqFilter > 0.0){
        filterParameters.estimateGenotypeFrequencies = true;
        filterParameters.epsilonF = Parameters.getParameterDoubleWithDefault("epsF", 0.0000001);
        logfile->list("Will filter on an allele frequency of " + toString(filterParameters.freqFilter) + ". (parameter 'minMAF')");
        filterCounter++;
    } else {
        filterParameters.estimateGenotypeFrequencies = false;
    }

    //filter on variant quality?
    filterParameters.minVariantQuality = Parameters.getParameterIntWithDefault("minVariantQuality", 0);
    if(filterParameters.minVariantQuality < 0) throw std::runtime_error("minVariantQuality must be >= 0!");
    if(filterParameters.minVariantQuality > 0){
        logfile->list("Will only keep sites with variant quality >= " + toString(filterParameters.minVariantQuality) + ". (parameter 'minVariantQuality')");
        filterCounter++;
    }

    // filter for specific chromosomes?
    if(Parameters.parameterExists("keepChromosomes")) {
        specifyChromosomesToKeep(Parameters);
        filterCounter++;
    }

    if (filterCounter == 0)
        logfile->list("No filter options set.");

    //set progress frequency
    progressFrequency = Parameters.getParameterIntWithDefault("reportFreq", 10000);
    logfile->endIndent();
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
    logfile->endIndent();
}

int TVcfFilters::filterOnDepth(TSampleLikelihoods* data, TVcfFileSingleLine & vcfFile, unsigned int numSamples){
    int numIndividualsWithData = 0;
    for(int s = 0; s < numSamples; ++s){
        // depth filter: if a locus has < minDepth reads, flag locus as missing (set all genotype likelihoods = 1)
        if (vcfFile.sampleDepth(s) < filterParameters.minDepth)
            vcfFile.setSampleMissing(s);
        else
            numIndividualsWithData++;

        //store if missing and haploid
        data[s].isMissing = vcfFile.sampleIsMissing(s);
        data[s].isHaploid = vcfFile.sampleIsHaploid(s);
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
    vcfFilters.init(Params, logfile);
    readOutputName(Params);
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

void TVcfConverter::readVcfAndWriteFile(){
    //set time at beginning
    logfile->startIndent("Start parsing VCF-file...");
    struct timeval startTime{};
    gettimeofday(&startTime, nullptr);

    // init data for one locus (will be overwritten every line)
    auto * data = new TSampleLikelihoods[vcfReader.numSamples];
    std::string locusName;

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
        if (vcfReader.readOneLineVcf(vcfFilters, data, locusName)){
            // filter passed
            writeData(data, locusName);
        }
    }

    // end of vcf file reached
    vcfFilters.concludeFilters(startTime);
    logfile->endIndent("Reached end of VCF file.");
}

void TVcfConverter::writeData(TSampleLikelihoods * data, const std::string & locusName){
    // will be overwritten by respective function in daughter class
    // because every daughter class will write output in different format
}

int TVcfConverter::baseToNumber(char base, const std::string & marker){
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
TVcfToBeagle::TVcfToBeagle(TParameters &Params, TLog *Logfile) : TVcfConverter(Logfile, Params) {
    beagleFile = nullptr;
}

void TVcfToBeagle::writeBeagleHeader(){
    //header string
    std::vector <std::string> header {"marker", "allele1", "allele2"};
    for(int i=0; i<vcfReader.vcfFile.numSamples(); ++i){
        for(int r=0; r<3; ++r)
            header.push_back(vcfReader.vcfFile.sampleName(i));
    }
    beagleFile->writeHeader(header);
}

void TVcfToBeagle::writeData(TSampleLikelihoods * data, const std::string & locusName){
    //write line
    (*beagleFile) << locusName; // marker
    (*beagleFile) << baseToNumber(vcfReader.vcfFile.getRefAllele(), locusName) << baseToNumber(vcfReader.vcfFile.getFirstAltAllele(), locusName); // ref and alt allele

    for (int s = 0; s < vcfReader.vcfFile.numSamples(); s++){
        if (data[s].isMissing)
            (*beagleFile) << 0.333 << 0.333 << 0.333; // need to do this manually, because otherwise missing data would be 1; but PCAngsd requires genotype likelihoods to sum to one
        else
            (*beagleFile) << glfConverter.toScaledLikelihood(data[s][0]) << glfConverter.toScaledLikelihood(data[s][1]) << glfConverter.toScaledLikelihood(data[s][2]);
    }

    beagleFile->endLine();
}

void TVcfToBeagle::vcfToBeagle(TParameters & Params){
    //open output files
    beagleFile = new TOutputFileZipped(outname + ".beagle.gz");
    writeBeagleHeader();

    // read Vcf and write output
    readVcfAndWriteFile();

    // clean up
    beagleFile->close();
}
/***************************************
 * 									   *
 *    Vcf to LFMM                      *
 * 									   *
 ***************************************/

TVcfToLFMM::TVcfToLFMM(TLog *Logfile, TParameters &Params) : TVcfConverter(Logfile, Params){
    postGeno = false;
    calledGeno = false;
    lfmmFile = nullptr;
    lociNamesFile = nullptr;

    getTask(Params, Logfile);
}

TVcfToLFMM::~TVcfToLFMM(){
    for (auto it = genotypes.begin(); it < genotypes.end(); it++)
        delete [] *it;
    delete lociNamesFile;
    delete lfmmFile;
}

void TVcfToLFMM::writeLFMMHeader(){
    // empty header
    lfmmFile->noHeader(genotypes.size());
}

void TVcfToLFMM::storeLocusNames(const std::string & locusName){
    loci_names.emplace_back(locusName);
}

void TVcfToLFMM::writeLFMM(){
    int numLoci = genotypes.size();
    for (int i = 0; i < vcfReader.vcfFile.numSamples(); i++){
        for (int l = 0; l < numLoci; l++){
            *(lfmmFile) << genotypes[l][i];
        }
        lfmmFile->endLine();
    }
}

void TVcfToLFMM::writeLociNames(){
    lociNamesFile->noHeader(loci_names.size());
    for (auto it = loci_names.begin(); it < loci_names.end(); it++)
        *(lociNamesFile) << *it;
    lociNamesFile->endLine();
}

void TVcfToLFMM::vcfToLFMM(TParameters & Params){
    //open output files
    lfmmFile = new TOutputFilePlain(outname + ".lfmm");
    lociNamesFile = new TOutputFilePlain(outname + ".lfmm.kept_loci");

    // read Vcf and store output
    readVcfAndWriteFile();

    // write actual lfmm
    writeLFMMHeader();
    writeLFMM();
    writeLociNames();

    // clean up
    lfmmFile->close();
    lociNamesFile->close();
}

void TVcfToLFMM::getTask(TParameters &Params, TLog *Logfile){
    std::string task = Params.getParameterString("geno");
    if (task == "postGeno"){
        postGeno = true;
        Logfile->list("Will store posterior genotypes.");
    }
    else if (task == "calledGeno"){
        calledGeno = true;
        Logfile->list("Will store called genotypes.");
    }
    else {
        throw std::runtime_error("Unknown LFMM format " + task + "! Please choose either 'postGeno' or 'calledGeno'.");
    }
}

/***************************************
 * 									   *
 *    Vcf to LFMM (called genotypes)   *
 * 									   *
 ***************************************/
TVcfToLFMMCalledGeno::TVcfToLFMMCalledGeno(TParameters &Params, TLog *Logfile) : TVcfToLFMM(Logfile, Params) {

}

/***************************************
 * 									   *
 * 	Vcf to LFMM (posterior genotype)   *
 * 									   *
 ***************************************/
TVcfToLFMMPostGeno::TVcfToLFMMPostGeno(TParameters &Params, TLog *Logfile) : TVcfToLFMM(Logfile, Params) {

}

void TVcfToLFMMPostGeno::writeData(TSampleLikelihoods * data, const std::string & locusName){
    // LFMM has individuals as rows and loci as columns -> we need to store these values first and then write
    storePosteriorGenotypes(data);
    storeLocusNames(locusName);
}

void TVcfToLFMMPostGeno::storePosteriorGenotypes(TSampleLikelihoods * data){
    auto * meanPostGenoForOneLocus = new double[vcfReader.vcfFile.numSamples()];
    for (int i = 0; i < vcfReader.vcfFile.numSamples(); i++){
        meanPostGenoForOneLocus[i] = computePosteriorGenotype(data, i);
    }
    genotypes.emplace_back(meanPostGenoForOneLocus);
}


double TVcfToLFMMPostGeno::computePosteriorGenotype(TSampleLikelihoods * data, int i){
    // first convert glf to genotype likelihood
    double llG0 = glfConverter.toScaledLikelihood(data[i][0]);
    double llG1 = glfConverter.toScaledLikelihood(data[i][1]);
    double llG2 = glfConverter.toScaledLikelihood(data[i][2]);

    // normalize by sum to get posterior genotype
    double denominator = llG0 + llG1 + llG2;
    double postG0 = llG0 / denominator;
    double postG1 = llG1 / denominator;
    double postG2 = llG2 / denominator;

    // take mean
    double meanPostGeno = postG0 * 0. + postG1 * 1. + postG2 * 2.;
    return meanPostGeno;
}
