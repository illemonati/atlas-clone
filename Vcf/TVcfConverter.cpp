//
// Created by caduffm on 2/24/20.
//

#include <TPopulationLikelihoods.h>
#include "TVcfConverter.h"


/***************************************
 * 									   *
 * 			 TVcfConverter			   *
 * 									   *
 ***************************************/

TVcfConverter::TVcfConverter(TLog * Logfile, TParameters & Params){
    logfile = Logfile;
    reader = nullptr;
    samples = nullptr;
}

TVcfConverter::~TVcfConverter(){
    delete reader;
    delete samples;
}

void TVcfConverter::readOutputName(TParameters & Params, std::string vcfFilename){
    //create out file
    std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
    _outname = Params.getParameterStringWithDefault("out", tmp);

    logfile->list("Writing output files with prefix '" + _outname + "'.");
}

void TVcfConverter::readVcfAndWriteFile(TParameters & Params){
    //read samples
    samples = new TPopulationSamples;
    if(Params.parameterExists("samples"))
        samples->readSamples(Params.getParameterString("samples"), logfile);

    //open VCF reader
    std::string vcfFilename = Params.getParameterString("vcf");
    logfile->startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
    reader = new TPopulationLikelihoodReader(Params, logfile, false);
    reader->openVCF(vcfFilename, logfile);
    logfile->endIndent();

    //Match samples
    if(samples->hasSamples())
        samples->fillVCFOrder(reader->getSampleVCFNames());
    else
        samples->readSamplesFromVCFNames(reader->getSampleVCFNames());
    // write header
    writeHeader();

    //open output file name
    readOutputName(Params, vcfFilename);

    // initialize variables for vcf-file
    struct timeval start{}; gettimeofday(&start, nullptr);
    TPopulationLikehoodLocus data(samples->numSamples());

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader->readDataFromVCF(data, *samples, glfConverter, logfile)){
        writeData(data);
    }

    // end of vcf file reached
    reader->concludeFilters(logfile);
    logfile->endIndent("Reached end of VCF file.");
}

void TVcfConverter::writeHeader(){
    // will be overwritten by respective function in daughter class
    // because every daughter class will write output in different format
}

void TVcfConverter::writeData(TPopulationLikehoodLocus & data){
    // will be overwritten by respective function in daughter class
    // because every daughter class will write output in different format
}

/***************************************
 * 									   *
 * 			 Vcf to Beagle			   *
 * 									   *
 ***************************************/
TVcfToBeagle::TVcfToBeagle(TParameters &Params, TLog *Logfile) : TVcfConverter(Logfile, Params) {
    beagleFile = nullptr;
}

int TVcfToBeagle::baseToNumber(char base){ // TODO: is this implemented somewhere in ATLAS???
    if(base == 'A') return 0;
    else if(base == 'C') return 1;
    else if(base == 'G') return 2;
    else if (base == 'T') return 3;
    else throw std::runtime_error("Unknown base " + toString(base) + "!");
}

void TVcfToBeagle::writeHeader(){
    //header string
    std::vector <std::string> header {"marker", "allele1", "allele2"};
    for(int p=0; p<samples->numPopulations(); p++){
        for(int r=0; r<3; ++r)
            header.push_back(samples->getPopulationName(p));
    }
    beagleFile->writeHeader(header);
}

void TVcfToBeagle::writePosition(){
    (*beagleFile) << toString(reader->position()) + ":" + reader->chr();
}

void TVcfToBeagle::writeRefAndAlt(){
    (*beagleFile) << baseToNumber(reader->refAllele()) << baseToNumber(reader->altAllele());
}

void TVcfToBeagle::writeData(TPopulationLikehoodLocus & data){
    writePosition();
    writeRefAndAlt();
    //write line
    for (int s = 0; s < samples->numSamples(); s++){
        if (data[s].isMissing)
            (*beagleFile) << 0.333 << 0.333 << 0.333; // need to do this manually, because otherwise missing data would be 1; but PCAngsd requires genotype likelihoods to sum to one
        else
            (*beagleFile) << glfConverter.toScaledLikelihood(data[s][0]) << glfConverter.toScaledLikelihood(data[s][1]) << glfConverter.toScaledLikelihood(data[s][2]);
    }

    beagleFile->endLine();
}

void TVcfToBeagle::vcfToBeagle(TParameters & Params){
    //open output files
    beagleFile = new TOutputFileZipped(_outname + ".beagle.gz");

    // read Vcf and write output
    readVcfAndWriteFile(Params);

    // clean up
    beagleFile->close();
}
/***************************************
 * 									   *
 *    Vcf to LFMM                      *
 * 									   *
 ***************************************/

TVcfToLFMM::TVcfToLFMM(TLog *Logfile, TParameters &Params) : TVcfConverter(Logfile, Params){
    lfmmFile = nullptr;
    lociNamesFile = nullptr;
}

TVcfToLFMM::~TVcfToLFMM(){
    for (auto it = genotypes.begin(); it < genotypes.end(); it++)
        delete [] *it;
    delete lociNamesFile;
    delete lfmmFile;
}

void TVcfToLFMM::writeHeader(){
    // empty header
}

void TVcfToLFMM::storeLocusNames(){
    loci_names.emplace_back(toString(reader->position()) + ":" + reader->chr());
}

void TVcfToLFMM::writeLFMM(){
    int numLoci = genotypes.size();
    for (int i = 0; i < samples->numSamples(); i++){
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
    lfmmFile = new TOutputFilePlain(_outname + ".lfmm");
    lociNamesFile = new TOutputFilePlain(_outname + ".lfmm.kept_loci");

    // read Vcf and store output
    readVcfAndWriteFile(Params);

    // write actual lfmm
    lfmmFile->noHeader(genotypes.size()); // we only know now how many loci there are
    writeLFMM();
    writeLociNames();

    // clean up
    lfmmFile->close();
    lociNamesFile->close();
}

/***************************************
 * 									   *
 *    Vcf to LFMM (called genotypes)   *
 * 									   *
 ***************************************/

TVcfToLFMMCalledGeno::TVcfToLFMMCalledGeno(TParameters &Params, TLog *Logfile) : TVcfToLFMM(Logfile, Params) {
    logfile->list("Will store the called genotype for each locus.");
}

void TVcfToLFMMCalledGeno::writeData(TPopulationLikehoodLocus & data){
    // LFMM has individuals as rows and loci as columns -> we need to store these values first and then write
    storeCalledGenotypes();
    storeLocusNames();
}

void TVcfToLFMMCalledGeno::storeCalledGenotypes(){
    auto * calledGenoForOneLocus = new float[samples->numSamples()];
    for (int i = 0; i < samples->numSamples(); i++){
        calledGenoForOneLocus[i] = reader->genotypes(*samples)[i];
    }
    genotypes.emplace_back(calledGenoForOneLocus);
}

/***************************************
 * 									   *
 * 	Vcf to LFMM (posterior genotype)   *
 * 									   *
 ***************************************/

TVcfToLFMMPostGeno::TVcfToLFMMPostGeno(TParameters &Params, TLog *Logfile) : TVcfToLFMM(Logfile, Params) {
    logfile->list("Will store the mean posterior genotype for each locus.");
}

void TVcfToLFMMPostGeno::writeData(TPopulationLikehoodLocus & data){
    // LFMM has individuals as rows and loci as columns -> we need to store these values first and then write
    storePosteriorGenotypes(data);
    storeLocusNames();
}

void TVcfToLFMMPostGeno::storePosteriorGenotypes(TPopulationLikehoodLocus & data){
    auto * meanPostGenoForOneLocus = new float[samples->numSamples()];
    for (int i = 0; i < samples->numSamples(); i++){
        meanPostGenoForOneLocus[i] = computePosteriorGenotype(data, i);
    }
    genotypes.emplace_back(meanPostGenoForOneLocus);
}


float TVcfToLFMMPostGeno::computePosteriorGenotype(TPopulationLikehoodLocus & data, int i){
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
    auto meanPostGeno = static_cast<float>(postG0 * 0. + postG1 * 1. + postG2 * 2.);
    return meanPostGeno;
}