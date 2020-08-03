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

    //open output file name (has to be done before readVcfAndWriteFile)
    readOutputName(Params);
}

TVcfConverter::~TVcfConverter(){
    delete reader;
}

void TVcfConverter::readOutputName(TParameters & Params){
    //create out file
    std::string vcfFilename = Params.getParameterString("vcf");
    std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
    _outname = Params.getParameterStringWithDefault("out", tmp);

    logfile->list("Writing output files with prefix '" + _outname + "'.");
}

void TVcfConverter::readVcfAndWriteFile(TParameters & Params){
    //read samples
    if(Params.parameterExists("samples"))
        samples.readSamples(Params.getParameterString("samples"), logfile);

    //open VCF reader
    std::string vcfFilename = Params.getParameterString("vcf");
    reader = new TPopulationLikelihoodReaderLocus(Params, logfile, false);
    reader->openVCF(vcfFilename);
    logfile->endIndent();

    //Match samples
    if(samples.hasSamples())
        samples.fillVCFOrder(reader->getSampleVCFNames());
    else
        samples.readSamplesFromVCFNames(reader->getSampleVCFNames());
    // write header
    initOutputFiles();
    writeHeader();

    // initialize variables for vcf-file
    struct timeval start{}; gettimeofday(&start, nullptr);
    TPopulationLikehoodLocus data(samples.numSamples());

    //run through VCF file
    logfile->list("Parsing VCF file...");
    while(reader->readDataFromVCF(data, samples, glfConverter)){
        writeData(data);
    }

    // end of vcf file reached
    reader->concludeFilters();
}

void TVcfConverter::writeHeader(){
    // will be overwritten by respective function in daughter class
    // because every daughter class will write output in different format
}

void TVcfConverter::initOutputFiles(){
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
    baseToNumber['A'] = 0;
    baseToNumber['C'] = 1;
    baseToNumber['G'] = 2;
    baseToNumber['T'] = 3;
}

TVcfToBeagle::~TVcfToBeagle() {
    delete beagleFile;
}

void TVcfToBeagle::initOutputFiles(){
    beagleFile = new TOutputFileZipped(_outname + ".beagle.gz");
}

void TVcfToBeagle::writeHeader(){
    //header string
    std::vector <std::string> header {"marker", "allele1", "allele2"};
    for(uint32_t s = 0; s < samples.numSamples(); s++){
        for(int r=0; r<3; ++r)
            header.push_back(samples.getNameFromOrderedIndex(s));
    }
    beagleFile->writeHeader(header);
}

void TVcfToBeagle::writePosition(){
    (*beagleFile) << reader->chr() + "_" + toString(reader->position());
}

void TVcfToBeagle::writeRefAndAlt(){
    (*beagleFile) << baseToNumber.find(reader->refAllele())->second << baseToNumber.find(reader->altAllele())->second;
}

void TVcfToBeagle::writeData(TPopulationLikehoodLocus & data){
    writePosition();
    writeRefAndAlt();
    //write line
    for (uint32_t s = 0; s < samples.numSamples(); s++){
        if (data[s].isMissing)
            (*beagleFile) << 0.333 << 0.333 << 0.333; // need to do this manually, because otherwise missing data would be 1; but PCAngsd requires genotype likelihoods to sum to one
        else if (data[s].isHaploid)
            (*beagleFile) << glfConverter.toScaledLikelihood(data[s][0]) << glfConverter.toScaledLikelihood(data[s][1]) << 0; // if haploid, the 3rd gtl should be 0 according to ANGSD (without doing this manually, it would just be extremely small, but not exactly 0).
        else
            (*beagleFile) << glfConverter.toScaledLikelihood(data[s][0]) << glfConverter.toScaledLikelihood(data[s][1]) << glfConverter.toScaledLikelihood(data[s][2]);
    }

    beagleFile->endLine();
}

void TVcfToBeagle::vcfToBeagle(TParameters & Params){
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
    sampleNamesFile = nullptr;
}

TVcfToLFMM::~TVcfToLFMM(){
    lfmmFile->close();
    lociNamesFile->close();
    delete lociNamesFile;
    delete lfmmFile;
    delete sampleNamesFile;
}

void TVcfToLFMM::writeHeader(){
    // empty header
}

void TVcfToLFMM::initOutputFiles() {
    lfmmFile = new TOutputFilePlain(_outname + ".lfmm");
    lociNamesFile = new TOutputFilePlain(_outname + ".lfmm.kept_loci");
    sampleNamesFile = new TOutputFilePlain(_outname + ".lfmm.kept_samples");
}

void TVcfToLFMM::storeLocusNames(){
    loci_names.emplace_back( reader->chr() + ":" + toString(reader->position()));
}

void TVcfToLFMM::writeLociNames(){
    lociNamesFile->noHeader(loci_names.size());
    for (auto it = loci_names.begin(); it < loci_names.end(); it++)
        *(lociNamesFile) << *it;
    lociNamesFile->endLine();
}

void TVcfToLFMM::writeSampleNames(){
    sampleNamesFile->noHeader(samples.numSamples());
    for (uint32_t i = 0; i < samples.numSamples(); i++)
        *(sampleNamesFile) << samples.getNameFromOrderedIndex(i);
    sampleNamesFile->endLine();
}

void TVcfToLFMM::prepareAndReadVcf(TParameters & Params){
    // read Vcf and store output
    readVcfAndWriteFile(Params);

    writeLociNames();
    writeSampleNames();
}

/***************************************
 * 									   *
 *    Vcf to LFMM (called genotypes)   *
 * 									   *
 ***************************************/

TVcfToLFMMCalledGeno::TVcfToLFMMCalledGeno(TParameters &Params, TLog *Logfile) : TVcfToLFMM(Logfile, Params) {
    logfile->list("Will store the called genotype for each locus.");
}

TVcfToLFMMCalledGeno::~TVcfToLFMMCalledGeno(){
    for (auto it = genotypes.begin(); it < genotypes.end(); it++)
        delete [] *it;
}

void TVcfToLFMMCalledGeno::writeData(TPopulationLikehoodLocus & data){
    // LFMM has individuals as rows and loci as columns -> we need to store these values first and then write
    storeCalledGenotypes();
    storeLocusNames();
}

void TVcfToLFMMCalledGeno::storeCalledGenotypes(){
    auto * calledGeno = new uint8_t[samples.numSamples()];
    reader->fillGenotypes(samples, calledGeno);
    for (uint32_t i = 0; i < samples.numSamples(); i++){
        if (calledGeno[i] == 3)
            calledGeno[i] = 9; // re-code missing genotypes to LFMM format
         // if locus was haploid -> is just 0 or 1 -> no need to treat in special way
    }
    genotypes.emplace_back(calledGeno);
}

void TVcfToLFMMCalledGeno::vcfToLFMM(TParameters & Params){
    prepareAndReadVcf(Params);

    // write actual lfmm
    lfmmFile->noHeader(genotypes.size()); // we only know now how many loci there are
    writeLFMM(genotypes);

    // clean up
    lfmmFile->close();
    lociNamesFile->close();
    sampleNamesFile->close();
}

/***************************************
 * 									   *
 * 	Vcf to LFMM (posterior genotype)   *
 * 									   *
 ***************************************/

TVcfToLFMMPostGeno::TVcfToLFMMPostGeno(TParameters &Params, TLog *Logfile) : TVcfToLFMM(Logfile, Params) {
    logfile->list("Will store the mean posterior genotype for each locus.");
}

TVcfToLFMMPostGeno::~TVcfToLFMMPostGeno(){
    for (auto it = genotypes.begin(); it < genotypes.end(); it++)
        delete [] *it;
}

void TVcfToLFMMPostGeno::vcfToLFMM(TParameters & Params){
    prepareAndReadVcf(Params);

    // write actual lfmm
    lfmmFile->noHeader(genotypes.size()); // we only know now how many loci there are
    writeLFMM(genotypes);

    // clean up
    lfmmFile->close();
    lociNamesFile->close();
    sampleNamesFile->close();
}

void TVcfToLFMMPostGeno::writeData(TPopulationLikehoodLocus & data){
    // LFMM has individuals as rows and loci as columns -> we need to store these values first and then write
    storePosteriorGenotypes(data);
}

void TVcfToLFMMPostGeno::storePosteriorGenotypes(TPopulationLikehoodLocus & data){
    auto * meanPostGenoForOneLocus = new float[samples.numSamples()];
    // first compute mean posterior genotype for all samples with data
    for (uint32_t i = 0; i < samples.numSamples(); i++){
        meanPostGenoForOneLocus[i] = computePosteriorGenotype(data, i);
    }
    // now go over all samples again and impute samples with missing data with the mean of the mean posterior genotype
    float meanPostGeno = calculateMeanOfMeanPosteriorGenotypes(data, meanPostGenoForOneLocus);
    for (uint32_t i = 0; i < samples.numSamples(); i++){
        if (data[i].isMissing)
            meanPostGenoForOneLocus[i] = meanPostGeno;
    }

    // check if all samples have the same meanPosteriorGenotype -> if yes, exclude this locus, as this causes trouble with correlations later
    bool allEqual = true;
    for (int i = 0; i < samples.numSamples() - 1; i++){
        if (meanPostGenoForOneLocus[i] != meanPostGenoForOneLocus[i + 1]) {
            allEqual = false;
            break;
        }
    }
    if (!allEqual){ // only store if not all are equal
        genotypes.emplace_back(meanPostGenoForOneLocus);
        storeLocusNames();
    }
}

float TVcfToLFMMPostGeno::calculateMeanOfMeanPosteriorGenotypes(TPopulationLikehoodLocus & data, const float * meanPostGenoForOneLocus){
    float sum = 0.;
    float c = 0.;
    for (uint32_t i = 0; i < samples.numSamples(); i++){
        if (!data[i].isMissing) {
            sum += meanPostGenoForOneLocus[i];
            c++;
        }
    }
    return sum / c;
}

float TVcfToLFMMPostGeno::computePosteriorGenotype(TPopulationLikehoodLocus & data, uint32_t i){
    if (data[i].isMissing){
        // skip for the moment, will be imputed in second round
        return 1.0;
    }
    // if locus is haploid -> llG2 will be ~0 -> gives same result as if I ignored it -> will not treat haploid data in a special way
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

/***************************************
 * 									   *
 * 	Vcf to posfile (STITCH)            *
 * 									   *
 ***************************************/

TVcfToPosFile::TVcfToPosFile(TParameters &Params, TLog *Logfile) : TVcfConverter(Logfile, Params) {
    // posFile is needed as input for STITCH
    // format:
    //   - tab-separated
    //   - no header
    //   - 4 columns: col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base
    //   - one row per SNP
    // Idea: one posfile per chromosome -> use option keepChromosomes
    posFile = nullptr;
}

TVcfToPosFile::~TVcfToPosFile() {
    delete posFile;
}

void TVcfToPosFile::writeHeader(){
    //no header
    posFile->noHeader(4);
}

void TVcfToPosFile::initOutputFiles() {
    posFile = new TOutputFilePlain(_outname + ".pos");
}

void TVcfToPosFile::writePosition(){
    (*posFile) << reader->chr() << reader->position();
}

void TVcfToPosFile::writeRefAndAlt(){
    (*posFile) << reader->refAllele() << reader->altAllele();
}

void TVcfToPosFile::writeData(TPopulationLikehoodLocus & data){
    writePosition();
    writeRefAndAlt();
    posFile->endLine();
}

void TVcfToPosFile::vcfToPosFile(TParameters & Params){
    // read Vcf and write output
    readVcfAndWriteFile(Params);

    // clean up
    posFile->close();
}

/***************************************
 * 									   *
 * Vcf to genotype truth set (STITCH)  *
 * 									   *
 ***************************************/

TVcfToGenotypeTruthSetFile::TVcfToGenotypeTruthSetFile(TParameters &Params, TLog *Logfile) : TVcfConverter(Logfile, Params) {
    // produces genotype truth set (genfile) for STITCH and bed files for samples
    // idea: first locus -> find 0-n samples that have depth higher than minDepth
    //                   -> write position of this locus into bed-files for these individuals
    //                   -> write genotypes of these individuals to genfile; write genotypes of all other individuals as NA to genfile
    //       next locus  -> is distance to previous locus more than x basepairs?
    //                   -> If yes: find 0-n samples that have depth higher than minDepth
    //                              -> write genotypes of these individuals to genfile; write genotypes of all other individuals as NA to genfile
    //                   -> Else: write genotypes of all individuals as NA to genfile
    // format:
    //   - produces a BED file (3 columns, col 1 = chromosome, col 2 = start (1-based), col 3 = stop) for each individual
    //   - produces a genfile
    //      - tab-separated
    //      - header = sample names from vcf
    //      - one row per SNP
    //      - genotypes encoded as 0,1,2 or NA
    genFile = nullptr;
    bedFiles = nullptr;
    bedFilesOpen = false;
    positionPreviousLocus = 0;
    minDistanceToPreviousLocus = 0;
    numSamplesPerLocus = 0;
    minMAF = 0.;
    curChr = "";
}

TVcfToGenotypeTruthSetFile::~TVcfToGenotypeTruthSetFile() {
    delete genFile;
    if (bedFilesOpen) {
        for (uint32_t s = 0; s < samples.numSamples(); s++)
            delete bedFiles[s];
        delete [] bedFiles;
    }
}

void TVcfToGenotypeTruthSetFile::writeHeader(){
    //header string
    std::vector <std::string> header;
    for(uint32_t s = 0; s < samples.numSamples(); s++){
        header.push_back(samples.getNameFromOrderedIndex(s));
    }
    genFile->writeHeader(header);
}

void TVcfToGenotypeTruthSetFile::initOutputFiles() {
    genFile = new TOutputFilePlain(_outname + ".gen");
    // initialize bed files (we now know how many samples there are)
    bedFiles = new TBed * [samples.numSamples()];
    for(uint32_t s = 0; s < samples.numSamples(); s++) {
        bedFiles[s] = new TBed;
    }
    bedFilesOpen = true;

    // check if minNumSamplesWithData-Filter of TPopulationLikelihoods is zero
    // (even if all samples have missing data, we still need to keep the locus, because otherwise
    // positions in posfile and genfile do not match)
    if (reader->getMinNumSamplesWithData() != 0)
        throw std::runtime_error("Parameter 'minSamplesWithData' must be 0 for this task! Please specify 'minSamplesWithData=0'.");
}

void TVcfToGenotypeTruthSetFile::mapIndividualsToDepth(std::vector<uint32_t> &samplesToKeep) {
    std::map< double, std::vector<uint32_t>, std::greater<double> > depthVsSampleIndexMap; // use depth as key in map -> automatically sorted in descending order
    for (auto & s : samplesToKeep){
        // does depth already exist in map?
        double depth = reader->depth(samples, s);
        auto it = depthVsSampleIndexMap.find(depth);
        if (it == depthVsSampleIndexMap.end()){ // new depth
            std::vector<uint32_t> v = {s};
            depthVsSampleIndexMap[depth] = v;
        } else { // depth already exists
            it->second.push_back(s);
        }
    }
    filterIndividualsWithHighestDepth(samplesToKeep, depthVsSampleIndexMap);
}

void TVcfToGenotypeTruthSetFile::filterIndividualsWithHighestDepth(std::vector<uint32_t> & samplesToKeep, const std::map< double, std::vector<uint32_t>, std::greater<double> > & depthVsSampleIndexMap){
    // only keep top numSamplesPerLocus
    samplesToKeep.clear();
    int c = 0;
    for (auto & it : depthVsSampleIndexMap) { // one depth
        for (auto & it2 : it.second){ // one sample at this depth
            if (c >= numSamplesPerLocus)
                return;
            samplesToKeep.push_back(it2);
            c++;
        }
    }
}

void TVcfToGenotypeTruthSetFile::filterIndividuals(TPopulationLikehoodLocus & data){
    std::vector<uint32_t> samplesToKeep;
    long distanceToPreviousLocus = reader->position() - positionPreviousLocus;

    // does current position also appear in posfile (if given)?
    // If not, ignore it (we only want positions in genfile that are also inside posfile)
    if (!positionsInPosFile.empty()){
        auto it = std::find(positionsInPosFile.begin(), positionsInPosFile.end(), toString(reader->position()));
        if (it == positionsInPosFile.end()) // not found -> don't write to genfile, just ignore
            return;
    }

    if (distanceToPreviousLocus >= minDistanceToPreviousLocus && reader->getMAF(data.samples(), samples, glfConverter) >= minMAF) { // check if distance is big enough and if MAF is ok
        // idea: TPopulationLikelihoods will filter on minDepth and set all samples with < minDepth as missing
        // here, we check how many individuals have > minDepth; we rank them and only keep numSamplesPerLocus of them
        for (uint32_t s = 0; s < samples.numSamples(); ++s) {
            if (!data[s].isMissing) {
                samplesToKeep.push_back(s);
            }
        }
        if (samplesToKeep.empty()) // no individuals have > minDepth
            writeToGenFile(samplesToKeep);
            // no need to write to bed file
        else if (samplesToKeep.size() <= static_cast<unsigned int>(numSamplesPerLocus)) { // keep all of them
            writeToGenFile(samplesToKeep);
            storeInBedFile(samplesToKeep);
        }
        else { // rank according to depth
            mapIndividualsToDepth(samplesToKeep);
            writeToGenFile(samplesToKeep);
            storeInBedFile(samplesToKeep);
        }
    } else {
        writeToGenFile(samplesToKeep);
        // no need to write to bed file
    }
}

void TVcfToGenotypeTruthSetFile::writeToGenFile(const std::vector<uint32_t> & samplesToKeep){
    for(uint32_t s = 0; s < samples.numSamples(); s++) {
        // should we write true genotype of sample?
        auto it = std::find(samplesToKeep.begin(), samplesToKeep.end(), s);
        if (it != samplesToKeep.end()) // sample found
            (*genFile) << static_cast<int>(reader->genotype(samples, s));
        else
            (*genFile) << "NA";
    }
    genFile->endLine();
}

void TVcfToGenotypeTruthSetFile::storeInBedFile(const std::vector<uint32_t> & samplesToKeep){
    for(uint32_t s = 0; s < samples.numSamples(); s++) {
        // should we write to bed of sample?
        auto it = std::find(samplesToKeep.begin(), samplesToKeep.end(), s);
        if (it != samplesToKeep.end()) // sample found
            bedFiles[s]->addSite(reader->position());
    }
    positionPreviousLocus = reader->position();
}

void TVcfToGenotypeTruthSetFile::writeData(TPopulationLikehoodLocus & data){
    if (curChr != reader->chr()){ // new chromosome
        curChr = reader->chr();
        for(uint32_t s = 0; s < samples.numSamples(); s++) {
            bedFiles[s]->setChrCreateIfNew(curChr);
        }
        resetDistance();
    }
    filterIndividuals(data);
}

void TVcfToGenotypeTruthSetFile::resetDistance() {
    positionPreviousLocus = -minDistanceToPreviousLocus - 1; // for first locus on chromosome
}

void TVcfToGenotypeTruthSetFile::readPosfileIntoMemory(std::string fileNamePosfile){
    std::ifstream file(fileNamePosfile.c_str());
    if(!file) throw "Failed to open file '" + fileNamePosfile + " for reading!";

    //tmp variables
    std::vector<std::string> vec;
    std::string line;

    //now parse and store
    while(file.good() && !file.eof()){
        std::getline(file, line);
        fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
        //skip empty lines
        if(vec.size() > 0)
            positionsInPosFile.push_back(vec.at(1));
    }

    //close file
    file.close();
}

void TVcfToGenotypeTruthSetFile::vcfToGenotypeTruthSetFile(TParameters & Params){
    // read params
    minDistanceToPreviousLocus = Params.getParameterIntWithDefault("minDistance", 100);
    logfile->list("Will keep loci that have a minimal distance of " + toString(minDistanceToPreviousLocus) + " to previous locus (parameter 'minDistance').");
    resetDistance();
    numSamplesPerLocus = Params.getParameterIntWithDefault("numSamples", 5);
    logfile->list("Will keep up to " + toString(numSamplesPerLocus) + " individuals per locus (parameter 'numSamples').");
    minMAF = Params.getParameterIntWithDefault("minMAFForGenFile", 0);
    logfile->list("Will keep loci that have a minimal minor allele frequency (MAF) of " + toString(minMAF) + " (parameter 'minMAFForGenFile').");
    std::string fileNamePosfile = Params.getParameterStringWithDefault("posfile", "");
    if (!fileNamePosfile.empty())
        readPosfileIntoMemory(fileNamePosfile);

    // read Vcf and write output
    readVcfAndWriteFile(Params);

    // write bed files (one per sample)
    for(uint32_t s = 0; s < samples.numSamples(); s++) {
        // check if sample name contains / (would be interpreted as path)
        std::string sample_name = samples.getNameFromOrderedIndex(s);
        if (stringContains(sample_name, '/'))
            sample_name = stringReplace("/", "_", sample_name);
        bedFiles[s]->write(_outname + "_" + sample_name + ".bed");
    }

    // clean up
    genFile->close();
}

/****************************************
 * 					                    *
 * 	Vcf to Vcf (filter)                 *
 * 					                    *
 ****************************************/

TVcfToVcf::TVcfToVcf(TParameters &Params, TLog *Logfile) : TVcfConverter(Logfile, Params){
};

void TVcfToVcf::initOutputFiles(){
    //open output file
    std::string filename = _outname + "_filtered.vcf.gz";
    reader->openOutputStream(filename, true);
}

void TVcfToVcf::writeData(TPopulationLikehoodLocus & data){
	reader->writeVCFLine();
};

void TVcfToVcf::writeHeader(){
    // TODO: order of header lines differ from input vcf, as lines are stored in map and are sorted alphabetically -> is this a problem?
	reader->writeUnknownHeader();
};

void TVcfToVcf::vcfToVcf(TParameters & Params){
    // read Vcf and write output
    readVcfAndWriteFile(Params);
}

/***************************************
 * 									   *
 *    Extract from VCF                 *
 * 									   *
 ***************************************/

TVcfExtract::TVcfExtract(TParameters &Params, TLog *Logfile) : TVcfConverter(Logfile, Params) {
    file = nullptr;
}

TVcfExtract::~TVcfExtract() {
    delete file;
}

void TVcfExtract::initOutputFiles(){
    std::string end;
    switch (whatToExtract) {
        case genotypes: end = ".genotypes"; break;
        case depth: end = ".depth"; break;
        case refAndAlt: end = ".refAndAlt"; break;
    }
    file = new TOutputFilePlain(_outname + end);
}

void TVcfExtract::_fillHeaderIndividuals(std::vector <std::string> & header){
    for (uint32_t s = 0; s < samples.numSamples(); s++)
        header.push_back(samples.getNameFromOrderedIndex(s));
}

void TVcfExtract::writeHeader(){
    //header string
    std::vector <std::string> header = {"position"};
    switch (whatToExtract) {
        // individual-specific? -> write header of sample names
        case genotypes:
            _fillHeaderIndividuals(header);
            break;
        case depth:
            _fillHeaderIndividuals(header);
            break;
         // same for all individuals? --> add more cases here, if ever needed
        case refAndAlt:
            header.emplace_back("refAllele"); header.emplace_back("altAllele");
    }
    file->writeHeader(header);
}

void TVcfExtract::writePosition(){
    (*file) << reader->chr() + "_" + toString(reader->position());
}

void TVcfExtract::_writeDepth(uint32_t s){
    (*file) << reader->depth(samples, s);
}

void TVcfExtract::_writeGenotype(uint32_t s){
    int genotype = static_cast<int>(reader->genotype(samples, s));
    if (genotype == 3) // missing
        (*file) << "NA";
    else
        (*file) << genotype;
}

void TVcfExtract::_writeRefAndAlt(){
    (*file) << reader->refAllele() << reader->altAllele();
}

void TVcfExtract::writeData(TPopulationLikehoodLocus & data){
    writePosition();
    switch (whatToExtract) {
        case genotypes:
            for (uint32_t s = 0; s < samples.numSamples(); s++)
                _writeGenotype(s);
            break;
        case depth:
            for (uint32_t s = 0; s < samples.numSamples(); s++)
                _writeDepth(s);
            break;
        case refAndAlt:
            _writeRefAndAlt();
            break;
    }
    file->endLine();
}

void TVcfExtract::vcfExtract(TParameters & Params){
    // what to extract?
    std::string type = Params.getParameterString("extract");
    if (type == "genotypes") {
        whatToExtract = genotypes;
    } else if (type == "depth"){
        whatToExtract = depth;
    } else if (type == "refAndAlt") {
        whatToExtract = refAndAlt;
    } // ... add more stuff here, if ever needed
    else throw "Unknown type '" + type + "' for parameter 'extract'! Specify either 'genotypes' or 'depth' or 'refAndAlt'.";

    // read Vcf and write output
    readVcfAndWriteFile(Params);

    // clean up
    file->close();
}

/***************************************
 * 					                   *
 * 	STITCH output Vcf reader           *
 * 			                    	   *
 ***************************************/

// this is a hack and not very nice - can't use VCF reader without bigger refactoring, because STITCH VCF has
// GP field (Posterior genotype probability) which can't be read by VCF reader...

TStitchVcfReader::TStitchVcfReader(){
    _vcfOpen = false;
    _vcf = nullptr;
    _numSamples = 0;
}

TStitchVcfReader::~TStitchVcfReader() {
    if (_vcfOpen) {
        _vcf->close();
        delete _vcf;
    }
}

void TStitchVcfReader::openVcf(const std::string &vcfFilename) {
    if (!_vcfOpen) {
        _vcf = new gz::igzstream(vcfFilename.c_str());
        _vcfOpen = true;
    }
    if(!_vcf || _vcf->fail() || !_vcf->good()) throw "Failed to open file '" + vcfFilename + "'!";
}

void TStitchVcfReader::parseVCFHeader(std::vector<std::string>& sampleNames){
    std::string temp, buf;
    std::vector<std::string> tmp;

    // read header and stop after that
    while(!_vcf->eof() && _vcf->peek()=='#'){
        getline(*_vcf, temp);
        if(stringContains(temp, "#CHROM")){ // here are the samples
            //analyze header: save samples
            trimString(temp);
            int i=0;
            while(!temp.empty()){
                buf=extractBeforeWhiteSpace(temp);
                trimString(buf);
                temp.erase(0,1);
                if(i >= 9) tmp.push_back(buf); // col 9 is the first col with sample names
                ++i;
            }
        }
    }

    _numSamples = tmp.size();
    sampleNames.insert(sampleNames.end(), tmp.begin(), tmp.end()); // append tmp to sampleNames
}

bool TStitchVcfReader::read(){
    std::string temp;
    do{
        if(_vcf->eof()){
            return false;
        }
        getline(*_vcf, temp);
    } while(temp.empty());

    fillVectorFromStringWhiteSpaceSkipEmpty(temp, _oneLine);
    return true;
}

std::string TStitchVcfReader::chr() {
    return _oneLine[0];
}

std::string TStitchVcfReader::pos() {
    return _oneLine[1];
}

std::string TStitchVcfReader::refAllele() {
    return _oneLine[3];
}

std::string TStitchVcfReader::altAllele() {
    return _oneLine[4];
}

double TStitchVcfReader::infoScore(){
    std::string infoField = _oneLine[7];
    std::string tmp = extractAfter(infoField, ";");
    std::string infoScore = extractBefore(tmp, ";");
    std::string actualScore = extractAfter(infoScore, "=");

    return stringToDouble(actualScore);
}

double TStitchVcfReader::HWE_pVal(){
    std::string infoField = _oneLine[7];
    std::string tmp1 = extractAfter(infoField, ";"); // get rid of first ;
    std::string tmp2 = extractAfter(tmp1, ";"); // get rid of second ;
    std::string pVal = extractBefore(tmp2, ";");
    std::string actualScore = extractAfter(pVal, "=");

    return stringToDouble(actualScore);
}

void TStitchVcfReader::genotypes(std::vector<std::string>& genotypes){
    // parse posterior genotypes
    auto it = _oneLine.begin(); it += 9; // jump to first data col
    for (; it < _oneLine.end(); it++){
        // get GT field (before first :) -> contains something like 0/1
        std::string tmp = extractBefore(*it, ":");

        // separate genotypes
        if (tmp == "./.") // missing
            genotypes.emplace_back("NA");
        else {
            std::string tmp2 = extractBefore(tmp, "/");
            tmp.erase(0, 1); // erase /
            int genotype = stringToInt(tmp) + stringToInt(tmp2);
            genotypes.push_back(toString(genotype));
        }
    }
}

void TStitchVcfReader::posteriorGenotypes(std::vector<std::string>& posteriorGenotypes){
    // parse posterior genotypes
    std::string gt0, gt1, gt2, buf;
    auto it = _oneLine.begin(); it += 9; // jump to first data col
    for (; it < _oneLine.end(); it++){
        // get GP field (in between : and :)
        std::string tmp = extractAfter(*it, ":");
        buf = extractBefore(tmp, ":");

        // separate genotype posteriors
        gt0 = extractBefore(buf, ",");
        buf.erase(0,1); // erase ,
        gt1 = extractBefore(buf, ",");
        buf.erase(0,1); // erase ,
        gt2 = extractBefore(buf, ",");
        buf.erase(0,1); // erase ,

        posteriorGenotypes.push_back(gt0);
        posteriorGenotypes.push_back(gt1);
        posteriorGenotypes.push_back(gt2);
    }
}

void TStitchVcfReader::maxPosteriorGenotypes(std::vector<std::string>& maxPosteriorGenotypes){
    // parse posterior genotypes
    std::vector<std::string> vecPosteriorGenotypes;
    posteriorGenotypes(vecPosteriorGenotypes);

    // now decide which is the maximum posterior genotype for each sample
    for (int i = 0; i < _numSamples*3; i+=3){
        double gt0 = stringToDouble(vecPosteriorGenotypes.at(i));
        double gt1 = stringToDouble(vecPosteriorGenotypes.at(i+1));
        double gt2 = stringToDouble(vecPosteriorGenotypes.at(i+2));

        double max;
        if (gt0 >= gt1 && gt0 >= gt2)
            max = gt0;
        else if (gt1 >= gt0 && gt1 >= gt2)
            max = gt1;
        else
            max = gt2;
        maxPosteriorGenotypes.push_back(toString(max));
    }
}

void TStitchVcfReader::meanPosteriorGenotypes(std::vector<std::string>& meanPosteriorGenotypes){
    // parse posterior genotypes
    std::vector<std::string> vecPosteriorGenotypes;
    posteriorGenotypes(vecPosteriorGenotypes);

    // now compute the mean posterior genotype for each sample
    for (int i = 0; i < _numSamples*3; i+=3){
        double gt0 = stringToDouble(vecPosteriorGenotypes.at(i));
        double gt1 = stringToDouble(vecPosteriorGenotypes.at(i+1));
        double gt2 = stringToDouble(vecPosteriorGenotypes.at(i+2));

        double mean = 0 * gt0 + 1. * gt1 + 2. * gt2;
        meanPosteriorGenotypes.push_back(toString(mean));
    }
}

/***************************************
 * 					                   *
 * 	STITCH output Vcf to beagle        *
 * 			                    	   *
 ***************************************/

TStitchVcfToBeagle::TStitchVcfToBeagle(TParameters &Params, TLog *logfile) {
    baseToNumber["A"] = 0;
    baseToNumber["C"] = 1;
    baseToNumber["G"] = 2;
    baseToNumber["T"] = 3;

    // open input vcf
    std::string vcfFilename = Params.getParameterString("vcf");
    logfile->startIndent("Reading vcf '" + vcfFilename + "'.");
    reader.openVcf(vcfFilename);

    // open output beagle
    std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
    std::string outname = Params.getParameterStringWithDefault("out", tmp);
    logfile->list("Writing output files with prefix '" + outname + "'.");
    beagleFile.open(outname + ".beagle.gz");

    // parse header
    parseVCFHeader();
    parseVCF();

    beagleFile.close();
}

void TStitchVcfToBeagle::parseVCFHeader(){
    std::vector <std::string> sampleNames;
    reader.parseVCFHeader(sampleNames);
    // now write header
    std::vector <std::string> header {"marker", "allele1", "allele2"};
    for (auto & it: sampleNames){
        for (int r = 0; r < 3; r++)
            header.push_back(it);
    }
    beagleFile.writeHeader(header);
}

void TStitchVcfToBeagle::parseVCF(){
    // parse rest
    while (reader.read()){
        beagleFile << reader.chr() + "_" + reader.pos();
        beagleFile << baseToNumber.find(reader.refAllele())->second << baseToNumber.find(reader.altAllele())->second;

        std::vector<std::string> postGeno;
        reader.posteriorGenotypes(postGeno);
        beagleFile.write(postGeno);
        beagleFile.endLine();
    }
}

/***************************************
 * 					                   *
 * 	STITCH output Vcf to posfile       *
 * 			                    	   *
 ***************************************/

TStitchVcfToPosfile::TStitchVcfToPosfile(TParameters &Params, TLog *logfile) {
    // open input vcf
    std::string vcfFilename = Params.getParameterString("vcf");
    logfile->startIndent("Reading vcf '" + vcfFilename + "'.");
    reader.openVcf(vcfFilename);

    // open output posfile
    std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
    std::string outname = Params.getParameterStringWithDefault("out", tmp);
    logfile->list("Writing output files with prefix '" + outname + "'.");
    posfile.open(outname + ".pos");

    // read filters
    minInfoScore = Params.getParameterDoubleWithDefault("minInfoScore", 0.4);
    minHWE_pval = Params.getParameterDoubleWithDefault("minHWE_pval", 0.000001);

    // parse header
    parseVCFHeader();
    parseVCF();

    posfile.close();
}

void TStitchVcfToPosfile::parseVCFHeader(){
    std::vector <std::string> sampleNames;
    reader.parseVCFHeader(sampleNames);
    // don't write to posfile
    posfile.noHeader(4);
}

void TStitchVcfToPosfile::parseVCF(){
    // parse rest
    while (reader.read()){
        if (reader.infoScore() >= minInfoScore && reader.HWE_pVal() >= minHWE_pval){
            posfile << reader.chr() << reader.pos() << reader.refAllele() << reader.altAllele();
            posfile.endLine();
        }
    }
}

/***************************************
 * 					                   *
 * 	STITCH output Vcf extract stuff    *
 * 			                    	   *
 ***************************************/

TStitchVcfExtract::TStitchVcfExtract(TParameters &Params, TLog *logfile) {
    // open input vcf
    std::string vcfFilename = Params.getParameterString("vcf");
    logfile->startIndent("Reading vcf '" + vcfFilename + "'.");
    reader.openVcf(vcfFilename);

    // what to extract?
    std::string type = Params.getParameterString("extract");
    if (type == "genotypes") {
        whatToExtract = genotypes;
    } else if (type == "posteriorGenotypes"){
        whatToExtract = posteriorGenotypes;
    } else if (type == "maxPostGenotype") {
        whatToExtract = maxPostGenotype;
    } // ... add more stuff here, if ever needed
    else throw "Unknown type '" + type + "' for parameter 'extract'! Specify either 'genotypes' or 'posteriorGenotypes' or 'maxPostGenotype'.";


    // open output file
    std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
    std::string outname = Params.getParameterStringWithDefault("out", tmp);
    logfile->list("Writing output files with prefix '" + outname + "'.");
    std::string end;
    switch (whatToExtract) {
        case genotypes: end = ".genotypes"; break;
        case posteriorGenotypes: end = ".posteriorGenotypes"; break;
        case maxPostGenotype: end = ".maxPostGenotype"; break;
    }
    file.open(outname + end);

    // parse header
    parseVCFHeader();
    parseVCF();

    file.close();
}

void TStitchVcfExtract::parseVCFHeader(){
    std::vector <std::string> header = {"position"};
    reader.parseVCFHeader(header);
    // write header
    switch (whatToExtract) {
        case genotypes: file.writeHeader(header); break;
        case maxPostGenotype: file.writeHeader(header); break;
        case posteriorGenotypes:
            // write header where each sample appears 3x
            std::vector<std::string> expandedHeader = {"position"};
            for (int i = 1; i < header.size(); i++){ // skip first element of header (= position)
                for (int r = 0; r < 3; r++)
                    expandedHeader.push_back(header[i]);
            }
            file.writeHeader(expandedHeader);
            break;
    }
}

void TStitchVcfExtract::_writeGenotypes(){
    // fill genotypes
    std::vector<std::string> genotypes;
    reader.genotypes(genotypes);
    // now write
    file.write(genotypes);
}

void TStitchVcfExtract::_writeMaxPosteriorGenotypes(){
    std::vector<std::string> postGeno;
    reader.maxPosteriorGenotypes(postGeno);
    file.write(postGeno);
}

void TStitchVcfExtract::_writePosteriorGenotypes(){
    std::vector<std::string> postGeno;
    reader.posteriorGenotypes(postGeno);
    file.write(postGeno);
}

void TStitchVcfExtract::parseVCF(){
    // parse rest
    while (reader.read()){
        // write position
        file << reader.chr() + "_" + reader.pos();

        switch (whatToExtract) {
            case genotypes:
                _writeGenotypes();
                break;
            case maxPostGenotype:
                _writeMaxPosteriorGenotypes();
                break;
            case posteriorGenotypes:
                _writePosteriorGenotypes();
                break;
        }
        file.endLine();
    }
}

/****************************************************
 * 									                *
 * 	STITCH output Vcf to LFMM (posterior genotypes) *
 * 									                *
 ****************************************************/


TStitchVcfToLFMMPostGeno::TStitchVcfToLFMMPostGeno(TParameters &Params, TLog *logfile) {
    // open input vcf
    std::string vcfFilename = Params.getParameterString("vcf");
    logfile->startIndent("Reading vcf '" + vcfFilename + "'.");
    reader.openVcf(vcfFilename);

    // open output files
    std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
    std::string outname = Params.getParameterStringWithDefault("out", tmp);
    logfile->list("Writing output files with prefix '" + outname + "'.");
    file.open(outname + ".lfmm");
    fileLociNames.open(outname + ".lfmm.kept_loci");
    fileSampleNames.open(outname + ".lfmm.kept_samples");

    // parse header
    parseVCFHeader();
    parseVCF();

    // now write!
    _write();

    file.close();
    fileLociNames.close();
    fileSampleNames.close();
}

void TStitchVcfToLFMMPostGeno::parseVCFHeader(){
    reader.parseVCFHeader(sample_names);
}

void TStitchVcfToLFMMPostGeno::_write() {
    _writeMeanPosteriorGenotypes();
    _writeLociNames();
    _writeSampleNames();
}

void TStitchVcfToLFMMPostGeno::_writeLociNames(){
    fileLociNames.noHeader(loci_names.size());
    for (auto it = loci_names.begin(); it < loci_names.end(); it++)
        fileLociNames << *it;
    fileLociNames.endLine();
}

void TStitchVcfToLFMMPostGeno::_writeSampleNames(){
    fileSampleNames.noHeader(reader.numSamples());

    for (auto it = sample_names.begin(); it < sample_names.end(); it++)
        fileSampleNames << *it;
    fileSampleNames.endLine();
}

void TStitchVcfToLFMMPostGeno::_writeMeanPosteriorGenotypes(){
    // we only know now how many loci there are
    file.noHeader(_genotypes.size());

    int numLoci = _genotypes.size();
    for (uint32_t i = 0; i < reader.numSamples(); i++){
        for (int l = 0; l < numLoci; l++) {
            file << _genotypes[l][i];
        }
        file.endLine();
    }
}

void TStitchVcfToLFMMPostGeno::parseVCF(){
    // parse
    while (reader.read()){
        // store mean posterior genotypes for later
        std::vector<std::string> meanPostGenoForOneLocus;
        reader.meanPosteriorGenotypes(meanPostGenoForOneLocus);
        // check if all samples have the same meanPosteriorGenotypes -> if yes, exclude locus, because this causes problems when calculating correlations later
        if (!std::equal(meanPostGenoForOneLocus.begin() + 1, meanPostGenoForOneLocus.end(), meanPostGenoForOneLocus.begin())){
            //not all equal -> store
            _genotypes.emplace_back(meanPostGenoForOneLocus);
            // write loci names to file
            loci_names.emplace_back(reader.chr() + ":" + reader.pos());
        }
    }
}