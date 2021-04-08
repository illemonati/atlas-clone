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
    beagleFile = new TOutputFile(_outname + ".beagle.gz");
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
}

TVcfToLFMM::~TVcfToLFMM(){
    lfmmFile->close();
    lociNamesFile->close();
    delete lociNamesFile;
    delete lfmmFile;
}

void TVcfToLFMM::writeHeader(){
    // empty header
}

void TVcfToLFMM::initOutputFiles() {
    lfmmFile = new TOutputFile(_outname + ".lfmm");
    lociNamesFile = new TOutputFile(_outname + ".lfmm.kept_loci");
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

void TVcfToLFMM::prepareAndReadVcf(TParameters & Params){
    // read Vcf and store output
    readVcfAndWriteFile(Params);

    writeLociNames();
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
}

void TVcfToLFMMPostGeno::writeData(TPopulationLikehoodLocus & data){
    // LFMM has individuals as rows and loci as columns -> we need to store these values first and then write
    storePosteriorGenotypes(data);
    storeLocusNames();
}

void TVcfToLFMMPostGeno::storePosteriorGenotypes(TPopulationLikehoodLocus & data){
    auto * meanPostGenoForOneLocus = new float[samples.numSamples()];
    for (uint32_t i = 0; i < samples.numSamples(); i++){
        meanPostGenoForOneLocus[i] = computePosteriorGenotype(data, i);
    }
    genotypes.emplace_back(meanPostGenoForOneLocus);
}

float TVcfToLFMMPostGeno::computePosteriorGenotype(TPopulationLikehoodLocus & data, uint32_t i){
    if (data[i].isMissing)
        throw std::runtime_error("Missing data at sample " + samples.getNameFromOrderedIndex(i) + " and locus " + reader->chr() + ":" + toString(reader->position())
        + "! LFMM2 does not accept missing genotypes, please impute your VCF file first.");
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
    posFile = new TOutputFile(_outname + ".pos");
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
    positionPreviousLocus = 0;
    minDistanceToPreviousLocus = 0;
    numSamplesPerLocus = 0;
    curChr = "";
}

TVcfToGenotypeTruthSetFile::~TVcfToGenotypeTruthSetFile() {
    delete genFile;
    for(uint32_t s = 0; s < samples.numSamples(); s++)
        delete bedFiles[s];
    delete [] bedFiles;
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
    genFile = new TOutputFile(_outname + ".gen");
    // initialize bed files (we now know how many samples there are)
    bedFiles = new TBed * [samples.numSamples()];
    for(uint32_t s = 0; s < samples.numSamples(); s++) {
        bedFiles[s] = new TBed;
    }

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
    if (distanceToPreviousLocus >= minDistanceToPreviousLocus) { // check if distance is big enough
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

void TVcfToGenotypeTruthSetFile::vcfToGenotypeTruthSetFile(TParameters & Params){
    // read params
    minDistanceToPreviousLocus = Params.getParameterIntWithDefault("minDistance", 100);
    logfile->list("Will keep loci that have a minimal distance of " + toString(minDistanceToPreviousLocus) + " to previous locus (parameter 'minDistance').");
    resetDistance();
    numSamplesPerLocus = Params.getParameterIntWithDefault("numSamples", 5);
    logfile->list("Will keep up to " + toString(numSamplesPerLocus) + " individuals per locus (parameter 'numSamples').");

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

/***************************************
 * 					*
 * 	Vcf to Vcf (filter, convert,...)*
 * 					*
 ***************************************/


TVcfToVcf::TVcfToVcf(TParameters &Params, TLog *Logfile) : TVcfConverter(Logfile, Params){
	//open output file
	std::string filename = _outname + (std::string) "_filtered.vcf.gz";
	gz::ogzstream out(filename.c_str());
	if(!out) throw "Failed to open outputfile '" + filename + "'!";
	reader->setOutStream(out);
};

void TVcfToVcf::writeData(TPopulationLikehoodLocus & data){
	reader->writeVCFLine();
};

void TVcfToVcf::writeHeader(){
	reader->writeUnknownHeader();
};
