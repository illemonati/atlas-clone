//
// Created by caduffm on 2/24/20.
//

#include "TVcfConverter.h"

namespace VCF{

//------------------------------------------
// TVcfConverter
//------------------------------------------
TVcfConverter::TVcfConverter(TLog * Logfile){
    _logfile = Logfile;
};

void TVcfConverter::_readOutputName(TParameters & Params, const std::string & VCFFilename){
    //create out file
    std::string tmp = coretools::str::readBeforeLast(VCFFilename, ".vcf");
    _outname = Params.getParameterWithDefault<std::string>("out", tmp);

    _logfile->list("Writing output files with prefix '" + _outname + "'.");
};

void TVcfConverter::readVcfAndWriteFile(TParameters & Params){
    // define outName
    std::string vcfFilename = Params.getParameterFilename("vcf");
    _readOutputName(Params, vcfFilename);

    // read samples
    if(Params.parameterExists("samples")) {
        _samples.readSamples(Params.getParameter<std::string>("samples"), _logfile);
    }

    // open VCF-reader
    _reader.initialize(Params, _logfile, false);
    _reader.openVCF(vcfFilename);

    // match samples
    if(_samples.hasSamples()) {
        _samples.fillVCFOrder(_reader.getSampleVCFNames());
    } else {
        _samples.readSamplesFromVCFNames(_reader.getSampleVCFNames());
    }

    // write header
    _initOutputFiles();
    _writeHeader();

    // initialize variables for vcf-file
    PopulationTools::TPopulationLikehoodLocus data(_samples.numSamples());

    // run through VCF file
    _logfile->list("Parsing VCF file...");
    while (_reader.readDataFromVCF(data, _samples)){
        _writeData(data);
    };

    // end of vcf file reached
    _reader.concludeFilters();
};

void TVcfConverter::_writeHeader(){
    // can stay empty if there is no header
};

//------------------------------------------
// TVcfToBeagle
//------------------------------------------
TVcfToBeagle::TVcfToBeagle(TLog *Logfile) : TVcfConverter(Logfile) {};

void TVcfToBeagle::_initOutputFiles(){
    _beagleFile.open(_outname + ".beagle.gz");
};

void TVcfToBeagle::_writeHeader(){
    //header string
    std::vector<std::string> header = {"marker", "allele1", "allele2"};
    for (size_t s = 0; s < _samples.numSamples(); s++){
        // each sample is represented three times (for the three _genotypes)
        std::string sampleName = _samples.sampleName(s);
        header.push_back(sampleName);
        header.push_back(sampleName);
        header.push_back(sampleName);
    }
    _beagleFile.writeHeader(header);
};

void TVcfToBeagle::_writePosition(){
    _beagleFile << _reader.chr() + "_" + coretools::str::toString(_reader.position());
};

void TVcfToBeagle::_writeRefAndAlt(){
    _beagleFile << _reader.refAllele() << _reader.altAllele();
};

void TVcfToBeagle::_writeData(PopulationTools::TPopulationLikehoodLocus & data){
    _writePosition();
    _writeRefAndAlt();

    //write line
    for (size_t s = 0; s < _samples.numSamples(); s++){
        if (data[s].isMissing()){
            _beagleFile << 0.333 << 0.333 << 0.333; // need to do this manually, because otherwise missing data would be 1; but PCAngsd requires genotype likelihoods to sum to one
    	} else if (data[s].isHaploid()){
            _beagleFile << (coretools::Probability) data[s][genometools::homoFirst] << (coretools::Probability) data[s][genometools::homoSecond] << 0; // if haploid, the 3rd gtl should be 0 according to ANGSD (without doing this manually, it would just be extremely small, but not exactly 0).
        } else {
            _beagleFile << (coretools::Probability) data[s][genometools::homoFirst] << (coretools::Probability) data[s][genometools::het] << (coretools::Probability) data[s][genometools::homoSecond];
        }
    };

    _beagleFile.endLine();
};

void TVcfToBeagle::vcfToBeagle(TParameters & Params){
    // read Vcf and write output
    readVcfAndWriteFile(Params);

    // clean up
    _beagleFile.close();
};


//------------------------------------------
// TVcfToGeno
//------------------------------------------
TVcfToGeno::TVcfToGeno(TLog *Logfile) : TVcfConverter(Logfile) {
    // geno format: used in LEA package (https://www.rdocumentation.org/packages/LEA/versions/1.4.0/topics/geno)
    // rows = SNPs, cols = individuals (without delim, just pasted together)
    // Each data point represents the number of copies of reference allele. Missing data encoded by 9.
    // Note: we usually define genotype as number of copies of alternative allele -> might be confusing, I stick to
    //       the file format here, i.e. I use number of copies of reference allele.
};

void TVcfToGeno::_initOutputFiles(){
    _genoFile.open(_outname + ".geno");
    _lociNamesFile.open(_outname + ".geno.kept_loci");
};

void TVcfToGeno::_closeOutputFiles() {
    _genoFile.close();
    _lociNamesFile.endLine();
    _lociNamesFile.close();
};

void TVcfToGeno::_writePosition(){
    _lociNamesFile << _reader.chr() + ":" + coretools::toString(_reader.position());
};

void TVcfToGeno::_writeData(PopulationTools::TPopulationLikehoodLocus & data){
    _writePosition();

    //write line
    std::string line;
    for (size_t s = 0; s < _samples.numSamples(); s++){
        if (data[s].isMissing()){
            line += "9";
        } else if (data[s].isHaploid()){
            // get counts of reference allele from counts of alternative allele
            line += coretools::str::toString(1 - _reader.biallelicGenotype(_samples, s).altAlleleCounts());
        } else {
            // get counts of reference allele from counts of alternative allele
            line += coretools::str::toString(2 - _reader.biallelicGenotype(_samples, s).altAlleleCounts());
        }
    }
    _genoFile << line << std::endl;
};

void TVcfToGeno::vcfToGeno(TParameters &Params) {
    // read Vcf and write output
    readVcfAndWriteFile(Params);

    // clean up
    _closeOutputFiles();
};

//------------------------------------------
// TVcfToLFMM
//------------------------------------------
TVcfToLFMM::TVcfToLFMM(TLog *Logfile) : TVcfConverter(Logfile){};

TVcfToLFMM::~TVcfToLFMM(){
    _closeOutputFiles();
};

void TVcfToLFMM::_initOutputFiles() {
    _lfmmFile.open(_outname + ".lfmm");
    _lociNamesFile.open(_outname + ".lfmm.kept_loci");
};

void TVcfToLFMM::_closeOutputFiles() {
    _lfmmFile.close();
    _lociNamesFile.close();
};

void TVcfToLFMM::_storeLocusNames(){
    _loci_names.emplace_back(_reader.chr() + ":" + coretools::toString(_reader.position()));
};

void TVcfToLFMM::_writeLociNames(){
    _lociNamesFile.noHeader(_loci_names.size());
    for (auto & it : _loci_names)
        _lociNamesFile << it;
    _lociNamesFile.endLine();
};

void TVcfToLFMM::_prepareAndReadVcf(TParameters & Params){
    // read Vcf and store output
    readVcfAndWriteFile(Params);

    _writeLociNames();
};

//------------------------------------------
// TVcfToLFMMCalledGeno
//------------------------------------------
TVcfToLFMMCalledGeno::TVcfToLFMMCalledGeno(TLog *Logfile) : TVcfToLFMM(Logfile) {
    _logfile->list("Will store the called genotype for each locus.");
};

TVcfToLFMMCalledGeno::~TVcfToLFMMCalledGeno(){
    for (auto it = _genotypes.begin(); it < _genotypes.end(); it++)
        delete [] *it;
};

void TVcfToLFMMCalledGeno::_writeData(PopulationTools::TPopulationLikehoodLocus &){
    // LFMM has individuals as rows and loci as columns -> we need to store these values first and then write
    storeCalledGenotypes();
    _storeLocusNames();
};

void TVcfToLFMMCalledGeno::storeCalledGenotypes(){
    auto tmp = _reader.biallelicGenotypes(_samples);
    auto* calledGeno = new uint8_t[_samples.numSamples()];
    for (size_t i = 0; i < _samples.numSamples(); i++){
        if (tmp[i].isMissing()){
            calledGeno[i] = 9; // re-code missing _genotypes to LFMM format
        } else {
        	calledGeno[i] = tmp[i].altAlleleCounts();
        }
    }
    _genotypes.emplace_back(calledGeno);
};

void TVcfToLFMMCalledGeno::vcfToLFMM(TParameters & Params){
    _prepareAndReadVcf(Params);

    // write actual lfmm
    _lfmmFile.noHeader(_genotypes.size()); // we only know now how many loci there are
    _writeLFMM(_genotypes);

    _closeOutputFiles();
};

//------------------------------------------
// TVcfToLFMMPostGeno
//------------------------------------------
TVcfToLFMMPostGeno::TVcfToLFMMPostGeno(TLog *Logfile) : TVcfToLFMM(Logfile) {
    _logfile->list("Will store the mean posterior genotype for each locus.");
};

TVcfToLFMMPostGeno::~TVcfToLFMMPostGeno(){
    for (auto it = _genotypes.begin(); it < _genotypes.end(); it++)
        delete [] *it;
};

void TVcfToLFMMPostGeno::vcfToLFMM(TParameters & Params){
    _prepareAndReadVcf(Params);

    // write actual lfmm
    _lfmmFile.noHeader(_genotypes.size()); // we only know now how many loci there are
    _writeLFMM(_genotypes);

    _closeOutputFiles();
};

void TVcfToLFMMPostGeno::_writeData(PopulationTools::TPopulationLikehoodLocus & data){
    // LFMM has individuals as rows and loci as columns -> we need to store these values first and then write
    _storePosteriorGenotypes(data);
    _storeLocusNames();
};

void TVcfToLFMMPostGeno::_storePosteriorGenotypes(PopulationTools::TPopulationLikehoodLocus & data){
    auto * meanPostGenoForOneLocus = new double[_samples.numSamples()];
    for (size_t i = 0; i < _samples.numSamples(); i++){
        meanPostGenoForOneLocus[i] = _computePosteriorGenotype(data, i);
    }
    _genotypes.emplace_back(meanPostGenoForOneLocus);
};

double TVcfToLFMMPostGeno::_computePosteriorGenotype(PopulationTools::TPopulationLikehoodLocus & data, uint32_t i){
    if(data[i].isMissing()) {
        throw "Missing data at sample " + _samples.sampleName(i) + " and locus " + _reader.chr() + ":" + coretools::str::toString(_reader.position()) + "! LFMM2 does not accept missing genotypes, please impute your VCF file first.";
    }
    return data[i].meanPosteriorGenotype();
};

//------------------------------------------
// TVcfToPosFile
//------------------------------------------
TVcfToPosFile::TVcfToPosFile(TLog *Logfile) : TVcfConverter(Logfile) {
    // posFile is needed as input for STITCH
    // format:
    //   - tab-separated
    //   - no header
    //   - 4 columns: col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base
    //   - one row per SNP
    // Idea: one posfile per chromosome -> use option keepChromosomes
};

void TVcfToPosFile::_writeHeader(){
    //no header
    _posFile.noHeader(4);
};

void TVcfToPosFile::_initOutputFiles() {
    _posFile.open(_outname + ".pos");
};

void TVcfToPosFile::_writePosition(){
    _posFile << _reader.chr() << _reader.position();
};

void TVcfToPosFile::_writeRefAndAlt(){
    _posFile << _reader.refAllele() << _reader.altAllele();
};

void TVcfToPosFile::_writeData(PopulationTools::TPopulationLikehoodLocus &){
    _writePosition();
    _writeRefAndAlt();
    _posFile.endLine();
};

void TVcfToPosFile::vcfToPosFile(TParameters & Params){
    // read Vcf and write output
    readVcfAndWriteFile(Params);

    // clean up
    _posFile.close();
};

//------------------------------------------
// TVcfToGenotypeTruthSetFile
//------------------------------------------
TVcfToGenotypeTruthSetFile::TVcfToGenotypeTruthSetFile(TLog *Logfile) : TVcfConverter(Logfile) {
    // produces genotype truth set (genfile) for STITCH and bed files for samples
    // idea: first locus -> find 0-n samples that have depth higher than minDepth
    //                   -> write position of this locus into bed-files for these individuals
    //                   -> write _genotypes of these individuals to genfile; write _genotypes of all other individuals as NA to genfile
    //       next locus  -> is distance to previous locus more than x basepairs?
    //                   -> If yes: find 0-n samples that have depth higher than minDepth
    //                              -> write _genotypes of these individuals to genfile; write _genotypes of all other individuals as NA to genfile
    //                   -> Else: write _genotypes of all individuals as NA to genfile
    // format:
    //   - produces a BED file (3 columns, col 1 = chromosome, col 2 = start (1-based), col 3 = stop) for each individual
    //   - produces a genfile
    //      - tab-separated
    //      - header = sample names from vcf
    //      - one row per SNP
    //      - genotypes encoded as 0,1,2 or NA
    _bedFiles = nullptr;
    _positionPreviousLocus = 0;
    _minDistanceToPreviousLocus = 0;
    _numSamplesPerLocus = 0;
};

TVcfToGenotypeTruthSetFile::~TVcfToGenotypeTruthSetFile() {
    for (size_t s = 0; s < _samples.numSamples(); s++)
        delete _bedFiles[s];
    delete [] _bedFiles;
};

void TVcfToGenotypeTruthSetFile::_writeHeader(){
    //header string
    std::vector <std::string> header;
    for (size_t s = 0; s < _samples.numSamples(); s++){
        header.push_back(_samples.sampleName(s));
    }
    _genFile.writeHeader(header);
};

void TVcfToGenotypeTruthSetFile::_initOutputFiles() {
    _genFile.open(_outname + ".gen");
    // initialize bed files (we know how many samples there are)
    _bedFiles = new BAM::TBed * [_samples.numSamples()];
    for (size_t s = 0; s < _samples.numSamples(); s++) {
        _bedFiles[s] = new BAM::TBed;
    };

    // check if minNumSamplesWithData-Filter of TPopulationLikelihoods is zero
    // (even if all samples have missing data, we still need to keep the locus, because otherwise
    // positions in posfile and genfile do not match)
    if (_reader.getMinNumSamplesWithData() != 0)
        throw "Parameter 'minSamplesWithData' must be 0 for this task! Please specify 'minSamplesWithData=0'.";
};

void TVcfToGenotypeTruthSetFile::_mapIndividualsToDepth(std::vector<uint32_t> &samplesToKeep) {
    std::map< double, std::vector<uint32_t>, std::greater<> > depthVsSampleIndexMap; // use depth as key in map -> automatically sorted in descending order
    for (auto & s : samplesToKeep){
        // does depth already exist in map?
        double depth = _reader.depth(_samples, s);
        auto it = depthVsSampleIndexMap.find(depth);
        if (it == depthVsSampleIndexMap.end()){ // new depth
            std::vector<uint32_t> v = {s};
            depthVsSampleIndexMap[depth] = v;
        } else { // depth already exists
            it->second.push_back(s);
        }
    }
    _filterIndividualsWithHighestDepth(samplesToKeep, depthVsSampleIndexMap);
};

void TVcfToGenotypeTruthSetFile::_filterIndividualsWithHighestDepth(std::vector<uint32_t> & samplesToKeep, const std::map< double, std::vector<uint32_t>, std::greater<> > & depthVsSampleIndexMap){
    // only keep top _numSamplesPerLocus
    samplesToKeep.clear();
    int c = 0;
    for (auto & it : depthVsSampleIndexMap) { // one depth
        for (auto & it2 : it.second){ // one sample at this depth
            if (c >= _numSamplesPerLocus)
                return;
            samplesToKeep.push_back(it2);
            c++;
        }
    }
};

void TVcfToGenotypeTruthSetFile::_filterIndividuals(PopulationTools::TPopulationLikehoodLocus & data){
    std::vector<uint32_t> samplesToKeep;
    long distanceToPreviousLocus = _reader.position() - _positionPreviousLocus;
    if (distanceToPreviousLocus >= _minDistanceToPreviousLocus) { // check if distance is big enough
        // idea: TPopulationLikelihoods will filter on minDepth and set all samples with < minDepth as missing
        // here, we check how many individuals have > minDepth; we rank them and only keep _numSamplesPerLocus of them
        for (size_t s = 0; s < _samples.numSamples(); ++s) {
            if (!data[s].isMissing()) {
                samplesToKeep.push_back(s);
            }
        }
        if (samplesToKeep.empty()) // no individuals have > minDepth
            _writeToGenFile(samplesToKeep);
            // no need to write to bed file
        else if (samplesToKeep.size() <= static_cast<unsigned int>(_numSamplesPerLocus)) { // keep all of them
            _writeToGenFile(samplesToKeep);
            _storeInBedFile(samplesToKeep);
        }
        else { // rank according to depth
            _mapIndividualsToDepth(samplesToKeep);
            _writeToGenFile(samplesToKeep);
            _storeInBedFile(samplesToKeep);
        }
    } else {
        _writeToGenFile(samplesToKeep);
        // no need to write to bed file
    }
};

void TVcfToGenotypeTruthSetFile::_writeToGenFile(const std::vector<uint32_t> & samplesToKeep){
    for (size_t s = 0; s < _samples.numSamples(); s++) {
        // should we write true genotype of sample?
        auto it = std::find(samplesToKeep.begin(), samplesToKeep.end(), s);
        if (it != samplesToKeep.end()){
        	// sample found
        	_genFile << (std::string) _reader.genotype(_samples, s);

        } else {
            _genFile << "NA";
        }
    }
    _genFile.endLine();
};

void TVcfToGenotypeTruthSetFile::_storeInBedFile(const std::vector<uint32_t> & samplesToKeep){
    for (size_t s = 0; s < _samples.numSamples(); s++) {
        // should we write to bed of sample?
        auto it = std::find(samplesToKeep.begin(), samplesToKeep.end(), s);
        if (it != samplesToKeep.end()) // sample found
            _bedFiles[s]->add(_reader.chr(), _reader.position());
    }
    _positionPreviousLocus = _reader.position();
};

void TVcfToGenotypeTruthSetFile::_writeData(PopulationTools::TPopulationLikehoodLocus & data){
    if (_curChr != _reader.chr()){ // new chromosome
        _curChr = _reader.chr();
        _resetDistance();
    }
    _filterIndividuals(data);
};

void TVcfToGenotypeTruthSetFile::_resetDistance() {
    _positionPreviousLocus = -_minDistanceToPreviousLocus - 1; // for first locus on chromosome
};

void TVcfToGenotypeTruthSetFile::vcfToGenotypeTruthSetFile(TParameters & Params){
    // read params
    _minDistanceToPreviousLocus = Params.getParameterWithDefault<int>("minDistance", 100);
    _logfile->list("Will keep loci that have a minimal distance of ", _minDistanceToPreviousLocus, " to previous locus (parameter 'minDistance').");
    _resetDistance();
    _numSamplesPerLocus = Params.getParameterWithDefault<int>("numSamples", 5);
    _logfile->list("Will keep up to ", _numSamplesPerLocus, " individuals per locus (parameter 'numSamples').");

    // read Vcf and write output
    readVcfAndWriteFile(Params);

    // write bed files (one per sample)
    for(uint32_t s = 0; s < _samples.numSamples(); s++) {
        // check if sample name contains / (would be interpreted as path)
        std::string sample_name = _samples.sampleName(s);
        if (coretools::str::stringContains(sample_name, '/'))
            sample_name = coretools::stringReplace("/", "_", sample_name);
        _bedFiles[s]->write(_outname + "_" + sample_name + ".bed");
    };

    // clean up
    _genFile.close();
};

//------------------------------------------
// TVcfToVcf
//------------------------------------------
TVcfToVcf::TVcfToVcf(TLog *Logfile) : TVcfConverter(Logfile){
	//open output file
	std::string filename = _outname + (std::string) "_filtered.vcf.gz";
	gz::ogzstream out(filename.c_str());
	if(!out) throw "Failed to open outputfile '" + filename + "'!";
	_reader.setOutStream(out);
};

void TVcfToVcf::_writeData(PopulationTools::TPopulationLikehoodLocus &){
	_reader.writeVCFLine();
};

void TVcfToVcf::_writeHeader(){
	_reader.writeUnknownHeader();
};

}; //end namespace
