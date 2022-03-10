/*
 * TPopulationLikelihoods.cpp
 *
 *  Created on: Dec 7, 2018
 *      Author: phaentu
 */

#include "TPopulationLikelihoods.h"
#include "TVCFReader.h"

namespace PopulationTools{

using coretools::str::toString;

////////////////////////////////////////////////////////////////////////////////////////////////
// TVcfFilter                                                                                 //
////////////////////////////////////////////////////////////////////////////////////////////////
/*
TVcfFilter::TVcfFilter(TParameters & Parameters, TLog* Logfile){
	vcfRead = false;
	_numLoci = 0;
	logfile = Logfile;
}

void TVcfFilter::filterVCF(TParameters & Parameters){
	if(vcfRead)
		throw "VCF already read!";

	//create reader
	bool saveAlleleFrequencies = false;
	TPopulationLikelihoodReader reader(Parameters, logfile, saveAlleleFrequencies);

	// open vcf file
	vcfFilename = Parameters.getParameter<std::string>("vcf");
	logfile->startIndent("Filtering VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename, logfile);

	//Match samples
	if(samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	 else
		 samples.readSamplesFromVCFNames(reader.getSampleVCFNames());


	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	uint8_t* curLocus = new uint8_t[samples.numSamples() * 3];
	bool* sampleIsMissing = new bool[samples.numSamples()];

	//output file name
	std::string tmp = extractBeforeLast(vcfFilename, ".vcf");
	std::string outputName = Parameters.getParameterWithDefault<std::string>("out", tmp) + "_filtered.vcf.gz";

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    _numLoci = 0;
    reader.filterVCF(curLocus, sampleIsMissing, samples, logfile, outputName);

    //clean up
	vcfRead = true;
	delete[] curLocus;

    //report final status
	logfile->endIndent();
	reader.concludeFilters(logfile);
	if(reader.numAcceptedLoci() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};
*/

////////////////////////////////////////////////////////////////////////////////////////////////
// TPopulationLikelihoods                                                                     //
////////////////////////////////////////////////////////////////////////////////////////////////
TPopulationLikelihoods::TPopulationLikelihoods(){
	init();
};


TPopulationLikelihoods::TPopulationLikelihoods(coretools::TParameters & Parameters, coretools::TLog* Logfile){
	init();
	readData(Parameters, Logfile);
};

TPopulationLikelihoods::~TPopulationLikelihoods(){
	clean();
};

void TPopulationLikelihoods::init(){
	vcfRead = false;
	saveAlleleFrequencies = false;
	saveTrueAlleleFrequencies = false;
    _positions = std::make_shared<genometools::TPositionsRaw>();
};

void TPopulationLikelihoods::clean(){
	for(TSampleLikelihoods* it : data){
		delete it;
	}
	data.clear();
	vcfRead = false;
};

void TPopulationLikelihoods::readData(coretools::TParameters & Parameters, coretools::TLog* Logfile){
	//check if we limit samples
	if(Parameters.parameterExists("samples"))
		samples.readSamples(Parameters.getParameter<std::string>("samples"), Logfile);

	//read Data
	readDataFromVCF(Parameters, Logfile);

	//make sure loop goes across all samples
	useAllSamples();
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read data from VCF-file                                                                      //
//////////////////////////////////////////////////////////////////////////////////////////////////
void TPopulationLikelihoods::readDataFromVCF(coretools::TParameters & Parameters, coretools::TLog* logfile){
	if(vcfRead)
		throw "VCF already read!";

	//create reader
	genometools::TPopulationLikelihoodReaderLocus reader(Parameters, logfile, saveAlleleFrequencies);
	if(saveAlleleFrequencies){
		reader.doEstimateGenotypeFrequencies();
	}
	if(saveTrueAlleleFrequencies)
		reader.doSaveTrueAlleleFrequencies();

	// open vcf file
	vcfFilename = Parameters.getParameter<std::string>("vcf");
	logfile->startIndent("Reading genotype likelihoods from VCF file '" + vcfFilename + "':");
	reader.openVCF(vcfFilename);

	//Match samples
	if (samples.hasSamples())
		samples.fillVCFOrder(reader.getSampleVCFNames());
	else
		samples.readSamplesFromVCFNames(reader.getSampleVCFNames());

	// initialize variables for vcf-file
	struct timeval start; gettimeofday(&start, NULL);
	data.emplace_back(new TSampleLikelihoods[samples.numSamples()]);

    //run through VCF file
    logfile->startIndent("Parsing VCF file:");
    while(reader.readDataFromVCF(data.back(), samples)){

        //store chromosome and position
        _positions->add(reader.position(), reader.chr());

		//store allele frequencies
		if(saveAlleleFrequencies)
			alleleFrequencies.emplace_back(reader.allelFrequency());
		if(saveTrueAlleleFrequencies)
			trueAlleleFrequencies.emplace_back(reader.trueAlleleFrequency());

		//update for next
		data.emplace_back(new TSampleLikelihoods[samples.numSamples()]);
    }

    //remove last, empty container
    delete data.back();
    data.pop_back();
    _positions->finalizeFilling();

    //clean up
	vcfRead = true;

    //report final status
	logfile->endIndent();
	reader.concludeFilters();
    if (_positions->size() < 1)
		throw "No usable loci in VCF file '" + vcfFilename + "'!";
	logfile->endIndent();
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// get-functions                                                                                //
//////////////////////////////////////////////////////////////////////////////////////////////////

size_t TPopulationLikelihoods::getNumIndividuals() const {
	return samples.numSamples();
};

size_t TPopulationLikelihoods::getNumLoci() const {
	return _positions->size();
};

TSampleLikelihoods* TPopulationLikelihoods::getDataAtLocus(long index){
	return data[index];
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// functions to loop over data                                                                  //
//////////////////////////////////////////////////////////////////////////////////////////////////
void TPopulationLikelihoods::useAllSamples(){
	individualStartIndex = 0;
	usedSampleSize = samples.numSamples();
};

void TPopulationLikelihoods::limitToSinglePopulation(int population){
	individualStartIndex = samples.startIndex(population);
	usedSampleSize = samples.numSamplesInPop(population);
};

void TPopulationLikelihoods::begin(){
	curLocusIndex = 0;
};

void TPopulationLikelihoods::beginAll(){
	useAllSamples();
	begin();
};

void TPopulationLikelihoods::beginOnePop(int population){
	limitToSinglePopulation(population);
	begin();
};

bool TPopulationLikelihoods::end(){
	return curLocusIndex == _positions->size();
};

void TPopulationLikelihoods::next(){
	++curLocusIndex;
};

TSampleLikelihoods* TPopulationLikelihoods::curData(){
	return data[curLocusIndex];
};

std::string TPopulationLikelihoods::curSampleName(int index){
	return samples.sampleName(individualStartIndex + index);
};

int TPopulationLikelihoods::curSampleSize(){
	return usedSampleSize;
};

std::string TPopulationLikelihoods::curChr(){
	return _positions->getJunkName(curLocusIndex);
};

long TPopulationLikelihoods::curPosition(){
	return _positions->getPosition(curLocusIndex);
};

const std::shared_ptr<genometools::TPositionsRaw> &TPopulationLikelihoods::positions() const {
    return _positions;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// print data                                                                                   //
//////////////////////////////////////////////////////////////////////////////////////////////////

void TPopulationLikelihoods::print(){
	begin();

	//write header
	std::cout << "Chr\tPos";
	for(int s=0; s<curSampleSize(); ++s){
		std::cout << "\t" << curSampleName(s);
	}
	std::cout << std::endl;

	//print data
	for(; !end(); next()){
		//print chromosome and position
		std::cout << curChr() << "\t" << curPosition();

		//print data
		TSampleLikelihoods* data = curData();
		for(int s=0; s<curSampleSize(); ++s){
			std::cout << "\t";
			data[s].print();
		}

		std::cout << std::endl;
	}
};

}; //end namespace
