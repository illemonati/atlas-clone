/*
 * THardyWeinbergTest.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#include "THardyWeinbergTest.h"

namespace PopulationTools{


//------------------------------------------------
//THWHetProbDB
//------------------------------------------------
THWHetProb::THWHetProb(const uint32_t & N, const uint32_t & n_A){
	//max num het
	_maxNumHet = n_A;
	if(_maxNumHet > N / 2.0){
		_maxNumHet = 2*N - n_A;
	}
	_odd = n_A % 2;

	//tmp variables
	uint32_t n_B = 2*N - n_A;
	double logNFactorial = TFactorial::factorialLog(N);
	double logTwoNFactorial = TFactorial::factorialLog(2*N);

	//fill probabilities
	_probs.resize(_maxNumHet + 1);
	bool calcProb = !_odd;
	for(uint32_t n_AB = 0; n_AB<=_maxNumHet; ++n_AB){
		if(calcProb){
			uint32_t n_AA = (n_A - n_AB) / 2;
			uint32_t n_BB = N - n_AB - n_AA;

			_probs[n_AB] = n_AB * 0.6931472
					     + logNFactorial - TFactorial::factorialLog(n_AA) - TFactorial::factorialLog(n_AB) - TFactorial::factorialLog(n_BB)
					     + TFactorial::factorialLog(n_A) + TFactorial::factorialLog(n_B) - logTwoNFactorial;
		} else {
			_probs[n_AB] = 0.0;
		}
		calcProb = !calcProb;
	}
};

THWHetProbVector::THWHetProbVector(const uint32_t & N){
	_N = N;
};

const THWHetProb& THWHetProbVector::getProbs(const uint32_t & n_A){
	//check if these probs already exist
	auto it = _probs.find(n_A);
	if(it == _probs.end()){
		//create
		it = _probs.emplace(std::piecewise_construct, std::make_tuple(n_A), std::make_tuple(_N, n_A)).first;
	}
	return it->second;
};

const THWHetProb& THWProbDB::getProbs(const uint32_t & N, const uint32_t & n_A){
	//check if prob vector already exists
	auto it = _probs.find(N);
	if(it == _probs.end()){
		//create
		it = _probs.emplace(N, N).first;
	}

	return it->second.getProbs(n_A);
};

//------------------------------------------------
//THWGenotypes
//------------------------------------------------
void THWGenotypes::clear(){
	_genoCounts.fill(0);
};

void THWGenotypes::add(const uint8_t & genotype){
	++_genoCounts[genotype];
};

uint32_t THWGenotypes::N(){
	return _genoCounts[0] + _genoCounts[1] + _genoCounts[2];
};

uint32_t THWGenotypes::n_A(){
	if(_genoCounts[0] < _genoCounts[2]){
		return 2*_genoCounts[0] + _genoCounts[1];
	} else {
		return 2*_genoCounts[2] + _genoCounts[1];
	}
};


//------------------------------------------------
//THWProbsMultiPop
//------------------------------------------------
THWProbsMultiPop::THWProbsMultiPop(){
	_maxNumHet = 0;
	_curMaxNumHetPlusOne = 0;
};

void THWProbsMultiPop::initialize(const uint32_t & maxNumHet, const THWHetProb & probs){
	_maxNumHet = maxNumHet;
	_probs.resize(_maxNumHet);
	std::fill(_probs.begin(), _probs.end(), 0.0);

	//copy from probs
	for(uint32_t i = probs.odd(); i<probs.maxNumHet(); i+=2){
		_probs[i] = probs[i];
	};
	_curMaxNumHetPlusOne = probs.maxNumHet() + 1;
};

void THWProbsMultiPop::extend(const THWHetProb & probs){
	//calculate P(N_AB = n | N_1, N_2, ..., n_A1.n_A2, ...) given these probs for one population less
	for(uint32_t i = probs.odd(); i<probs.maxNumHet(); i+=2){
		for(uint32_t j = 0; j<_curMaxNumHetPlusOne; ++i){

		}
	}
};

//------------------------------------------------
//THWPopulations
//------------------------------------------------
THWPopulations::THWPopulations(const uint16_t & numPops){
	resize(numPops);
};

void THWPopulations::resize(const uint16_t & numPops){
	_populations.resize(numPops);
};

void THWPopulations::clear(){
	for(auto& p : _populations){
		p.clear();
	}
};

void THWPopulations::add(const uint16_t & pop, const uint8_t & genotyp){
	_populations[pop].add(genotyp);
};

void THWPopulations::runTest(){
	// do dynamic programming to get probabilities

	//get max num het across populations
	uint32_t maxNumHet = 0;
	for(auto& p : _populations){
		maxNumHet += p.n_A();
	}

	//prepare storage
	std::vector<double> probs(maxNumHet+1, 0.0);

	//loop over populations to calculate P(N_AB = n | N_1, N_2, ..., n_A1.n_A2, ...) using dynamic programming
	//initialize with first population
	auto p = _populations.begin();
	const THWHetProb& popProbs = _probDB.getProbs(p->N(), p->n_A());

}

//------------------------------------------------
//THardyWeinbergTest
//------------------------------------------------
THardyWeinbergTest::THardyWeinbergTest(TParameters & Parameters, TLog* logfile, TRandomGenerator* RandonGenerator){
	_logfile = logfile;
	_randonGenerator = RandonGenerator;

	//read samples
	if(Parameters.parameterExists("samples")){
		_samples.readSamples(Parameters.getParameterString("samples"), logfile);
	}

	//open VCF
	_vcfFilename = Parameters.getParameterString("vcf");
	_logfile->startIndent("Reading vcf from file '" + _vcfFilename + "'.");
	_vcfFile.openStream(_vcfFilename);

	//enable parsers
	_vcfFile.enablePositionParsing();
	_vcfFile.enableVariantParsing();
	_vcfFile.enableVariantQualityParsing();
	_vcfFile.enableFormatParsing();
	_vcfFile.enableSampleParsing();

	//Match samples
	if(_samples.hasSamples()){
		_samples.fillVCFOrder(_vcfFile.parser.samples);
	} else {
		_samples.readSamplesFromVCFNames(_vcfFile.parser.samples);
	}
};

void THardyWeinbergTest::testForHardyWeinberg(){
	//traverse VCF
	while(_vcfFile.next()){
		//reset counts
		_populations.clear();

		//add data at current line
		for(uint32_t s = 0; s<_samples.numSamples(); ++s){
			uint32_t vcfIndex = _samples.VCF_order(s);
			if(!_vcfFile.sampleIsMissing(vcfIndex)){
				_populations.add(_samples.popIndex(s), _vcfFile.sampleGenotype(vcfIndex));
			}
		}

		//perform test

	}
};


}; //end namesapce
