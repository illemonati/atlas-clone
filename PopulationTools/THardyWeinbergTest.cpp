/*
 * THardyWeinbergTest.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#include "THardyWeinbergTest.h"

namespace PopulationTools{


//------------------------------------------------
//THWHetProb
//------------------------------------------------
THWHetProb::THWHetProb(){
	clear();
};

THWHetProb::THWHetProb(const uint32_t & N, const uint32_t & n_A){
	//max num het
	_maxNumHetPlusOne = n_A + 1;
	if(_maxNumHetPlusOne > N / 2.0){
		_maxNumHetPlusOne = 2*N - n_A + 1;
	}
	_onlyOdd = n_A % 2;
	_onlyEven = !_onlyOdd;

	//tmp variables
	uint32_t n_B = 2*N - n_A;
	double logNFactorial = TFactorial::factorialLog(N);
	double logTwoNFactorial = TFactorial::factorialLog(2*N);

	//fill probabilities
	_probs.resize(_maxNumHetPlusOne);
	bool calcProb = _onlyEven;
	for(uint32_t n_AB = 0; n_AB<_maxNumHetPlusOne; ++n_AB){
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

void THWHetProb::clear(){
	_maxNumHetPlusOne = 0;
	_onlyOdd = false;
	_onlyEven = false;
};

void THWHetProb::extend(const THWHetProb & other){
	uint32_t otherMaxNumHetPlusOne = other.maxNumHet() + 1;

	//create new storage
	uint32_t newMaxHetPlusOne = _maxNumHetPlusOne + other.maxNumHet();
	std::vector<double> newProbs(newMaxHetPlusOne);
	std::fill(newProbs.begin(), newProbs.end(), 0.0);

	//combine: depends on odd even in both
	if(_onlyOdd){
		if(other.onlyOdd()){
			for(uint32_t i = 1; i<otherMaxNumHetPlusOne; i+=2){
				for(uint32_t j = 1; j<_maxNumHetPlusOne; j+=2){
					newProbs[j+i] += _probs[j] * other[i];
				}
			}
		} else if(other.onlyEven()){
			for(uint32_t i = 0; i<otherMaxNumHetPlusOne; i+=2){
				for(uint32_t j = 1; j<_maxNumHetPlusOne; j+=2){
					newProbs[j+i] += _probs[j] * other[i];
				}
			}
		} else {
			for(uint32_t i = 0; i<otherMaxNumHetPlusOne; ++i){
				for(uint32_t j = 1; j<_maxNumHetPlusOne; j+=2){
					newProbs[j+i] += _probs[j] * other[i];
				}
			}
		}
	} else if(_onlyEven){
		if(other.onlyOdd()){
			for(uint32_t i = 1; i<otherMaxNumHetPlusOne; i+=2){
				for(uint32_t j = 0; j<_maxNumHetPlusOne; j+=2){
					newProbs[j+i] += _probs[j] * other[i];
				}
			}
		} else if(other.onlyEven()){
			for(uint32_t i = 0; i<otherMaxNumHetPlusOne; i+=2){
				for(uint32_t j = 0; j<_maxNumHetPlusOne; j+=2){
					newProbs[j+i] += _probs[j] * other[i];
				}
			}
		} else {
			for(uint32_t i = 0; i<otherMaxNumHetPlusOne; ++i){
				for(uint32_t j = 0; j<_maxNumHetPlusOne; j+=2){
					newProbs[j+i] += _probs[j] * other[i];
				}
			}
		}
	} else {
		if(other.onlyOdd()){
			for(uint32_t i = 1; i<otherMaxNumHetPlusOne; i+=2){
				for(uint32_t j = 0; j<_maxNumHetPlusOne; ++j){
					newProbs[j+i] += _probs[j] * other[i];
				}
			}
		} else if(other.onlyEven()){
			for(uint32_t i = 0; i<otherMaxNumHetPlusOne; i+=2){
				for(uint32_t j = 0; j<_maxNumHetPlusOne; ++j){
					newProbs[j+i] += _probs[j] * other[i];
				}
			}
		} else {
			for(uint32_t i = 0; i<otherMaxNumHetPlusOne; ++i){
				for(uint32_t j = 0; j<_maxNumHetPlusOne; ++j){
					newProbs[j+i] += _probs[j] * other[i];
				}
			}
		}
	}

	//update storage
	std::swap(_probs, newProbs);
	_maxNumHetPlusOne = newMaxHetPlusOne;
	_onlyOdd = _onlyOdd && other.onlyOdd();
	_onlyEven = _onlyEven && other.onlyEven();

	//probabilities should sum to one. Make sure they do.
	double sum = 0.0;
	for(auto& i : _probs){
		sum += i;
	}

	std::cout << "SUM = " << sum << std::endl;

	for(auto& i : _probs){
		i /= sum;
	}
};

double THWHetProb::sum(const uint32_t & upTo){
	double sum = 0.0;
	for(uint32_t i = 0; i<=upTo; ++i){
		sum += _probs[i];
	}
	return sum;
};

//------------------------------------------------
//THWHetProbDB
//------------------------------------------------


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

uint32_t THWGenotypes::N() const{
	return _genoCounts[0] + _genoCounts[1] + _genoCounts[2];
};

uint32_t THWGenotypes::n_A() const{
	if(_genoCounts[0] < _genoCounts[2]){
		return 2*_genoCounts[0] + _genoCounts[1];
	} else {
		return 2*_genoCounts[2] + _genoCounts[1];
	}
};

uint32_t THWGenotypes::n_AB() const{
	return _genoCounts[1];
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

void THWPopulations::addToHeader(std::vector<std::string> & header){
	header.push_back("N");
	header.push_back("n_AB");
	header.push_back("pLess");
	header.push_back("pMore");
};

void THWPopulations::runTest(TOutputFile & out){
	// do dynamic programming to get probabilities

	//get max num het across populations
	uint32_t maxNumHet = 0;
	for(auto& p : _populations){
		maxNumHet += p.n_A();
	}

	//loop over populations to calculate P(N_AB = n | N_1, N_2, ..., n_A1.n_A2, ...) using dynamic programming
	//first do those with even and odd num hets separately
	THWHetProb probs, probsEven;
	uint32_t obsNumHet = 0;
	uint32_t obsN = 0;

	//loop over populations
	for(auto& p : _populations){
		obsNumHet += p.n_AB();
		obsN += p.N();
		const THWHetProb& popProbs = _probDB.getProbs(p.N(), p.n_A());
		if(popProbs.onlyOdd()){
			probs.extend(popProbs);
		} else {
			probsEven.extend(popProbs);
		}
	}

	//combine even and odd
	probs.extend(probsEven);

	//calculate p-values
	double pLess = probs.sum(obsNumHet);
	out << obsN << obsNumHet << pLess << 1.0 - pLess + probs[obsNumHet];
};

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
	//open output file
	std::vector<std::string> header = {"chr", "position"};
	_populations.addToHeader(header);
	TOutputFile out;

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
