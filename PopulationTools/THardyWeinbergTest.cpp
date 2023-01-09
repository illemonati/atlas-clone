/*
 * THardyWeinbergTest.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */


#include "THardyWeinbergTest.h"

#include <math.h>
#include <algorithm>
#include <exception>
#include <iostream>
#include <memory>
#include <tuple>
#include <utility>

#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/TTimer.h"
#include "genometools/VCF/TVcfParser.h"

namespace PopulationTools{

using coretools::TOutputFile;
using coretools::instances::logfile;
using coretools::instances::parameters;

//------------------------------------------------
//THWHetProb
//------------------------------------------------
THWHetProb::THWHetProb(){
	clear();
};

THWHetProb::THWHetProb(uint32_t numInd_N, uint32_t alleleFreq_n_A){
	using namespace coretools::TFactorial;
	if(alleleFreq_n_A > 0 && alleleFreq_n_A < 2*numInd_N){
		//max num het
		_maxNumHetPlusOne = std::min(alleleFreq_n_A, 2*numInd_N - alleleFreq_n_A) + 1;
		_onlyOdd = alleleFreq_n_A % 2;
		_onlyEven = !_onlyOdd;

		//tmp variables
		uint32_t n_B = 2*numInd_N - alleleFreq_n_A;
		double logNFactorial = factorialLog(numInd_N);
		double logTwoNFactorial = factorialLog(2*numInd_N);

		//fill probabilities
		_probs.resize(_maxNumHetPlusOne);

		if(_maxNumHetPlusOne == 1){
			_probs[0] = 1.0;
		} else {
			bool calcProb = _onlyEven;
			for(uint32_t n_AB = 0; n_AB<_maxNumHetPlusOne; ++n_AB){
				if(calcProb){
					uint32_t n_AA = (alleleFreq_n_A - n_AB) / 2;
					uint32_t n_BB = numInd_N - n_AB - n_AA;

					_probs[n_AB] = exp(n_AB * 0.6931472
								 + logNFactorial - factorialLog(n_AA) - factorialLog(n_AB) - factorialLog(n_BB)
								 + factorialLog(alleleFreq_n_A) + factorialLog(n_B) - logTwoNFactorial);
				} else {
					_probs[n_AB] = 0.0;
				}
				calcProb = !calcProb;
			}
		}
	} else {
		clear();
	}
};

void THWHetProb::clear(){
	_maxNumHetPlusOne = 1;
	_onlyOdd = false;
	_onlyEven = false;
	_probs.resize(_maxNumHetPlusOne);
	_probs[0] = 1.0;
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

	for(auto& i : _probs){
		i /= sum;
	}
};

double THWHetProb::sum(uint32_t upTo){
	double sum = 0.0;
	for(uint32_t i = 0; i<=upTo; ++i){
		sum += _probs[i];
	}
	return sum;
};

void THWHetProb::print(){
	for(uint32_t i=0; i<_maxNumHetPlusOne; ++i){
		std::cout << "P(n_AB = " << i << ") = " << _probs[i] << std::endl;
	}
};

//------------------------------------------------
//THWHetProbDB
//------------------------------------------------
THWHetProbVector::THWHetProbVector(uint32_t N){
	_N = N;
};

const THWHetProb& THWHetProbVector::getProbs(uint32_t n_A){
	//check if these probs already exist
	auto it = _probs.find(n_A);
	if(it == _probs.end()){
		//create
		it = _probs.emplace(std::piecewise_construct, std::make_tuple(n_A), std::make_tuple(_N, n_A)).first;
	}

	return it->second;
};

const THWHetProb& THWProbDB::getProbs(uint32_t N, uint32_t n_A){
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
	_genoCounts[0] = 0;
	_genoCounts[1] = 0;
	_genoCounts[2] = 0;
};

void THWGenotypes::add(const genometools::BiallelicGenotype & genotype){
	using namespace genometools;
	//ToDo: extend test to haploid
	if(!isMissing(genotype) && isDiploid(genotype)){
		++_genoCounts[coretools::index(genotype)];
	}
};

uint32_t THWGenotypes::N() const{
	return _genoCounts[0] + _genoCounts[1] + _genoCounts[2];
};

uint32_t THWGenotypes::MAF() const{
	if(_genoCounts[0] < _genoCounts[2]){
		return 2*_genoCounts[0] + _genoCounts[1];
	} else {
		return 2*_genoCounts[2] + _genoCounts[1];
	}
};

uint32_t THWGenotypes::n_A() const{
	return 2*_genoCounts[0] + _genoCounts[1];
};

uint32_t THWGenotypes::n_AB() const{
	return _genoCounts[1];
};

//------------------------------------------------
//THWPopulations
//------------------------------------------------
THWPopulations::THWPopulations(uint16_t numPops){
	resize(numPops);
};

void THWPopulations::resize(uint16_t numPops){
	_populations.resize(numPops);
};

void THWPopulations::clear(){
	for(auto& p : _populations){
		p.clear();
	}
};

void THWPopulations::add(uint16_t pop, const genometools::BiallelicGenotype & genotyp){
	_populations[pop].add(genotyp);
};

void THWPopulations::addToHeader(std::vector<std::string> & header){
	header.push_back("N");
	header.push_back("n_A");
	header.push_back("MAF");
	header.push_back("n_AB");
	header.push_back("pLess");
	header.push_back("pMore");
};

void THWPopulations::runTest(TOutputFile & out){
	// do dynamic programming to get probabilities

	//loop over populations to calculate P(N_AB = n | N_1, N_2, ..., n_A1.n_A2, ...) using dynamic programming
	//first do those with even and odd num hets separately
	THWHetProb probs, probsEven;
	uint32_t obsNumHet = 0;
	uint32_t obsNumA = 0;
	uint32_t obsN = 0;

	//loop over populations
	for(auto& p : _populations){
		obsNumHet += p.n_AB();
		obsNumA += p.n_A();
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
	double MAF = std::min(obsNumA, 2*obsN - obsNumA) / (double) (obsN * 2.0);
	out.writeln(obsN, obsNumA, MAF, obsNumHet, pLess, 1.0 - pLess + probs[obsNumHet]);
};

//------------------------------------------------
//THardyWeinbergTest
//------------------------------------------------
THardyWeinbergTest::THardyWeinbergTest(){

	//read samples
	if(parameters().parameterExists("samples")){
		_samples.readSamples(parameters().getParameter<std::string>("samples"));
	}

	//open VCF
	_vcfFilename = parameters().getParameter<std::string>("vcf");
	logfile().list("Reading vcf from file '" + _vcfFilename + "'.");
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

	//resize populations
	_populations.resize(_samples.numPopulations());

	//get output name
	std::string tmp = coretools::str::extractBeforeLast(_vcfFilename, ".vcf");
	_outname = parameters().getParameterWithDefault<std::string>("out", tmp);

	//limit lines?
	_limitLines = parameters().parameterExists("limitLines");
	if(_limitLines){
		_maxNumLines = parameters().getParameter<long>("limitLines");
	} else {
		_maxNumLines = 0;
	}
};

void THardyWeinbergTest::testForHardyWeinberg(){
	//open output file
	std::string filename = _outname + "_HWTest.txt.gz";
	logfile().list("Writing HW test results to file '" + filename + "'.");
	std::vector<std::string> header = {"chr", "position"};
	_populations.addToHeader(header);
	TOutputFile out(filename, header);

	//progress
	coretools::TTimer timer;
	uint64_t lineCounter = 0;
	uint64_t numFiltered = 0;

	//traverse VCF
	logfile().startIndent("Traversing VCF file:");
	while(_vcfFile.next()){
		++lineCounter;

		//exclude multiallelic
		if(_vcfFile.getNumAlleles() == 2){
			//reset counts
			_populations.clear();

			//write position
			out << _vcfFile.chr() << _vcfFile.position();

			//add data at current line
			for(uint32_t s = 0; s<_samples.numSamples(); ++s){
				uint32_t vcfIndex = _samples.sampleIndexInVCF(s);
				if(!_vcfFile.sampleIsMissing(vcfIndex)){
					_populations.add(_samples.populationIndex(s), _vcfFile.sampleBiallelicGenotype(vcfIndex));
				}
			}

			//perform test
			_populations.runTest(out);

			//progress / limit lines
			if(lineCounter % 10000 == 0){
				logfile().list("Parsed ", lineCounter, " lines in " + timer.formattedTime());
			}
			if(_limitLines && lineCounter == _maxNumLines){
				break;
			}
		} else {
			++numFiltered;
		}
	}
	logfile().list("Reached end of VCf file.");
	logfile().conclude("Parsed ", lineCounter, " lines in ", timer.formattedTime());
	logfile().conclude("Ignored ", numFiltered, " sites that were not bi-allelic SNPs.");
	logfile().endIndent();

	//close file
	out.close();
};


}; //end namesapce
