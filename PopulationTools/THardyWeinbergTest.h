/*
 * THardyWeinbergTest.h
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_THARDYWEINBERGTEST_H_
#define POPULATIONTOOLS_THARDYWEINBERGTEST_H_

#include "TParameters.h"
#include "TPopulationLikelihoods.h"
#include "TRandomGenerator.h"
#include "mathFunctions.h"

namespace PopulationTools{

using coretools::TParameters;
using coretools::TLog;
using coretools::TRandomGenerator;
using coretools::TOutputFile;

//------------------------------------------------
//THWHetProb
// Stores the probabilities P(N_AB = n|N, n_A, n_B),
// i.e the probabilities of observing n hets among N individuals if the allele frequencies are n_A + n_B = N
//------------------------------------------------
class THWHetProb{
private:
	std::vector<double> _probs;
	uint32_t _maxNumHetPlusOne;
	bool _onlyOdd, _onlyEven;

public:
	THWHetProb();
	THWHetProb(uint32_t N, uint32_t n_A);

	void clear();
	THWHetProb& operator=(const THWHetProb & other);
	void extend(const THWHetProb & other);
	bool onlyOdd() const { return _onlyOdd; };
	bool onlyEven() const { return _onlyEven; };
	uint32_t maxNumHet() const { return _maxNumHetPlusOne - 1; };
	const double& operator[](uint32_t i) const { return _probs[i]; };
	double sum(uint32_t upTo);
	void print();
};

//------------------------------------------------
//THWHetProbDB
//------------------------------------------------
class THWHetProbVector{
private:
	uint32_t _N;
	std::map<uint32_t, THWHetProb> _probs;

public:
	THWHetProbVector(uint32_t N);
	const THWHetProb& getProbs(uint32_t n_A);
};

class THWProbDB{
private:
	std::map<uint32_t, THWHetProbVector> _probs;

public:
	const THWHetProb& getProbs(uint32_t N, uint32_t n_A);
};

//------------------------------------------------
//THWGenotypes
//------------------------------------------------
class THWGenotypes{
private:
	std::array<uint32_t,3> _genoCounts;

public:
	THWGenotypes(){ clear(); };

	void clear();
	void add(const genometools::BiallelicGenotype & genotype);
	uint32_t N() const;
	uint32_t MAF() const;
	uint32_t n_A() const;
	uint32_t n_AB() const;
};

//------------------------------------------------
//THWPopulations
//------------------------------------------------
class THWPopulations{
private:
	std::vector<THWGenotypes> _populations;
	THWProbDB _probDB;

public:
	THWPopulations(){};
	THWPopulations(uint16_t numPops);
	void clear();
	void resize(uint16_t numPops);
	void add(uint16_t pop, const genometools::BiallelicGenotype & genotyp);
	void addToHeader(std::vector<std::string> & header);
	void runTest(TOutputFile & out);
};

//------------------------------------------------
//THardyWeinbergTest
//------------------------------------------------
class THardyWeinbergTest{
private:
	TLog* _logfile;
	TRandomGenerator* _randonGenerator;
	std::string _outname;

	//vcf-file
	std::string _vcfFilename;
	VCF::TVcfFileSingleLine _vcfFile;
	bool _limitLines;
	uint64_t _maxNumLines;

	//samples
	TPopulationSamples _samples;

	//genotype data
	THWPopulations _populations;

	void _openVCF(TParameters & Parameters);
	void _closeVCF();

public:
	THardyWeinbergTest(TParameters & Parameters, TLog* logfile, TRandomGenerator* RandonGenerator);
	void testForHardyWeinberg();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_testHardyWeinberg:public coretools::TTask{
public:
	TTask_testHardyWeinberg(){ _explanation = "Testing for Hardy-Weinberg equilibrium across multiple populations"; };

	void run(TParameters & Parameters, TLog* Logfile){
		THardyWeinbergTest HW_test(Parameters, Logfile, _randomGenerator);
		HW_test.testForHardyWeinberg();
	};
};

}; //end namespace


#endif /* POPULATIONTOOLS_THARDYWEINBERGTEST_H_ */
