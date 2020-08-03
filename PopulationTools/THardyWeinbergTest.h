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
#include "TQualityMap.h"
#include "TRandomGenerator.h"
#include "mathFunctions.h"

namespace PopulationTools{



//------------------------------------------------
//THWHetProbDB
// Stores the probabilities P(N_AB = n|N, n_A, n_B),
// i.e the probabilities of observing n hets among N individuals if the allele frequencies are n_A + n_B = N
//------------------------------------------------
class THWHetProb{
private:
	uint32_t _maxNumHet;

	bool _odd; //indicates whether the numbe ro fpossible hets are even or odd
	std::vector<double> _probs;

public:
	THWHetProb(const uint32_t & N, const uint32_t & n_A);

	const uint32_t& maxNumHet() const { return _maxNumHet; };
	const bool& odd() const { return _odd; };
	const double& operator[](const uint32_t & i) const { return _probs[i]; };
};


class THWHetProbVector{
private:
	uint32_t _N;
	std::map<uint32_t, THWHetProb> _probs;

public:
	THWHetProbVector(const uint32_t & N);
	const THWHetProb& getProbs(const uint32_t & n_A);
};

class THWProbDB{
private:
	std::map<uint32_t, THWHetProbVector> _probs;

public:
	const THWHetProb& getProbs(const uint32_t & N, const uint32_t & n_A);
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
	void add(const uint8_t & genotype);
	uint32_t N();
	uint32_t n_A();
};

//------------------------------------------------
//THWProbsMultiPop
//------------------------------------------------
class THWProbsMultiPop{
private:
	std::vector<double> _probs;
	uint32_t _maxNumHet;
	uint32_t _curMaxNumHetPlusOne;

public:
	THWProbsMultiPop();

	void initialize(const uint32_t & maxNumHet, const THWHetProb & probs);
	void extend(const THWHetProb & probs);
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
	THWPopulations(const uint16_t & numPops);
	void clear();
	void resize(const uint16_t & numPops);
	void add(const uint16_t & pop, const uint8_t & genotyp);

	void runTest();
};

//------------------------------------------------
//THardyWeinbergTest
//------------------------------------------------
class THardyWeinbergTest{
private:
	TLog* _logfile;
	TRandomGenerator* _randonGenerator;

	//vcf-file
	std::string _vcfFilename;
	VCF::TVcfFileSingleLine _vcfFile;

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
class TTask_testHardyWeinberg:public TTask{
public:
	TTask_testHardyWeinberg(){ _explanation = "Testing for HardyWeinberg (Multi-population)"; };

	void run(TParameters & Parameters, TLog* Logfile){
		THardyWeinbergTest HW_test(Parameters, Logfile, _randomGenerator);
		HW_test.testForHardyWeinberg();
	};
};

}; //end namespace


#endif /* POPULATIONTOOLS_THARDYWEINBERGTEST_H_ */
