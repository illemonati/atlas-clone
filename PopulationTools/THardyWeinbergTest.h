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

namespace PopulationTools{

//------------------------------------------------
//THWGenotypes
//------------------------------------------------
class THWGenotypes{
private:
	std::array<uint32_t,3> _genoCounts;

public:
	THWGenotypes(){ clear(); };

	void clear(){ _genoCounts.fill(0); };
	void add(uint8_t genotype){ ++_genoCounts[genotype]; };
};

//------------------------------------------------
//THWPopulations
//------------------------------------------------
class THWPopulations{
private:
	std::vector<THWGenotypes> _populations;

public:
	THWPopulations(uint16_t numPops){
		_populations.resize(numPops);
	};

	void clear(){
		for(auto& p : _populations){
			p.clear();
		}
	};

	void add(uint16_t pop, uint8_t genotyp){
		_populations[pop].add(genotyp);
	};
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
