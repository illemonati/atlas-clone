/*
 * THardyWeinbergTest.h
 *
 *  Created on: Jul 31, 2020
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_THARDYWEINBERGTEST_H_
#define POPULATIONTOOLS_THARDYWEINBERGTEST_H_

#include <array>
#include <map>
#include <string>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "genometools/VCF/TPopulation.h"
#include "genometools/VCF/TVcfFile.h"

namespace coretools { class TOutputFile; }

namespace PopulationTools{


//------------------------------------------------
//THWHetProbDB
//------------------------------------------------
class THWHetProbVector{
private:
	size_t _N;
	std::map<size_t, THWHetProb> _probs;

public:
	THWHetProbVector(size_t N);
	const THWHetProb& getProbs(size_t n_A);
};

class THWProbDB{
private:
	std::map<size_t, THWHetProbVector> _probs;

public:
	const THWHetProb& getProbs(size_t N, size_t n_A);
};

//------------------------------------------------
//THWGenotypes
//------------------------------------------------
class THWGenotypes{
private:
	std::array<size_t,3> _genoCounts;

public:
	THWGenotypes(){ clear(); };

	void clear();
	void add(const genometools::BiallelicGenotype & genotype);
	size_t N() const;
	size_t MAF() const;
	size_t n_A() const;
	size_t n_AB() const;
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
	THWPopulations(size_t numPops);
	void clear();
	void resize(size_t numPops);
	void add(size_t pop, const genometools::BiallelicGenotype & genotyp);
	void addToHeader(std::vector<std::string> & header);
	void runTest(coretools::TOutputFile & out);
};

//------------------------------------------------
//THardyWeinbergTest
//------------------------------------------------
class TUTest{
private:
	std::string _outname;

	//vcf-file
	std::string _vcfFilename;
    genometools::TVcfFileSingleLine _vcfFile;
	bool _limitLines;
	size_t _maxNumLines;

	//samples
    genometools::TPopulationSamples _samples;

	//genotype data
	THWPopulations _populations;

	void _openVCF();
	void _closeVCF();

public:
	TUTest();
	void run();
};

} // namespace PopulationTools

#endif /* POPULATIONTOOLS_THARDYWEINBERGTEST_H_ */
