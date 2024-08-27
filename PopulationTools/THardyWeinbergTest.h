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

#include "genometools/VCF/TPopulation.h"
#include "genometools/VCF/TVcfFile.h"

namespace coretools { class TOutputFile; }

namespace PopulationTools{

//------------------------------------------------
//THWHetProb
// Stores the probabilities P(N_AB = n|N, n_A, n_B),
// i.e the probabilities of observing n hets among N individuals if the allele frequencies are n_A + n_B = N
//------------------------------------------------
class THWHetProb{
private:
	std::vector<double> _probs;
	size_t _maxNumHetPlusOne;
	bool _onlyOdd, _onlyEven;

public:
	THWHetProb();
	THWHetProb(size_t N, size_t n_A);

	void clear();
	THWHetProb& operator=(const THWHetProb & other);
	void extend(const THWHetProb & other);
	bool onlyOdd() const { return _onlyOdd; };
	bool onlyEven() const { return _onlyEven; };
	size_t maxNumHet() const { return _maxNumHetPlusOne - 1; };
	double operator[](size_t i) const { return _probs[i]; };
	double sum(size_t upTo);
	void print();
};

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
class THardyWeinbergTest{
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
	THardyWeinbergTest();
	void run();
};

} // namespace PopulationTools

#endif /* POPULATIONTOOLS_THARDYWEINBERGTEST_H_ */
