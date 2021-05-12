/*
 * TBase.h
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#ifndef TSITE_H_
#define TSITE_H_

#include <vector>
#include "../BAM/TSequencedBase.h"
#include "Types.h"
#include "TGenotypeData.h"
#include "TSubsamplePicker.h"

namespace GenotypeLikelihoods{

#define maxQualToPrint 1000
#define maxQualToPrintNaturalScale 1E-100

//TODO: write as templated classes
//----------------------------------------------------------------------------------------------------------------------------------
// TSite
// Class that stores bases.
//----------------------------------------------------------------------------------------------------------------------------------
class TSite {
protected:
    BAM::Base _referenceBase;
    BAM::Genotype _genotype;

	std::vector<BAM::TSequencedBase> _bases;

	void normalizeGenotypeLikelihoods(double* emissionProbabilitiesPhredScaled, uint8_t* normalizedGL, uint32_t & maxLL, const int nGenotypes);

public:
	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;

    TSite(){
        _referenceBase = BAM::N;
        _genotype = BAM::NN;
    };

	void clear();

	// access
	BAM::TSequencedBase& operator[](size_t i){ return _bases[i]; };
	const BAM::TSequencedBase& operator[](size_t i) const { return _bases[i]; };

	// reference base
    void setRefBase(const BAM::Base & ref){ _referenceBase = ref; };
    BAM::Base refBase() const {return _referenceBase;};

    // genotype
    void setGenotype(const BAM::Genotype & genotype){ _genotype = genotype; };
    BAM::Genotype genotype() const{ return _genotype; };

    // add
    void add(const BAM::TSequencedBase & base);
    void addToBaseFrequencies(TBaseData & frequencies) const;
	void downsample(const uint32_t & maxDepth, const TSubsamplePicker & picker);

	// getters
	bool empty() const{ return _bases.empty(); };
	uint32_t depth() const;
	uint32_t refDepth() const;
	std::string getBases() const;
	std::string getQualities() const;

	void countAlleles(TBaseCounts & alleleCounts) const;
	void countMates(int* mateCounts) const;
	void countFwdRev(int* frCounts) const;

	//loop
	std::vector<BAM::TSequencedBase>::iterator begin(){ return _bases.begin(); };
	std::vector<BAM::TSequencedBase>::iterator end(){ return _bases.end(); };
	std::vector<BAM::TSequencedBase>::const_iterator cbegin() const{ return _bases.cbegin(); };
	std::vector<BAM::TSequencedBase>::const_iterator cend() const{ return _bases.cend(); };
};

}; //end namespace

#endif /* TSITE_H_ */
