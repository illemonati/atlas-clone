/*
 * TBase.h
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#ifndef TSITE_H_
#define TSITE_H_

#include <TBase.h>
#include <vector>
#include "TQualityMap.h"
#include "TGenotypeMap.h"
#include "TGenotypeData.h"
#include "TSubsamplePicker.h"

namespace GenotypeLikelihoods{

#define maxQualToPrint 1000
#define maxQualToPrintNaturalScale 1E-100

//TODO: write as templated classes
//----------------------------------------------------------------------------------------------------------------------------------
// TSite
// Class that stores bases.
// Used in routines based in windows that store alignments with bases (no need to copy memory).
//----------------------------------------------------------------------------------------------------------------------------------
class TSite {
protected:
    Base _referenceBase;
    Genotype _genotype;

	std::vector<BAM::TBase> _bases;

	void normalizeGenotypeLikelihoods(double* emissionProbabilitiesPhredScaled, uint8_t* normalizedGL, uint32_t & maxLL, const int nGenotypes);

public:
	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;

    TSite(){
        _referenceBase = N;
        _genotype = NN;
    };

	void clear();

	// access
	const BAM::TBase& at(size_t i) const{ return _bases[i]; };
	BAM::TBase& operator[](size_t i){ return _bases[i]; };

	// reference base
    void setRefBase(const Base ref){ _referenceBase = ref; };
    Base refBase() const {return _referenceBase;};

    // genotype
    void setGenotype(const Genotype genotype){ _genotype = genotype; };
    Genotype genotype() const{ return _genotype; };

    // add
	void add(const BAM::TBase * base);
    void add(const BAM::TBase & base);
    void addToBaseFrequencies(TBaseData & frequencies) const;
	void downsample(const uint32_t & maxDepth, const TSubsamplePicker & picker);

	// getters
	bool empty() const{ return _bases.empty(); };
	uint32_t depth() const;
	uint32_t refDepth() const;
	std::string getBases(const TGenotypeMap & genoMap) const;
	std::string getQualities(const BAM::TQualityMap & qualMap) const;

	void countAlleles(TBaseCounts & alleleCounts) const;
	void countMates(int* mateCounts) const;
	void countFwdRev(int* frCounts) const;

	//loop
	std::vector<BAM::TBase>::iterator begin(){ return _bases.begin(); };
	std::vector<BAM::TBase>::iterator end(){ return _bases.end(); };
	std::vector<BAM::TBase>::const_iterator cbegin() const{ return _bases.cbegin(); };
	std::vector<BAM::TBase>::const_iterator cend() const{ return _bases.cend(); };
};

}; //end namespace

#endif /* TSITE_H_ */
