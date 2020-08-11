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
// TSite_base
// Base class not meant to be used.
//----------------------------------------------------------------------------------------------------------------------------------
class TSite_base{
protected:
	Base _referenceBase;
	Genotype _genotype;

public:
	TSite_base(){
		_referenceBase = N;
		_genotype = NN;
	};

	void setRefBase(const Base ref){ _referenceBase = ref; };
	Base refBase() const {return _referenceBase;};

	void setGenotype(const Genotype genotype){ _genotype = genotype; };
	Genotype genotype() const{ return _genotype; };
};

//----------------------------------------------------------------------------------------------------------------------------------
// TSite
// Class that used pointer to bases.
// Used in routines based in windows that store alignments with bases (no need to copy memory).
//----------------------------------------------------------------------------------------------------------------------------------
class TSite:public TSite_base{
protected:
	std::vector<BAM::TBase*> _bases;

	void normalizeGenotypeLikelihoods(double* emissionProbabilitiesPhredScaled, uint8_t* normalizedGL, uint32_t & maxLL, const int nGenotypes);

public:
	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;

	void clear();
	const BAM::TBase& at(size_t i) const{ return *_bases[i]; };
	BAM::TBase& operator[](size_t i){ return *_bases[i]; };

	void add(BAM::TBase * base);
	void addToBaseFrequencies(TBaseData & frequencies) const;
	void downsample(const uint32_t & maxDepth, const TSubsamplePicker & picker);

	bool empty() const{ return _bases.empty(); };
	uint32_t depth() const;
	uint32_t refDepth() const;
	std::string getBases(const TGenotypeMap & genoMap) const;
	std::string getQualities(const BAM::TQualityMap & qualMap) const;

	void countAlleles(TBaseCounts & alleleCounts) const;
	void countMates(int* mateCounts) const;
	void countFwdRev(int* frCounts) const;

	//loop
	std::vector<BAM::TBase*>::iterator begin(){ return _bases.begin(); };
	std::vector<BAM::TBase*>::iterator end(){ return _bases.end(); };
	std::vector<BAM::TBase*>::const_iterator cbegin() const{ return _bases.cbegin(); };
	std::vector<BAM::TBase*>::const_iterator cend() const{ return _bases.cend(); };
};

//----------------------------------------------------------------------------------------------------------------------------------
// TSiteStorage
// Class that stores bases (not pointers)
// Used in routines that need to keep sites in memory beyond scope of windows.
// Standard way of using is to first read data into TSites, then to filter sites and "copy" the data of TSite into TSiteStorage
//----------------------------------------------------------------------------------------------------------------------------------
class TSiteStorage:public TSite_base{
protected:
	std::vector<BAM::TBase> _bases;

public:
	TSiteStorage(const TSite & site);
	void clear();

	const BAM::TBase& at(size_t i) const{ return _bases[i]; };
	BAM::TBase& operator[](size_t i){ return _bases[i]; };

	void add(const BAM::TBase * base);
	void add(const BAM::TBase & base);

	void addToBaseFrequencies(TBaseData & frequencies) const;
	bool empty() const{ return _bases.empty(); };
	uint32_t depth() const;

	//loop
	std::vector<BAM::TBase>::iterator begin(){ return _bases.begin(); };
	std::vector<BAM::TBase>::iterator end(){ return _bases.end(); };
	std::vector<BAM::TBase>::const_iterator cbegin() const{ return _bases.cbegin(); };
	std::vector<BAM::TBase>::const_iterator cend() const{ return _bases.cend(); };
};


}; //end namespace

#endif /* TSITE_H_ */
