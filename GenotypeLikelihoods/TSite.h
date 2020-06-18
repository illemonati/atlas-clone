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

namespace GenotypeLikelihoods{

#define maxQualToPrint 1000
#define maxQualToPrintNaturalScale 1E-100


//---------------------------------------------------------------
//TSite
//---------------------------------------------------------------
class TSite{
protected:
	void normalizeGenotypeLikelihoods(double* emissionProbabilitiesPhredScaled, uint8_t* normalizedGL, uint32_t & maxLL, const int nGenotypes);

	std::vector<BAM::TBase*> bases;
	bool hasData;
	Base referenceBase; //optional

public:
	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;


	TSite(){
		hasData = false;
		referenceBase = N;
	};

	//TSite(TSite* other):TSite(){stealFromOther(other);};
	virtual ~TSite(){ clear(); };
	void clear();
	//void stealFromOther(TSite* other);

	const BAM::TBase& at(size_t i) const{ return *bases[i]; };
	BAM::TBase& operator[](size_t i){ return *bases[i]; };

	void add(const BAM::TBase * base);
	void setRefBase(const Base ref){ referenceBase = ref; };
	Base getRefBase() const {return referenceBase;};
	void addToBaseFrequencies(TBaseData & frequencies) const;
	std::string getBases(const TGenotypeMap & genoMap) const;
	std::string getQualities(const TQualityMap & qualMap) const;
	bool empty() const{ return bases.empty(); };
	uint32_t depth() const;
	uint32_t refDepth() const;

	void countAlleles(TBaseCounts & alleleCounts) const;
	void countMates(int* mateCounts) const;
	void countFwdRev(int* frCounts) const;

	//loop
	std::vector<BAM::TBase*>::iterator begin(){ return bases.begin(); };
	std::vector<BAM::TBase*>::iterator end(){ return bases.end(); };
	std::vector<BAM::TBase*>::const_iterator cbegin() const{ return bases.cbegin(); };
	std::vector<BAM::TBase*>::const_iterator cend() const{ return bases.cend(); };
};

}; //end namespace

#endif /* TSITE_H_ */
