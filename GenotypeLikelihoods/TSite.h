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

public:
	std::vector<TBase*> bases;
	bool hasData;
	Base referenceBase; //optional

	GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoods;


	TSite(){
		hasData = false;
		referenceBase = 'N';
	};

	//TSite(TSite* other):TSite(){stealFromOther(other);};
	virtual ~TSite(){ clear(); };
	void clear();
	//void stealFromOther(TSite* other);

	void add(const TBase * base);
	void setRefBase(const Base ref){ referenceBase = ref; };
	Base getRefBase() const {return referenceBase;};
	void addToBaseFrequencies(TBaseData & frequencies) const;
	std::string getBases(const TGenotypeMap & genoMap) const;
	uint32_t depth() const;
	uint32_t refDepth() const;

	void countAlleles(TBaseCounts & alleleCounts) const;
	void countMates(int* mateCounts) const;
	void countFwdRev(int* frCounts) const;
};

}; //end namespace

#endif /* TSITE_H_ */
