/*
 * TBase.h
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#ifndef TSITE_H_
#define TSITE_H_

#include <math.h>
#include <TBase.h>
#include <vector>
#include "TParameters.h"
#include "TGenotypeMap.h"
#include "TQualityMap.h"
#include "TGenotypeData.h"
#include "gzstream.h"
#include <algorithm>
#include "TRandomGenerator.h"
#include "../TAllelicDepthCounts.h"
#include "auxiliaryTools.h"

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
	void stealFromOther(TSite* other);

	void add(const TBase * base);
	void setRefBase(const Base ref){
		referenceBase = ref;
	};
	Base getRefBase() const {return referenceBase;};
	uint32_t depth();
	uint32_t refDepth();
	void addToBaseFrequencies(TBaseFrequencies & frequencies);
	std::string getBases(const TGenotypeMap & genoMap);

	void countAlleles(GenotypeLikelihoods::TBaseData & alleleCounts) const;
	void countAllelesForImbalance(TAllelicDepthCounts & counts);
	void countMates(int* mateCounts);
	void countFwdRev(int* frCounts);
};

#endif /* TSITE_H_ */
