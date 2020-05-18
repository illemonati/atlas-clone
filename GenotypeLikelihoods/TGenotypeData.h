/*
 * TGenotypePrior.h
 *
 *  Created on: Nov 20, 2018
 *      Author: phaentu
 */

#ifndef TGENOTYPEDATA_H_
#define TGENOTYPEDATA_H_

#include "TGenotypeMap.h"
#include "TFile.h"
#include <algorithm>

namespace GenotypeLikelihoods{


//--------------------------------------------------------------------
// TBaseLikelihoods
//--------------------------------------------------------------------
class TBaseLikelihoods{
private:
	double likelihoods[4];

public:
	TBaseLikelihoods();
	void operator=(const TBaseLikelihoods & other);
	double& operator[](const Base base){ return likelihoods[base];};
	double at(const Base base) const { return likelihoods[base]; };
	void reset();
};

//--------------------------------------------------------------------
// TGenotypeData
// base class for likelihoods, prior and posterior
//--------------------------------------------------------------------
class TGenotypeData{
protected:
	double data[10];

	void _copyFrom(const TGenotypeData & other);

public:
	TGenotypeData(){};
	virtual ~TGenotypeData(){};

	void operator=(const TGenotypeData & other);
	double& operator[](const Genotype genotype){ return data[genotype]; };
	double& operator[](const uint8_t genotype){ return data[genotype]; };
	double at(const Genotype genotype) const { return data[genotype]; };
	double at(const uint8_t genotype) const { return data[genotype]; };
	double* pointerToData(){ return data; };

	void set(const double val);
	virtual void reset();
	void add(const TGenotypeData & other);
	double sum();
	double weightedSum(const TGenotypeData & weights);
	void normalize();

	void write(TOutputFileZipped & out) const;
};


//--------------------------------------------------------------------
// TGenotypeLikelihoods
//--------------------------------------------------------------------
class TGenotypeLikelihoods:public TGenotypeData{
public:
	TGenotypeLikelihoods();

	void operator=(const TGenotypeLikelihoods & other);

	void fill(const std::vector<TBaseLikelihoods> & bases);
	void fill(const std::vector<TBaseLikelihoods> & bases, const size_t size);
};


//--------------------------------------------------------------------
// TGenotypePosteriorProbabilities
//--------------------------------------------------------------------
class TGenotypePosteriorProbabilities:public TGenotypeData{
public:
	TGenotypePosteriorProbabilities();
	void reset();
	void fill(const TGenotypeData & likelihoods, const TGenotypeData & prior);
	double probHomozygous();
	double probHeterozygous();
};

}; //end namespace


#endif /* TGENOTYPEDATA_H_ */
