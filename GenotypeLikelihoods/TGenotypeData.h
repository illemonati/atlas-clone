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
class TBaseData{
private:
	double data[4];

public:
	TBaseData();
	TBaseData(const double val);
	TBaseData(const Base trueBase, const double error);
	void operator=(const TBaseData & other);
	void operator+=(const TBaseData & other);
	void operator*=(const TBaseData & other);
	double& operator[](const Base base){ return data[base];};
	double at(const Base base) const { return data[base]; };

	void set(const double val);
	void set(const Base trueBase, const double error);
	void reset();
	void add(const TBaseData & other);
	double sum() const;
	double weightedSum(const TBaseData & weights) const;
	uint8_t numAlleles() const;
	void normalize();
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

	virtual void addNames(std::vector<std::string> & vec, const TGenotypeMap & genoMap) const;
	void write(TOutputFile & out) const;
};


//--------------------------------------------------------------------
// TGenotypeLikelihoods
//--------------------------------------------------------------------
class TGenotypeLikelihoods:public TGenotypeData{
public:
	TGenotypeLikelihoods();

	void operator=(const TGenotypeLikelihoods & other);

	void fill(const std::vector<TBaseData> & bases);
	void fill(const std::vector<TBaseData> & bases, const size_t size);

	void addNames(std::vector<std::string> & vec, const TGenotypeMap & genoMap) const;
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

	void addNames(std::vector<std::string> & vec, const TGenotypeMap & genoMap) const;
};

}; //end namespace


#endif /* TGENOTYPEDATA_H_ */
