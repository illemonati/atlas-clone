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
#include "TRandomGenerator.h"
#include <algorithm>

namespace GenotypeLikelihoods{

#define _MINLIKELIHOODVALUE 1.0E-200

//--------------------------------------------------------------------
// TBaseData
//--------------------------------------------------------------------
class TBaseData{
private:
	double _data[4];

public:
	TBaseData();
	TBaseData(const double val);
	TBaseData(const Base trueBase, const double error);
	void operator=(const TBaseData & other);
	void operator+=(const TBaseData & other);
	void operator*=(const TBaseData & other);

	double* data(){ return _data; };
	double& operator[](const Base base){ return _data[base];};
	double& operator[](const uint8_t base){ return _data[base];};
	double at(const Base base) const{ return _data[base]; };
	double at(const uint8_t base) const{ return _data[base]; };

	void set(const double val);
	void set(const Base trueBase, const double error);
	void reset();
	void add(const TBaseData & other);
	void add(const Base base, const double value);
	double sum() const;
	double weightedSum(const TBaseData & weights) const;
	void normalize();
};

std::ostream& operator<<(std::ostream& os, const TBaseData & baseData);


//--------------------------------------------------------------------
// TBaseCounts
//TODO:: merge with base frequencies?
//--------------------------------------------------------------------
class TBaseCounts{
private:
	uint32_t _counts[5];
public:
	TBaseCounts();

	void reset();
	uint32_t& operator[](const Base base){ return _counts[base]; };
	uint32_t& operator[](const uint8_t base){ return _counts[base]; };
	uint32_t at(const Base base) const{ return _counts[base]; };
	uint32_t at(const uint8_t base) const{ return _counts[base]; };
	uint32_t* data(){ return _counts; };

	void add(const Base base){ ++_counts[base]; };
	uint32_t size() const;
	uint8_t numAlleles() const;
	void fillFrequencies(TBaseData & freq);
	void fillCumulativeFrequencies(TBaseData & freq);
	void downsample(const uint32_t & max, TRandomGenerator & RandomGenerator);
};

//--------------------------------------------------------------------
// TGenotypeData
// base class for likelihoods, prior and posterior
//--------------------------------------------------------------------
class TGenotypeData{
protected:
	double _data[10];

	void _copyFrom(const TGenotypeData & other);

public:
	TGenotypeData(){};
	virtual ~TGenotypeData(){};

	void operator=(const TGenotypeData & other);
	double& operator[](const Genotype genotype){ return _data[genotype]; };
	double& operator[](const uint8_t genotype){ return _data[genotype]; };
	double at(const Genotype genotype) const { return _data[genotype]; };
	double at(const uint8_t genotype) const { return _data[genotype]; };
	double* pointerToData(){ return _data; };

	void set(const double val);
	virtual void reset();
	void add(const TGenotypeData & other);
	double sum();
	virtual double weightedSum(const TGenotypeData & weights);
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

	virtual void fill(const std::vector<TBaseData> & bases);
	virtual void fill(const std::vector<TBaseData> & bases, const size_t size);

	void addNames(std::vector<std::string> & vec, const TGenotypeMap & genoMap) const;
};

//--------------------------------------------------------------------
// TGenotypeLikelihoodsHaploid
//--------------------------------------------------------------------
class TGenotypeLikelihoodsHaploid:public TGenotypeLikelihoods{
public:
	TGenotypeLikelihoodsHaploid();

	void reset();
	void fill(const std::vector<TBaseData> & bases, const size_t size);
	double weightedSum(const TGenotypeData & weights);
};

//--------------------------------------------------------------------
// TGenotypeProbabilities
//--------------------------------------------------------------------
class TGenotypeProbabilities:public TGenotypeData{
public:
	TGenotypeProbabilities();
	void reset();
	void fill(const TGenotypeData & likelihoods, const TGenotypeData & prior);
	double probHomozygous();
	double probHeterozygous();

	void addNames(std::vector<std::string> & vec, const TGenotypeMap & genoMap) const;
};

}; //end namespace


#endif /* TGENOTYPEDATA_H_ */
