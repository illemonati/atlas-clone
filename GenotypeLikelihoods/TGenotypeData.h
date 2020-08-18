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

//TODO: use template base class for all

//--------------------------------------------------------------------
// TData_base
//--------------------------------------------------------------------

/*
template <typename T, > class TData_base{
protected:
	T _data;

public:

	template <typename T> double& operator[](const T & index){ return _data[index];};
	template <typename T> const double& operator[](const T & index){ return _data[index];};

	void operator=(const TData_base & other);
	void operator+=(const TData_base & other);
	void operator*=(const TData_base & other);

	size_t size() const { return _data.size(); };
	double min() const;
	double max() const;
};
*/

//--------------------------------------------------------------------
// TBaseData
//--------------------------------------------------------------------
class TBaseData{
private:
	std::array<double, 4> _data;

public:
	TBaseData();
	TBaseData(const double val);
	TBaseData(const Base trueBase, const double error);
	void operator=(const TBaseData & other);
	void operator+=(const TBaseData & other);
	void operator*=(const TBaseData & other);

	size_t size() const { return 4; };
	double* data(){ return _data.data(); };
	double& operator[](const Base base){ return _data[base];};
	double& operator[](const uint8_t base){ return _data[base];};
	const double& operator[](const Base base) const { return _data[base];};
	const double& operator[](const uint8_t base) const { return _data[base];};
	double at(const Base base) const{ return _data[base]; };
	double at(const uint8_t base) const{ return _data[base]; };

	void set(const double val);
	void set(const Base trueBase, const double error);
	void reset();
	void add(const TBaseData & other);
	void add(const Base base, const double value);
	double sum() const;
	double min() const;
	double max() const;
	double weightedSum(const TBaseData & weights) const;
	void normalize();

	std::array<double, 4>::iterator begin(){ return _data.begin(); };
	const std::array<double, 4>::const_iterator cbegin() const { return _data.cbegin(); };
	std::array<double, 4>::iterator end(){ return _data.end(); };
	const std::array<double, 4>::const_iterator cend() const { return _data.cend(); };
};

std::ostream& operator<<(std::ostream& os, const TBaseData & baseData);


//--------------------------------------------------------------------
// TBaseCounts
//TODO:: merge with base frequencies?
//--------------------------------------------------------------------
class TBaseCounts{
private:
	std::array<uint32_t, 5> _counts;
public:
	TBaseCounts();

	void reset();
	uint32_t& operator[](const Base base){ return _counts[base]; };
	uint32_t& operator[](const uint8_t base){ return _counts[base]; };
	const uint32_t& operator[](const Base base) const { return _counts[base]; };
	const uint32_t& operator[](const uint8_t base) const { return _counts[base]; };
	uint32_t at(const Base base) const{ return _counts[base]; };
	uint32_t at(const uint8_t base) const{ return _counts[base]; };
	uint32_t* data(){ return _counts.data(); };

	void add(const Base base){ ++_counts[base]; };
	uint32_t size() const;
	uint8_t numAlleles() const;
	void fillFrequencies(TBaseData & freq);
	void fillCumulativeFrequencies(TBaseData & freq);
	void downsample(const uint32_t & max, TRandomGenerator & RandomGenerator);

	std::array<uint32_t, 5>::iterator begin(){ return _counts.begin(); };
	const std::array<uint32_t, 5>::const_iterator cbegin() const { return _counts.cbegin(); };
	std::array<uint32_t, 5>::iterator end(){ return _counts.end(); };
	const std::array<uint32_t, 5>::const_iterator cend() const { return _counts.cend(); };
};

//--------------------------------------------------------------------
// TGenotypeData
// base class for likelihoods, prior and posterior
//--------------------------------------------------------------------
class TGenotypeData{
protected:
	std::array<double, 10> _data;

	void _copyFrom(const TGenotypeData & other);

public:
	TGenotypeData(){};
	virtual ~TGenotypeData(){};

	size_t size() const { return 10; };
	void operator=(const TGenotypeData & other);
	double& operator[](const Genotype genotype){ return _data[genotype]; };
	double& operator[](const uint8_t genotype){ return _data[genotype]; };
	const double& operator[](const Genotype genotype) const { return _data[genotype]; };
	const double& operator[](const uint8_t genotype) const { return _data[genotype]; };
	double at(const Genotype genotype) const { return _data[genotype]; };
	double at(const uint8_t genotype) const { return _data[genotype]; };
	double* pointerToData(){ return _data.data(); };

	void set(const double val);
	virtual void reset();
	void add(const TGenotypeData & other);
	double sum() const;
	double min() const;
	double max() const;
	virtual double weightedSum(const TGenotypeData & weights);
	void normalize();
	std::array<double, 10>::iterator begin(){ return _data.begin(); };
	std::array<double, 10>::const_iterator cbegin() const{ return _data.cbegin(); };
	std::array<double, 10>::iterator end(){ return _data.end(); };
	std::array<double, 10>::const_iterator cend() const{ return _data.cend(); };

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
