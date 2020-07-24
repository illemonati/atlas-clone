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

/*
//---------------------------------------------------------------
//TBaseFrequencies
//---------------------------------------------------------------
class TBaseFrequencies{
public:
	double freq[4];
	bool wasNormalized;

	TBaseFrequencies(){
		for(int i = 0; i < 4; ++i) freq[i] = 0.0;
		wasNormalized = false;
	};
	void add(Base B, double & weight){
		freq[B] += weight;
	};
	void addNoRef(Base B, double weight){
		freq[B] += weight;
	};
	void normalize(){
		if(!wasNormalized){
			double sum = 0.0;
			for(int i = 0; i < 4; ++i) sum += freq[i];
			sum += 4.0;
			for(int i = 0; i < 4; ++i) freq[i] = (freq[i] + 1.0) / sum;
			wasNormalized = true;
		}
	};
	void setEqualBaseFreq(){
		for(int i = 0; i < 4; ++i) freq[i] = 0.25;
	};
	void clear(){
		for(int i = 0; i < 4; ++i) freq[i] = 0.0;
		wasNormalized = false;
	};
	void print() const{
		std::cout << "freq(A) = " << freq[0] << ", freq(C) = " << freq[1] << ", freq(G) = " << freq[2] << ", freq(T) = " << freq[3] << std::endl;
	};
	double& operator[](int pos){
		return freq[pos];
	};
};
*/

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
	uint32_t* pointerToCounts(){ return _counts; };
	void add(const Base base){ ++_counts[base]; };
	uint32_t size() const;
	uint8_t numAlleles() const;
	void fillFrequencies(TBaseData & freq);
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
