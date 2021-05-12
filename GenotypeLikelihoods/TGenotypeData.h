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
// TData_base
//--------------------------------------------------------------------
template <typename T, size_t S> class TData_base{
protected:
	std::array<T, S> _data;

public:
	TData_base() = default;
	TData_base(const T & val){ set(val); };
	virtual ~TData_base(){};
	template <typename U> T& operator[](const U & index){ return _data[index]; };
	template <typename U> const T& operator[](const U & index) const { return _data[index]; };
	size_t size() const { return _data.size(); };
	T* pointerToData(){ return _data.data(); };
	size_t sizeOf() const { return sizeof(T) * S; };

	virtual void set(const T & val){
		std::fill(begin(), end(), val);
	};

	virtual void reset(){
		set({});
	};

	virtual void operator+=(const TData_base & other){
		for(size_t i = 0; i<_data.size(); ++i){
			_data[i] += other[i];
		}
	};

	virtual void operator*=(const TData_base & other){
		for(size_t i = 0; i<_data.size(); ++i){
			_data[i] *= other[i];
		}
	};

	virtual void add(const TData_base & other){
		*this += other;
	};

	T min() const{
		return *std::min_element(_data.begin(), _data.end());
	};

	T max() const{
		return *std::max_element(_data.begin(), _data.end());
	};

	virtual T sum() const{
		T s{};
		for(auto& i : _data){
			s += i;
		}
		return s;
	};

	auto begin(){ return _data.begin(); };
	auto cbegin() const { return _data.cbegin(); };
	auto end(){ return _data.end(); };
	auto cend() const { return _data.cend(); };
};


//--------------------------------------------------------------------
// TBaseData
//--------------------------------------------------------------------
class TBaseData:public TData_base<double, 4>{
public:
	TBaseData();
	TBaseData(const double val);
	TBaseData(const BAM::Base trueBase, const double error);

	void operator=(const TBaseData & other);
	void operator+=(const TBaseData & other);
	void operator*=(const TBaseData & other);

	void set(const double val);
	void set(const BAM::Base trueBase, const double error);
	void reset();
	void add(const TBaseData & other);
	void add(const BAM::Base base, const double value);
	double sum() const;
	double weightedSum(const TBaseData & weights) const;
	void normalize();
};

std::ostream& operator<<(std::ostream& os, const TBaseData & baseData);


//--------------------------------------------------------------------
// TBaseCounts
//TODO:: merge with base frequencies?
//--------------------------------------------------------------------
class TBaseCounts:public TData_base<uint32_t, 5>{
public:
	TBaseCounts();

	void reset();

	void add(const BAM::Base base){ ++_data[static_cast<uint8_t>( base.get() )]; };
	inline uint32_t sum() const{
		return _data[BAM::A] + _data[BAM::C] + _data[BAM::G] + _data[BAM::T] + _data[BAM::N];
	};
	uint8_t numAlleles() const;
	void fillFrequencies(TBaseData & freq);
	void fillCumulativeFrequencies(TBaseData & freq);
	void downsample(const uint32_t & max, TRandomGenerator & RandomGenerator);
};

//--------------------------------------------------------------------
// TGenotypeData
// base class for likelihoods, prior and posterior
//--------------------------------------------------------------------
class TGenotypeData:public TData_base<double, 10>{
protected:
	//void _copyFrom(const TGenotypeData & other);

public:
	TGenotypeData(){};
	virtual ~TGenotypeData(){};

	void operator=(const TGenotypeData & other);
	void set(const double val);
	virtual void reset();
	void add(const TGenotypeData & other);
	inline double sum() const{
		return _data[BAM::AA] + _data[BAM::AC] + _data[BAM::AG] + _data[BAM::AT] + _data[BAM::CC] + _data[BAM::CG] + _data[BAM::CT] + _data[BAM::GG] + _data[BAM::GT] + _data[BAM::TT];
	};
	virtual double weightedSum(const TGenotypeData & weights);
	void normalize();

	virtual void addNames(std::vector<std::string> & vec) const;
	void write(TOutputFile & out) const;
};


//--------------------------------------------------------------------
// TGenotypeLikelihoods
//--------------------------------------------------------------------
class TGenotypeLikelihoods:public TGenotypeData{
public:
	TGenotypeLikelihoods();

	virtual void fill(const std::vector<TBaseData> & bases);
	virtual void fill(const std::vector<TBaseData> & bases, const size_t size);

	void addNames(std::vector<std::string> & vec) const;
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

	void addNames(std::vector<std::string> & vec) const;
};

}; //end namespace


#endif /* TGENOTYPEDATA_H_ */
