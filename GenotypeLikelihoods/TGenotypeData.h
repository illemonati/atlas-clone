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
#include "probability.h"

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
/*
	virtual void operator+=(const TData_base & other){
		for(size_t i = 0; i<_data.size(); ++i){
			_data[i] += other[i];
		}
	};
	*/

	virtual void operator*=(const TData_base & other){
		for(size_t i = 0; i<_data.size(); ++i){
			_data[i] *= other[i];
		}
	};

	/*
	virtual void add(const TData_base & other){
		*this += other;
	};
	*/

	T min() const{
		return *std::min_element(_data.begin(), _data.end());
	};

	T max() const{
		return *std::max_element(_data.begin(), _data.end());
	};

	auto begin(){ return _data.begin(); };
	auto cbegin() const { return _data.cbegin(); };
	auto end(){ return _data.end(); };
	auto cend() const { return _data.cend(); };
};


//--------------------------------------------------------------------
// TBaseDataq
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
class TGenotypeData:public TData_base<Probability, 10>{
protected:
	//void _copyFrom(const TGenotypeData & other);

	constexpr double _sum() const{
		return _data[BAM::AA].get() + _data[BAM::AC].get() + _data[BAM::AG].get() + _data[BAM::AT].get() + _data[BAM::CC].get() + _data[BAM::CG].get() + _data[BAM::CT].get() + _data[BAM::GG].get() + _data[BAM::GT].get() + _data[BAM::TT].get();
	};
public:
	TGenotypeData(){};
	virtual ~TGenotypeData(){};

	void operator=(const TGenotypeData & other);
	void set(const double val);
	virtual void reset();
	void add(const TGenotypeData & other);
	virtual double weightedSum(const TGenotypeData & weights);
	void normalize(const Probability & theSum);
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
