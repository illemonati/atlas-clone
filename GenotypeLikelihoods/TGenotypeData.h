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
template <typename Type, typename IndexType, typename EnumType, size_t Size>
class TData_base{
protected:
	std::array<Type, Size> _data;

	//keep protected so the base class can not be used
	TData_base() = default;
	TData_base(const Type & val){ set(val); };
	 virtual ~TData_base(){};

public:
	size_t size() const { return _data.size(); };
	Type* pointerToData(){ return _data.data(); };
	size_t sizeOf() const { return sizeof(Type) * Size; };

	template <typename T>
	void set(const T & val){
		std::fill(begin(), end(), Type(val));
	};

	virtual void reset(){ set(Type{}); };

	void operator=(const TData_base & other){
		_data = other._data;
	};

	Type min() const{
		return *std::min_element(_data.begin(), _data.end());
	};

	Type max() const{
		return *std::max_element(_data.begin(), _data.end());
	};

	IndexType pickIndexAtMax(coretools::TRandomGenerator & RandomGenerator) const {
		//find maximum
		Type m = max();

		//get vec of all indexes at maximum
		std::vector<IndexType> vec;
		for(IndexType i=IndexType::min(); i<IndexType::max(); ++i){
			if(_data[i.get()] == m){
				vec.push_back(i);
			}
		}

		//return random index among those at max
		return vec[RandomGenerator.sample(vec.size())];
	};

	IndexType pickIndexAtMax(const IndexType & ExcludeIndex, coretools::TRandomGenerator & RandomGenerator) const {
		//find max
		IndexType i = IndexType::min();
		if(i == ExcludeIndex){
			++i;
		}
		Type m = _data[i.get()];
		++i;

		for(; i<IndexType::max(); ++i){
			if(i!= ExcludeIndex && _data[i.get()] > m){
				m = _data[i.get()];
			}
		}

		//get vec of all index at maximum
		std::vector<IndexType> vec;
		for(IndexType i=IndexType::min(); i<IndexType::max(); ++i){
			if(i != ExcludeIndex && _data[i.get()] == m){
				vec.push_back(i);
			}
		}

		//return random index among those at max
		return vec[RandomGenerator.sample(vec.size())];
	};

	Type& operator[](const IndexType & index){ return _data[index.get()]; };
	const Type& operator[](const IndexType & index) const { return _data[index.get()]; };
	Type& operator[](const EnumType & index){ return _data[index]; };
	const Type& operator[](const EnumType & index) const { return _data[index]; };

	auto begin(){ return _data.begin(); };
	auto cbegin() const { return _data.cbegin(); };
	auto end(){ return _data.end(); };
	auto cend() const { return _data.cend(); };

	void operator+=(const TData_base & other){
		std::transform(_data.begin(), _data.end(), other._data.begin(), _data.begin(), std::plus<Type>() );
	};

	void addConstant(const IndexType & index, const Type & value){
		_data[index] += value;
	};

	void operator*=(const TData_base & other){
		for(uint8_t i=1; i<Size; ++i){
			_data[i] *= other._data[i];
		}
	};

	Type sum() const{
		return std::accumulate(_data.begin(), _data.end(), Type{});
	};

	void normalize(){
		Type tot = sum();
		std::for_each(_data.begin(), _data.end(), [tot](Type &x){ x /= tot; });
	};

	void normalize(const Type & tot){
		std::for_each(_data.begin(), _data.end(), [tot](Type &x){ x /= tot; });
	};

	template <typename U>
	double weightedSum(const U & weights) const {
		auto wsum = [](const double & s, const Type & v) { return s + (double) v; };
		return std::accumulate(_data.begin(), _data.end(), double{}, wsum);
	};

	virtual void addNames(std::vector<std::string> & vec) const{
		for(uint8_t i=0; i<Size; ++i){
			vec.push_back( (std::string) IndexType(static_cast<EnumType>(i)));
		}
	};

	void write(coretools::TOutputFile & out) const {
		for(uint8_t i=0; i<Size; ++i){
			out << coretools::str::toString(_data[i]);
		}
	};

	friend std::ostream& operator<<(std::ostream& os, const TData_base & data){
		os << (std::string) IndexType(static_cast<EnumType>(0)) + ": " << data._data[0];
		for(uint8_t i=1; i<Size; ++i){
			os << ", " + (std::string) IndexType(static_cast<EnumType>(i)) + ": " << data._data[i];
		}
		return os;
	};
};

//--------------------------------------------------------------------
// TBaseData_base
//--------------------------------------------------------------------
template <typename T>
class TBaseData_base : public TData_base<T, BAM::Base, BAM::BaseEnum, 4>{
protected:
	using TData_base<T, BAM::Base, BAM::BaseEnum, 4>::_data;

	//keep constructors protected so the base class can not be used!
	TBaseData_base() = default;
	TBaseData_base(const T& val){ set(val); };
public:
	using TData_base<T, BAM::Base, BAM::BaseEnum, 4>::set;
};

//--------------------------------------------------------------------
// TBaseProbabilities
//--------------------------------------------------------------------
class TBaseProbabilities : public TBaseData_base<coretools::Probability>{
public:
	TBaseProbabilities(){ reset(); };
	TBaseProbabilities(const coretools::Probability & val) : TBaseData_base(val) {};
};

//--------------------------------------------------------------------
// TBaseLikelihoods
//--------------------------------------------------------------------
class TBaseLikelihoods : public TBaseProbabilities{
public:
	TBaseLikelihoods(){ reset(); };
	TBaseLikelihoods(const coretools::Probability & val) : TBaseProbabilities(val) {};

	TBaseLikelihoods(const BAM::Base & trueBase, const coretools::Probability & error);

	void reset() override { set(1.0); };

	void setFromError(const BAM::Base & trueBase, const coretools::Probability & error);
};

//--------------------------------------------------------------------
// TBaseCounts
//--------------------------------------------------------------------
class TBaseCounts:public TData_base<uint32_t, BAM::Base, BAM::BaseEnum, 5>{
private:
	void _fillCumulativeFrequencies(std::array<double, 4> & freq);

public:
	TBaseCounts(){ reset(); };

	void add(const BAM::Base base){ ++_data[static_cast<uint8_t>( base.get() )]; };

	uint8_t numAlleles() const;
	void fillFrequencies(TBaseProbabilities & freq);
	void downsample(const uint32_t & max, coretools::TRandomGenerator & RandomGenerator);
};

//--------------------------------------------------------------------
// TBaseData
//--------------------------------------------------------------------
class TBaseData:public TBaseData_base<double>{
public:
	TBaseData(){ reset(); };
	TBaseData(const double & val) : TBaseData_base(val) {};

	[[nodiscard]] TBaseProbabilities asFrequencies();

	void operator+=(const TBaseProbabilities & probs);
};

//--------------------------------------------------------------------
// TGenotypeData_base
// base class for TGenotypeData, likelihoods, prior and posterior
//--------------------------------------------------------------------
template <typename T>
class TGenotypeData_base : public TData_base<T, BAM::Genotype, BAM::GenotypeEnum, 10>{
protected:
	using TData_base<T, BAM::Genotype, BAM::GenotypeEnum, 10>::_data;

public:
	TGenotypeData_base(){};
	virtual ~TGenotypeData_base(){};
};

//-------------------------------------
// TGenotypeProbability_base
// base class for TGenotypeLikelihoods, TGenotypeLikelihoodsHaploid and TGenotypeProbabilities
//-------------------------------------
class TGenotypeProbability_base : public TGenotypeData_base<coretools::Probability>{
protected:
	TGenotypeProbability_base() {};
};

//--------------------------------------------------------------------
// TGenotypeLikelihoods
//--------------------------------------------------------------------
class TGenotypeLikelihoods : public TGenotypeProbability_base{
public:
	TGenotypeLikelihoods();

	virtual void fill(const std::vector<TBaseLikelihoods> & bases);
	virtual void fill(const std::vector<TBaseLikelihoods> & bases, const size_t size);

	void reset() override { set(coretools::Probability(1.0)); };
	void addNames(std::vector<std::string> & vec) const override;
};

//--------------------------------------------------------------------
// TGenotypeProbabilities
//--------------------------------------------------------------------
class TGenotypeProbabilities : public TGenotypeProbability_base{
public:
	TGenotypeProbabilities() { reset(); };

	void reset() override { set(coretools::Probability(0.1)); };

	void fillBayesian(const TGenotypeLikelihoods & likelihoods, const TGenotypeProbabilities & prior);

	coretools::Probability probHomozygous();
	coretools::Probability probHeterozygous();
};

//--------------------------------------------------------------------
// TGenotypeData
//--------------------------------------------------------------------
class TGenotypeData : public TGenotypeData_base<double>{
public:
	TGenotypeData(){ reset(); };

	void operator+=(const TGenotypeProbability_base & other);
};

//--------------------------------------------------------------------
// TGenotypeLikelihoodsHaploid
//--------------------------------------------------------------------
class TGenotypeLikelihoodsHaploid:public TGenotypeLikelihoods{
public:
	TGenotypeLikelihoodsHaploid(){ reset(); };

	void reset() override;
	void fill(const std::vector<TBaseData> & bases, const size_t size);
};


}; //end namespace


#endif /* TGENOTYPEDATA_H_ */
