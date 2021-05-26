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

	void set(const Type & val){
		std::fill(begin(), end(), val);
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
		std::transform(_data.begin(), _data.end(), other._data.begin(), _data.begin(), std::multiplies<Type>() );
	};

	Type sum(){
		return std::accumulate(_data , _data+Size, Type{});
	};

	void normalize(){
		Type tot = accumulate(_data , _data+Size, Type{});
		std::for_each(_data.begin(), _data.end(), [tot](Type &x){ x /= tot; });
	};

	void normalize(const Type & tot){
		std::for_each(_data.begin(), _data.end(), [tot](Type &x){ x /= tot; });
	};

	template <typename U>
	double weightedSum(const U & weights) const{
		auto wsum = [](const double & s, const Type & v) { return s + (double) v; };
		return std::accumulate(_data.begin(), _data.end(), double{}, wsum);
	};

	virtual void addNames(std::vector<std::string> & vec) const{
		for(uint8_t i=0; i<Size; ++i){
			vec.push_back( (std::string) IndexType(static_cast<EnumType>(i)));
		}
	};

	void write(TOutputFile & out) const{
		std::for_each(_data.begin(), _data.end(), [out](const Type & v){ out << v; });
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
	using TData_base<T, BAM::Base, BAM::BaseEnum, 4>::set;

	//keep constructors protected so the base class can not be used!
	TBaseData_base() = default;
	TBaseData_base(const T& val){ set(val); };
};

//--------------------------------------------------------------------
// TBaseData
//--------------------------------------------------------------------
class TBaseData:public TBaseData_base<double>{
public:
	TBaseData(){ reset(); };
	TBaseData(const double & val) : TBaseData_base(val) {};
};

//--------------------------------------------------------------------
// TBaseProbabilities
//--------------------------------------------------------------------
class TBaseProbabilities:public TBaseData_base<Probability>{
public:
	TBaseProbabilities(){ reset(); };
	TBaseProbabilities(const Probability & val) : TBaseData_base(val) {};
};

//--------------------------------------------------------------------
// TBaseLikelihoods
//--------------------------------------------------------------------
class TBaseLikelihoods:public TBaseData_base<Probability>{
public:
	TBaseLikelihoods(){ reset(); };
	TBaseLikelihoods(const Probability & val) : TBaseData_base(val) {};

	TBaseLikelihoods(const BAM::Base & trueBase, const Probability & error);

	void reset() override { set(Probability(1.0)); };

	void setFromError(const BAM::Base trueBase, const Probability & error);
};

//--------------------------------------------------------------------
// TBaseCounts
//TODO:: merge with base frequencies?
//--------------------------------------------------------------------
class TBaseCounts:public TData_base<uint32_t, BAM::Base, BAM::BaseEnum, 5>{
public:
	TBaseCounts(){ reset(); };

	void add(const BAM::Base base){ ++_data[static_cast<uint8_t>( base.get() )]; };

	uint8_t numAlleles() const;
	void fillFrequencies(TBaseData & freq);
	void fillCumulativeFrequencies(TBaseProbabilities & freq);
	void downsample(const uint32_t & max, TRandomGenerator & RandomGenerator);
};

//--------------------------------------------------------------------
// TGenotypeData_base
// base class for TGenotypeData, likelihoods, prior and posterior
//--------------------------------------------------------------------
template <typename T>
class TGenotypeData_base : public TData_base<T, BAM::Genotype, BAM::GenotypeEnum, 10>{
protected:
	using TData_base<T, BAM::Genotype, BAM::GenotypeEnum, 10>::_data;

	//keep constructors protected so the base class can not be used!
	TGenotypeData_base(){};
	virtual ~TGenotypeData_base(){};
};

//--------------------------------------------------------------------
// TGenotypeData
// base class for TGenotypeLikelihoods, TGenotypeLikelihoodsHaploid and TGenotypeProbabilities
//--------------------------------------------------------------------
class TGenotypeData : public TGenotypeData_base<double>{
public:
	TGenotypeData(){ reset(); };
};

//--------------------------------------------------------------------
// TGenotypeLikelihoods
//--------------------------------------------------------------------
class TGenotypeLikelihoods : public TGenotypeData_base<Probability>{
public:
	TGenotypeLikelihoods();

	virtual void fill(const std::vector<TBaseData> & bases);
	virtual void fill(const std::vector<TBaseData> & bases, const size_t size);

	void reset() override { set(Probability(1.0)); };
	void addNames(std::vector<std::string> & vec) const override;
};

//--------------------------------------------------------------------
// TGenotypeProbabilities
//--------------------------------------------------------------------
class TGenotypeProbabilities : public TGenotypeData_base<Probability>{
public:
	TGenotypeProbabilities() { reset(); };

	void reset() override { set(Probability(0.1)); };

	void fill(const TGenotypeLikelihoods & likelihoods, const TGenotypeData & prior);

	double probHomozygous();
	double probHeterozygous();
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
