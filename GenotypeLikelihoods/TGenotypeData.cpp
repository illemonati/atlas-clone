/*
 * TGenotypeData.cpp
 *
 *  Created on: May 15, 2020
 *      Author: phaentu
 */

#include "TGenotypeData.h"
#include "GenotypeTypes.h"

namespace GenotypeLikelihoods{

//Note: for speed, most loops are unrolled

using coretools::Probability;

//--------------------------------------------------------------------
// TBaseLikelihoods (can also be used as haploid genotype likelihoods)
//--------------------------------------------------------------------
TBaseLikelihoods::TBaseLikelihoods(const genometools::Base & trueBase, const Probability & error){
	setFromError(trueBase, error);
};

void TBaseLikelihoods::setFromError(const genometools::Base & trueBase, const Probability & error){
	set(error / 3.0);
	_data[genometools::index(trueBase)] = error.complement();
};

//--------------------------------------------------------------------
// TBaseCounts
//--------------------------------------------------------------------
uint8_t TBaseCounts::numAlleles() const{
	using genometools::Base;
	using genometools::index;
	uint8_t n = 0;
	if(_data[index(Base::A)] > 0) ++n;
	if(_data[index(Base::C)] > 0) ++n;
	if(_data[index(Base::G)] > 0) ++n;
	if(_data[index(Base::T)] > 0) ++n;
	return n;
};

void TBaseCounts::fillFrequencies(TBaseProbabilities & freq){
	using genometools::Base;
	using genometools::index;
	double tot = _data[index(Base::A)] + _data[index(Base::C)] + _data[index(Base::G)] + _data[index(Base::T)];
	freq[Base::A] = _data[index(Base::A)] / tot;
	freq[Base::C] = _data[index(Base::C)] / tot;
	freq[Base::G] = _data[index(Base::G)] / tot;
	freq[Base::T] = _data[index(Base::T)] / tot;
};

void TBaseCounts::_fillCumulativeFrequencies(std::array<double, 4> & freq){
	using genometools::Base;
	using genometools::index;
	double tot = sum();
	freq[index(Base::A)] = _data[index(Base::A)] / tot;
	freq[index(Base::C)] = freq[index(Base::A)] + _data[index(Base::C)] / tot;
	freq[index(Base::G)] = freq[index(Base::C)] + _data[index(Base::G)] / tot;
	freq[index(Base::T)] = 1.0;
};

void TBaseCounts::downsample(uint32_t max, coretools::TRandomGenerator & RandomGenerator){
	std::array<double, 4> cumulFreqs;
	std::array<uint32_t, 5> newCounts{};

	for(uint32_t i=0; i<max; ++i){
		_fillCumulativeFrequencies(cumulFreqs);
		uint8_t b = RandomGenerator.pickOne(cumulFreqs);
		++newCounts[b];
		--_data[b];
	}

	//set counts
	_data = newCounts;
};

//--------------------------------------------------------------------
// TBaseData
//--------------------------------------------------------------------
void TBaseData::operator+=(const TBaseProbabilities & probs){
	using genometools::Base;
	using genometools::index;
	_data[index(Base::A)] += probs[Base::A].get();
	_data[index(Base::C)] += probs[Base::C].get();
	_data[index(Base::G)] += probs[Base::G].get();
	_data[index(Base::T)] += probs[Base::T].get();
};

TBaseProbabilities TBaseData::asFrequencies(){
	using genometools::Base;
	using genometools::index;
	TBaseProbabilities freq;
	double tot = sum();
	freq[Base::A] = _data[index(Base::A)] / tot;
	freq[Base::C] = _data[index(Base::C)] / tot;
	freq[Base::G] = _data[index(Base::G)] / tot;
	freq[Base::T] = _data[index(Base::T)] / tot;

	return freq;
};

//--------------------------------------------------------------------
// TGenotypeProbability_base
//--------------------------------------------------------------------
void TGenotypeProbability_base::operator=(const TGenotypeData & Other){
	using GT = genometools::Genotype;
	using genometools::index;
	_data[index(GT::AA)] = Other[GT::AA];
	_data[index(GT::AC)] = Other[GT::AC];
	_data[index(GT::AG)] = Other[GT::AG];
	_data[index(GT::AT)] = Other[GT::AT];
	_data[index(GT::CC)] = Other[GT::CC];
	_data[index(GT::CG)] = Other[GT::CG];
	_data[index(GT::CT)] = Other[GT::CT];
	_data[index(GT::GG)] = Other[GT::GG];
	_data[index(GT::GT)] = Other[GT::GT];
	_data[index(GT::TT)] = Other[GT::TT];
};

//--------------------------------------------------------------------
// TGenotypeLikelihoods
//--------------------------------------------------------------------
TGenotypeLikelihoods::TGenotypeLikelihoods(){
	reset();
};

void TGenotypeLikelihoods::fill(const std::vector<TBaseLikelihoods> & bases){
	fill(bases, bases.size());
};

void TGenotypeLikelihoods::fill(const std::vector<TBaseLikelihoods> & bases, const size_t size){
	using genometools::Base;
	using GT = genometools::Genotype;
	using genometools::index;
	//allows for vector to be longer than what is to be used
	//do in log if depth is high
	if(bases.size() > 50){
		//initialize tmp to zero
		std::array<double, 10> tmp{};

		//add to log genotype data
		for(size_t i=0; i<size; ++i){
		    tmp[index(GT::AA)] += log(bases[i][Base::A]);
		    tmp[index(GT::AC)] += log(0.5 * (bases[i][Base::A] + bases[i][Base::C]));
		    tmp[index(GT::AG)] += log(0.5 * (bases[i][Base::A] + bases[i][Base::G]));
		    tmp[index(GT::AT)] += log(0.5 * (bases[i][Base::A] + bases[i][Base::T]));
		    tmp[index(GT::CC)] += log(bases[i][Base::C].get());
		    tmp[index(GT::CG)] += log(0.5 * (bases[i][Base::C] + bases[i][Base::G]));
		    tmp[index(GT::CT)] += log(0.5 * (bases[i][Base::C] + bases[i][Base::T]));
		    tmp[index(GT::GG)] += log(bases[i][Base::G].get());
		    tmp[index(GT::GT)] += log(0.5 * (bases[i][Base::G] + bases[i][Base::T]));
		    tmp[index(GT::TT)] += log(bases[i][Base::T].get());
		}

		//standardize and de-log
		double max = *std::max_element(tmp.begin(), tmp.end());
		_data[index(GT::AA)] = exp(tmp[index(GT::AA)] - max);
		_data[index(GT::AC)] = exp(tmp[index(GT::AC)] - max);
		_data[index(GT::AG)] = exp(tmp[index(GT::AG)] - max);
		_data[index(GT::AT)] = exp(tmp[index(GT::AT)] - max);
		_data[index(GT::CC)] = exp(tmp[index(GT::CC)] - max);
		_data[index(GT::CG)] = exp(tmp[index(GT::CG)] - max);
		_data[index(GT::CT)] = exp(tmp[index(GT::CT)] - max);
		_data[index(GT::GG)] = exp(tmp[index(GT::GG)] - max);
		_data[index(GT::GT)] = exp(tmp[index(GT::GT)] - max);
		_data[index(GT::TT)] = exp(tmp[index(GT::TT)] - max);
	} else { //on natural scale
		//initialize tmp to 1.0
		TGenotypeData tmp(1.0);

		for(size_t i=0; i<size; ++i){
			tmp[GT::AA] *= bases[i][Base::A];
			tmp[GT::AC] *= 0.5 * (bases[i][Base::A] + bases[i][Base::C]);
			tmp[GT::AG] *= 0.5 * (bases[i][Base::A] + bases[i][Base::G]);
			tmp[GT::AT] *= 0.5 * (bases[i][Base::A] + bases[i][Base::T]);
			tmp[GT::CC] *= bases[i][Base::C].get();
			tmp[GT::CG] *= 0.5 * (bases[i][Base::C] + bases[i][Base::G]);
			tmp[GT::CT] *= 0.5 * (bases[i][Base::C] + bases[i][Base::T]);
			tmp[GT::GG] *= bases[i][Base::G].get();
			tmp[GT::GT] *= 0.5 * (bases[i][Base::G] + bases[i][Base::T]);
			tmp[GT::TT] *= bases[i][Base::T];
		}
		*this = tmp;
	}
};

void TGenotypeLikelihoods::addNames(std::vector<std::string> & vec) const{
	using GT = genometools::Genotype;
	using genometools::index;
	for(auto g = GT::AA; g < GT::NN; ++g){
		vec.push_back( "P(D|" + genometools::toString(g) + ")");
	}
};

//--------------------------------------------------------------------
// TGenotypeProbabilities
//--------------------------------------------------------------------
void TGenotypeProbabilities::fillBayesian(const TGenotypeLikelihoods & likelihoods, const TGenotypeProbabilities & prior){
	using GT = genometools::Genotype;
	using genometools::index;
	//calculate normalized genotype probabilities according to Bayes rule
	_data[index(GT::AA)] = likelihoods[GT::AA] * prior[GT::AA];
	_data[index(GT::AC)] = likelihoods[GT::AC] * prior[GT::AC];
	_data[index(GT::AG)] = likelihoods[GT::AG] * prior[GT::AG];
	_data[index(GT::AT)] = likelihoods[GT::AT] * prior[GT::AT];
	_data[index(GT::CC)] = likelihoods[GT::CC] * prior[GT::CC];
	_data[index(GT::CG)] = likelihoods[GT::CG] * prior[GT::CG];
	_data[index(GT::CT)] = likelihoods[GT::CT] * prior[GT::CT];
	_data[index(GT::GG)] = likelihoods[GT::GG] * prior[GT::GG];
	_data[index(GT::GT)] = likelihoods[GT::GT] * prior[GT::GT];
	_data[index(GT::TT)] = likelihoods[GT::TT] * prior[GT::TT];

	normalize();
};

Probability TGenotypeProbabilities::probHomozygous(){
	using GT = genometools::Genotype;
	using genometools::index;
	return _data[index(GT::AA)] + _data[index(GT::CC)] + _data[index(GT::GG)] + _data[index(GT::TT)];
};

Probability TGenotypeProbabilities::probHeterozygous(){
	return probHomozygous().complement();
};

//--------------------------------------------------------------------
// TGenotypeLikelihoodsHaploid
//--------------------------------------------------------------------
void TGenotypeLikelihoodsHaploid::reset(){
	using GT = genometools::Genotype;
	using genometools::index;
	//initialize to 1.0
	_data[index(GT::AA)] = 1.0; _data[index(GT::CC)] = 1.0; _data[index(GT::GG)] = 1.0; _data[index(GT::TT)] = 1.0;

	//initialize het to minimum
	_data[index(GT::AC)] = _MINLIKELIHOODVALUE; _data[index(GT::AG)] = _MINLIKELIHOODVALUE; _data[index(GT::AT)] = _MINLIKELIHOODVALUE;
	_data[index(GT::CG)] = _MINLIKELIHOODVALUE; _data[index(GT::CT)] = _MINLIKELIHOODVALUE; _data[index(GT::GT)] = _MINLIKELIHOODVALUE;
};

void TGenotypeLikelihoodsHaploid::fill(const std::vector<TBaseLikelihoods> & bases, const size_t size){
	using GT = genometools::Genotype;
	using genometools::index;
	using genometools::Base;
	//allows for vector to be longer than what is to be used
	//initialize het to minimum
	_data[index(GT::AC)] = _MINLIKELIHOODVALUE; _data[index(GT::AG)] = _MINLIKELIHOODVALUE; _data[index(GT::AT)] = _MINLIKELIHOODVALUE;
	_data[index(GT::CG)] = _MINLIKELIHOODVALUE; _data[index(GT::CT)] = _MINLIKELIHOODVALUE; _data[index(GT::GT)] = _MINLIKELIHOODVALUE;

	//do in log if depth is high
	if(bases.size() > 50){
		//initialize tmp to zero
		std::array<double, 4> tmp{};

		//add to log genotype data
		for(size_t i=0; i<size; ++i){
			tmp[0] += log(bases[i][Base::A]);
			tmp[1] += log(bases[i][Base::C]);
			tmp[2] += log(bases[i][Base::G]);
			tmp[3] += log(bases[i][Base::T]);
		}

		//find max
		double max = *std::max_element(tmp.begin(), tmp.end());

		//standardize and de-log
		_data[index(GT::AA)] = exp(tmp[0] - max);
		_data[index(GT::CC)] = exp(tmp[1] - max);
		_data[index(GT::GG)] = exp(tmp[2] - max);
		_data[index(GT::TT)] = exp(tmp[3] - max);
	} else { //on natural scale
		//initialize
		_data[index(GT::AA)] = 1.0; _data[index(GT::CC)] = 1.0; _data[index(GT::GG)] = 1.0; _data[index(GT::TT)] = 1.0;

		for(size_t i=0; i<size; ++i){
			_data[index(GT::AA)] *= bases[i][Base::A];
			_data[index(GT::CC)] *= bases[i][Base::C];
			_data[index(GT::GG)] *= bases[i][Base::G];
			_data[index(GT::TT)] *= bases[i][Base::T];
		}
	}
};

//--------------------------------------------------------------------
// TGenotypeData
//--------------------------------------------------------------------
void TGenotypeData::operator+=(const TGenotypeProbability_base & other){
	using GT = genometools::Genotype;
	using genometools::index;
	_data[index(GT::AA)] += other[GT::AA].get();
	_data[index(GT::AC)] += other[GT::AC].get();
	_data[index(GT::AG)] += other[GT::AG].get();
	_data[index(GT::AT)] += other[GT::AT].get();
	_data[index(GT::CC)] += other[GT::CC].get();
	_data[index(GT::CG)] += other[GT::CG].get();
	_data[index(GT::CT)] += other[GT::CT].get();
	_data[index(GT::GG)] += other[GT::GG].get();
	_data[index(GT::GT)] += other[GT::GT].get();
	_data[index(GT::TT)] += other[GT::TT].get();
};

}; // end namespace


