/*
 * TGenotypeData.cpp
 *
 *  Created on: May 15, 2020
 *      Author: phaentu
 */

#include "TGenotypeData.h"

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
	_data[trueBase.get()] = error.complement();
};

//--------------------------------------------------------------------
// TBaseCounts
//--------------------------------------------------------------------
uint8_t TBaseCounts::numAlleles() const{
	uint8_t n = 0;
	if(_data[genometools::A] > 0) ++n;
	if(_data[genometools::C] > 0) ++n;
	if(_data[genometools::G] > 0) ++n;
	if(_data[genometools::T] > 0) ++n;
	return n;
};

void TBaseCounts::fillFrequencies(TBaseProbabilities & freq){
	double tot = _data[genometools::A] + _data[genometools::C] + _data[genometools::G] + _data[genometools::T];
	freq[genometools::A] = _data[genometools::A] / tot;
	freq[genometools::C] = _data[genometools::C] / tot;
	freq[genometools::G] = _data[genometools::G] / tot;
	freq[genometools::T] = _data[genometools::T] / tot;
};

void TBaseCounts::_fillCumulativeFrequencies(std::array<double, 4> & freq){
	double tot = sum();
	freq[genometools::A] = _data[genometools::A] / tot;
	freq[genometools::C] = freq[genometools::A] + _data[genometools::C] / tot;
	freq[genometools::G] = freq[genometools::C] + _data[genometools::G] / tot;
	freq[genometools::T] = 1.0;
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
	_data[genometools::A] += probs[genometools::A].get();
	_data[genometools::C] += probs[genometools::C].get();
	_data[genometools::G] += probs[genometools::G].get();
	_data[genometools::T] += probs[genometools::T].get();
};

TBaseProbabilities TBaseData::asFrequencies(){
	TBaseProbabilities freq;
	double tot = sum();
	freq[genometools::A] = _data[genometools::A] / tot;
	freq[genometools::C] = _data[genometools::C] / tot;
	freq[genometools::G] = _data[genometools::G] / tot;
	freq[genometools::T] = _data[genometools::T] / tot;

	return freq;
};

//--------------------------------------------------------------------
// TGenotypeProbability_base
//--------------------------------------------------------------------
void TGenotypeProbability_base::operator=(const TGenotypeData & Other){
	_data[genometools::AA] = Other[genometools::AA];
	_data[genometools::AC] = Other[genometools::AC];
	_data[genometools::AG] = Other[genometools::AG];
	_data[genometools::AT] = Other[genometools::AT];
	_data[genometools::CC] = Other[genometools::CC];
	_data[genometools::CG] = Other[genometools::CG];
	_data[genometools::CT] = Other[genometools::CT];
	_data[genometools::GG] = Other[genometools::GG];
	_data[genometools::GT] = Other[genometools::GT];
	_data[genometools::TT] = Other[genometools::TT];
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
	//allows for vector to be longer than what is to be used
	//do in log if depth is high
	if(bases.size() > 50){
		//initialize tmp to zero
		std::array<double, 10> tmp{};

		//add to log genotype data
		for(size_t i=0; i<size; ++i){
			tmp[genometools::AA] += log(bases[i][genometools::A]);
			tmp[genometools::AC] += log(0.5 * (bases[i][genometools::A] + bases[i][genometools::C]));
			tmp[genometools::AG] += log(0.5 * (bases[i][genometools::A] + bases[i][genometools::G]));
			tmp[genometools::AT] += log(0.5 * (bases[i][genometools::A] + bases[i][genometools::T]));
			tmp[genometools::CC] += log(bases[i][genometools::C].get());
			tmp[genometools::CG] += log(0.5 * (bases[i][genometools::C] + bases[i][genometools::G]));
			tmp[genometools::CT] += log(0.5 * (bases[i][genometools::C] + bases[i][genometools::T]));
			tmp[genometools::GG] += log(bases[i][genometools::G].get());
			tmp[genometools::GT] += log(0.5 * (bases[i][genometools::G] + bases[i][genometools::T]));
			tmp[genometools::TT] += log(bases[i][genometools::T].get());
		}

		//standardize and de-log
		double max = *std::max_element(tmp.begin(), tmp.end());
		_data[genometools::AA] = exp(tmp[genometools::AA] - max);
		_data[genometools::AC] = exp(tmp[genometools::AC] - max);
		_data[genometools::AG] = exp(tmp[genometools::AG] - max);
		_data[genometools::AT] = exp(tmp[genometools::AT] - max);
		_data[genometools::CC] = exp(tmp[genometools::CC] - max);
		_data[genometools::CG] = exp(tmp[genometools::CG] - max);
		_data[genometools::CT] = exp(tmp[genometools::CT] - max);
		_data[genometools::GG] = exp(tmp[genometools::GG] - max);
		_data[genometools::GT] = exp(tmp[genometools::GT] - max);
		_data[genometools::TT] = exp(tmp[genometools::TT] - max);
	} else { //on natural scale
		//initialize tmp to 1.0
		TGenotypeData tmp(1.0);

		for(size_t i=0; i<size; ++i){
			tmp[genometools::AA] *= bases[i][genometools::A];
			tmp[genometools::AC] *= 0.5 * (bases[i][genometools::A] + bases[i][genometools::C]);
			tmp[genometools::AG] *= 0.5 * (bases[i][genometools::A] + bases[i][genometools::G]);
			tmp[genometools::AT] *= 0.5 * (bases[i][genometools::A] + bases[i][genometools::T]);
			tmp[genometools::CC] *= bases[i][genometools::C].get();
			tmp[genometools::CG] *= 0.5 * (bases[i][genometools::C] + bases[i][genometools::G]);
			tmp[genometools::CT] *= 0.5 * (bases[i][genometools::C] + bases[i][genometools::T]);
			tmp[genometools::GG] *= bases[i][genometools::G].get();
			tmp[genometools::GT] *= 0.5 * (bases[i][genometools::G] + bases[i][genometools::T]);
			tmp[genometools::TT] *= bases[i][genometools::T];
		}
		*this = tmp;
	}
};

void TGenotypeLikelihoods::addNames(std::vector<std::string> & vec) const{
	for(uint16_t g = genometools::AA; g < genometools::NN; g++){
		vec.push_back( "P(D|" + (std::string) genometools::Genotype(static_cast<genometools::GenotypeEnum>(g)) + ")");
	}
};

//--------------------------------------------------------------------
// TGenotypeProbabilities
//--------------------------------------------------------------------
void TGenotypeProbabilities::fillBayesian(const TGenotypeLikelihoods & likelihoods, const TGenotypeProbabilities & prior){
	//calculate normalized genotype probabilities according to Bayes rule
	_data[genometools::AA] = likelihoods[genometools::AA] * prior[genometools::AA];
	_data[genometools::AC] = likelihoods[genometools::AC] * prior[genometools::AC];
	_data[genometools::AG] = likelihoods[genometools::AG] * prior[genometools::AG];
	_data[genometools::AT] = likelihoods[genometools::AT] * prior[genometools::AT];
	_data[genometools::CC] = likelihoods[genometools::CC] * prior[genometools::CC];
	_data[genometools::CG] = likelihoods[genometools::CG] * prior[genometools::CG];
	_data[genometools::CT] = likelihoods[genometools::CT] * prior[genometools::CT];
	_data[genometools::GG] = likelihoods[genometools::GG] * prior[genometools::GG];
	_data[genometools::GT] = likelihoods[genometools::GT] * prior[genometools::GT];
	_data[genometools::TT] = likelihoods[genometools::TT] * prior[genometools::TT];

	normalize();
};

Probability TGenotypeProbabilities::probHomozygous(){
	return _data[genometools::AA] + _data[genometools::CC] + _data[genometools::GG] + _data[genometools::TT];
};

Probability TGenotypeProbabilities::probHeterozygous(){
	return probHomozygous().complement();
};

//--------------------------------------------------------------------
// TGenotypeLikelihoodsHaploid
//--------------------------------------------------------------------
void TGenotypeLikelihoodsHaploid::reset(){
	//initialize to 1.0
	_data[genometools::AA] = 1.0; _data[genometools::CC] = 1.0; _data[genometools::GG] = 1.0; _data[genometools::TT] = 1.0;

	//initialize het to minimum
	_data[genometools::AC] = _MINLIKELIHOODVALUE; _data[genometools::AG] = _MINLIKELIHOODVALUE; _data[genometools::AT] = _MINLIKELIHOODVALUE;
	_data[genometools::CG] = _MINLIKELIHOODVALUE; _data[genometools::CT] = _MINLIKELIHOODVALUE; _data[genometools::GT] = _MINLIKELIHOODVALUE;
};

void TGenotypeLikelihoodsHaploid::fill(const std::vector<TBaseLikelihoods> & bases, const size_t size){
	//allows for vector to be longer than what is to be used
	//initialize het to minimum
	_data[genometools::AC] = _MINLIKELIHOODVALUE; _data[genometools::AG] = _MINLIKELIHOODVALUE; _data[genometools::AT] = _MINLIKELIHOODVALUE;
	_data[genometools::CG] = _MINLIKELIHOODVALUE; _data[genometools::CT] = _MINLIKELIHOODVALUE; _data[genometools::GT] = _MINLIKELIHOODVALUE;

	//do in log if depth is high
	if(bases.size() > 50){
		//initialize tmp to zero
		std::array<double, 4> tmp{};

		//add to log genotype data
		for(size_t i=0; i<size; ++i){
			tmp[0] += log(bases[i][genometools::A]);
			tmp[1] += log(bases[i][genometools::C]);
			tmp[2] += log(bases[i][genometools::G]);
			tmp[3] += log(bases[i][genometools::T]);
		}

		//find max
		double max = *std::max_element(tmp.begin(), tmp.end());

		//standardize and de-log
		_data[genometools::AA] = exp(tmp[0] - max);
		_data[genometools::CC] = exp(tmp[1] - max);
		_data[genometools::GG] = exp(tmp[2] - max);
		_data[genometools::TT] = exp(tmp[3] - max);
	} else { //on natural scale
		//initialize
		_data[genometools::AA] = 1.0; _data[genometools::CC] = 1.0; _data[genometools::GG] = 1.0; _data[genometools::TT] = 1.0;

		for(size_t i=0; i<size; ++i){
			_data[genometools::AA] *= bases[i][genometools::A];
			_data[genometools::CC] *= bases[i][genometools::C];
			_data[genometools::GG] *= bases[i][genometools::G];
			_data[genometools::TT] *= bases[i][genometools::T];
		}
	}
};

//--------------------------------------------------------------------
// TGenotypeData
//--------------------------------------------------------------------
void TGenotypeData::operator+=(const TGenotypeProbability_base & other){
	_data[genometools::AA] += other[genometools::AA].get();
	_data[genometools::AC] += other[genometools::AC].get();
	_data[genometools::AG] += other[genometools::AG].get();
	_data[genometools::AT] += other[genometools::AT].get();
	_data[genometools::CC] += other[genometools::CC].get();
	_data[genometools::CG] += other[genometools::CG].get();
	_data[genometools::CT] += other[genometools::CT].get();
	_data[genometools::GG] += other[genometools::GG].get();
	_data[genometools::GT] += other[genometools::GT].get();
	_data[genometools::TT] += other[genometools::TT].get();
};

}; // end namespace


