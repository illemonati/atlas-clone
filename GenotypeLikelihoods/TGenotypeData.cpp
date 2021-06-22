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
TBaseLikelihoods::TBaseLikelihoods(const BAM::Base & trueBase, const Probability & error){
	setFromError(trueBase, error);
};

void TBaseLikelihoods::setFromError(const BAM::Base & trueBase, const Probability & error){
	set(error / 3.0);
	_data[trueBase.get()] = error.complement();
};

//--------------------------------------------------------------------
// TBaseCounts
//--------------------------------------------------------------------
uint8_t TBaseCounts::numAlleles() const{
	uint8_t n = 0;
	if(_data[BAM::A] > 0) ++n;
	if(_data[BAM::C] > 0) ++n;
	if(_data[BAM::G] > 0) ++n;
	if(_data[BAM::T] > 0) ++n;
	return n;
};

void TBaseCounts::fillFrequencies(TBaseProbabilities & freq){
	double tot = _data[BAM::A] + _data[BAM::C] + _data[BAM::G] + _data[BAM::T];
	freq[BAM::A] = _data[BAM::A] / tot;
	freq[BAM::C] = _data[BAM::C] / tot;
	freq[BAM::G] = _data[BAM::G] / tot;
	freq[BAM::T] = _data[BAM::T] / tot;
};

void TBaseCounts::_fillCumulativeFrequencies(std::array<double, 4> & freq){
	double tot = sum();
	freq[BAM::A] = _data[BAM::A] / tot;
	freq[BAM::C] = freq[BAM::A] + _data[BAM::C] / tot;
	freq[BAM::G] = freq[BAM::C] + _data[BAM::G] / tot;
	freq[BAM::T] = 1.0;
};

void TBaseCounts::downsample(const uint32_t & max, coretools::TRandomGenerator & RandomGenerator){
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
	_data[BAM::A] += probs[BAM::A].get();
	_data[BAM::C] += probs[BAM::C].get();
	_data[BAM::G] += probs[BAM::G].get();
	_data[BAM::T] += probs[BAM::T].get();
};

TBaseProbabilities TBaseData::asFrequencies(){
	TBaseProbabilities freq;
	double tot = sum();
	freq[BAM::A] = _data[BAM::A] / tot;
	freq[BAM::C] = _data[BAM::C] / tot;
	freq[BAM::G] = _data[BAM::G] / tot;
	freq[BAM::T] = _data[BAM::T] / tot;

	return freq;
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
			tmp[BAM::AA] += log(bases[i][BAM::A].get());
			tmp[BAM::AC] += log(0.5 * (bases[i][BAM::A].get() + bases[i][BAM::C].get()));
			tmp[BAM::AG] += log(0.5 * (bases[i][BAM::A].get() + bases[i][BAM::G].get()));
			tmp[BAM::AT] += log(0.5 * (bases[i][BAM::A].get() + bases[i][BAM::T].get()));
			tmp[BAM::CC] += log(bases[i][BAM::C].get());
			tmp[BAM::CG] += log(0.5 * (bases[i][BAM::C].get() + bases[i][BAM::G].get()));
			tmp[BAM::CT] += log(0.5 * (bases[i][BAM::C].get() + bases[i][BAM::T].get()));
			tmp[BAM::GG] += log(bases[i][BAM::G].get());
			tmp[BAM::GT] += log(0.5 * (bases[i][BAM::G].get() + bases[i][BAM::T].get()));
			tmp[BAM::TT] += log(bases[i][BAM::T].get());
		}

		//standardize and de-log
		double max = *std::max_element(tmp.begin(), tmp.end());
		_data[BAM::AA] = exp(tmp[BAM::AA] - max);
		_data[BAM::AC] = exp(tmp[BAM::AC] - max);
		_data[BAM::AG] = exp(tmp[BAM::AG] - max);
		_data[BAM::AT] = exp(tmp[BAM::AT] - max);
		_data[BAM::CC] = exp(tmp[BAM::CC] - max);
		_data[BAM::CG] = exp(tmp[BAM::CG] - max);
		_data[BAM::CT] = exp(tmp[BAM::CT] - max);
		_data[BAM::GG] = exp(tmp[BAM::GG] - max);
		_data[BAM::GT] = exp(tmp[BAM::GT] - max);
		_data[BAM::TT] = exp(tmp[BAM::TT] - max);
	} else { //on natural scale
		//initialize
		set(1.0);

		for(size_t i=0; i<size; ++i){
			_data[BAM::AA] *= bases[i][BAM::A].get();
			_data[BAM::AC] *= 0.5 * (bases[i][BAM::A].get() + bases[i][BAM::C].get());
			_data[BAM::AG] *= 0.5 * (bases[i][BAM::A].get() + bases[i][BAM::G].get());
			_data[BAM::AT] *= 0.5 * (bases[i][BAM::A].get() + bases[i][BAM::T].get());
			_data[BAM::CC] *= bases[i][BAM::C].get();
			_data[BAM::CG] *= 0.5 * (bases[i][BAM::C].get() + bases[i][BAM::G].get());
			_data[BAM::CT] *= 0.5 * (bases[i][BAM::C].get() + bases[i][BAM::T].get());
			_data[BAM::GG] *= bases[i][BAM::G].get();
			_data[BAM::GT] *= 0.5 * (bases[i][BAM::G].get() + bases[i][BAM::T].get());
			_data[BAM::TT] *= bases[i][BAM::T].get();
		}
	}
};

void TGenotypeLikelihoods::addNames(std::vector<std::string> & vec) const{
	for(uint16_t g = BAM::AA; g < BAM::NN; g++){
		vec.push_back( "P(D|" + (std::string) BAM::Genotype(static_cast<BAM::GenotypeEnum>(g)) + ")");
	}
};

//--------------------------------------------------------------------
// TGenotypeProbabilities
//--------------------------------------------------------------------
void TGenotypeProbabilities::fillBayesian(const TGenotypeLikelihoods & likelihoods, const TGenotypeProbabilities & prior){
	//calculate normalized genotype probabilities according to Bayes rule
	_data[BAM::AA] = likelihoods[BAM::AA] * prior[BAM::AA];
	_data[BAM::AC] = likelihoods[BAM::AC] * prior[BAM::AC];
	_data[BAM::AG] = likelihoods[BAM::AG] * prior[BAM::AG];
	_data[BAM::AT] = likelihoods[BAM::AT] * prior[BAM::AT];
	_data[BAM::CC] = likelihoods[BAM::CC] * prior[BAM::CC];
	_data[BAM::CG] = likelihoods[BAM::CG] * prior[BAM::CG];
	_data[BAM::CT] = likelihoods[BAM::CT] * prior[BAM::CT];
	_data[BAM::GG] = likelihoods[BAM::GG] * prior[BAM::GG];
	_data[BAM::GT] = likelihoods[BAM::GT] * prior[BAM::GT];
	_data[BAM::TT] = likelihoods[BAM::TT] * prior[BAM::TT];

	normalize();
};

Probability TGenotypeProbabilities::probHomozygous(){
	return _data[BAM::AA] + _data[BAM::CC] + _data[BAM::GG] + _data[BAM::TT];
};

Probability TGenotypeProbabilities::probHeterozygous(){
	return probHomozygous().complement();
};

//--------------------------------------------------------------------
// TGenotypeLikelihoodsHaploid
//--------------------------------------------------------------------
void TGenotypeLikelihoodsHaploid::reset(){
	//initialize to 1.0
	_data[BAM::AA] = 1.0; _data[BAM::CC] = 1.0; _data[BAM::GG] = 1.0; _data[BAM::TT] = 1.0;

	//initialize het to minimum
	_data[BAM::AC] = _MINLIKELIHOODVALUE; _data[BAM::AG] = _MINLIKELIHOODVALUE; _data[BAM::AT] = _MINLIKELIHOODVALUE;
	_data[BAM::CG] = _MINLIKELIHOODVALUE; _data[BAM::CT] = _MINLIKELIHOODVALUE; _data[BAM::GT] = _MINLIKELIHOODVALUE;
};

void TGenotypeLikelihoodsHaploid::fill(const std::vector<TBaseData> & bases, const size_t size){
	//allows for vector to be longer than what is to be used
	//initialize het to minimum
	_data[BAM::AC] = _MINLIKELIHOODVALUE; _data[BAM::AG] = _MINLIKELIHOODVALUE; _data[BAM::AT] = _MINLIKELIHOODVALUE;
	_data[BAM::CG] = _MINLIKELIHOODVALUE; _data[BAM::CT] = _MINLIKELIHOODVALUE; _data[BAM::GT] = _MINLIKELIHOODVALUE;

	//do in log if depth is high
	if(bases.size() > 50){
		//initialize tmp to zero
		std::array<double, 4> tmp{};

		//add to log genotype data
		for(size_t i=0; i<size; ++i){
			tmp[0] += log(bases[i][BAM::A]);
			tmp[1] += log(bases[i][BAM::C]);
			tmp[2] += log(bases[i][BAM::G]);
			tmp[3] += log(bases[i][BAM::T]);
		}

		//find max
		double max = *std::max_element(tmp.begin(), tmp.end());

		//standardize and de-log
		_data[BAM::AA] = exp(tmp[0] - max);
		_data[BAM::CC] = exp(tmp[1] - max);
		_data[BAM::GG] = exp(tmp[2] - max);
		_data[BAM::TT] = exp(tmp[3] - max);
	} else { //on natural scale
		//initialize
		_data[BAM::AA] = 1.0; _data[BAM::CC] = 1.0; _data[BAM::GG] = 1.0; _data[BAM::TT] = 1.0;

		for(size_t i=0; i<size; ++i){
			_data[BAM::AA] *= bases[i][BAM::A];
			_data[BAM::CC] *= bases[i][BAM::C];
			_data[BAM::GG] *= bases[i][BAM::G];
			_data[BAM::TT] *= bases[i][BAM::T];
		}
	}
};

//--------------------------------------------------------------------
// TGenotypeData
//--------------------------------------------------------------------
void TGenotypeData::operator+=(const TGenotypeProbability_base & other){
	_data[BAM::AA] += other[BAM::AA].get();
	_data[BAM::AC] += other[BAM::AC].get();
	_data[BAM::AG] += other[BAM::AG].get();
	_data[BAM::AT] += other[BAM::AT].get();
	_data[BAM::CC] += other[BAM::CC].get();
	_data[BAM::CG] += other[BAM::CG].get();
	_data[BAM::CT] += other[BAM::CT].get();
	_data[BAM::GG] += other[BAM::GG].get();
	_data[BAM::GT] += other[BAM::GT].get();
	_data[BAM::TT] += other[BAM::TT].get();
};

}; // end namespace


