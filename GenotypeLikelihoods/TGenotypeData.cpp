/*
 * TGenotypeData.cpp
 *
 *  Created on: May 15, 2020
 *      Author: phaentu
 */

#include "TGenotypeData.h"

namespace GenotypeLikelihoods{

//Note: for speed, most loops are unrolled


//--------------------------------------------------------------------
// TBaseLikelihoods (can also be used as haploid genotype likelihoods)
//--------------------------------------------------------------------
TBaseLikelihoods::TBaseLikelihoods(const BAM::Base & trueBase, const Probability & error){
	setFromError(trueBase, error);
};

void TBaseLikelihoods::setFromError(const BAM::Base trueBase, const Probability & error){
	set(error / 3.0);
	_data[ (BAM::BaseEnum) trueBase] = 1.0 - error;
};

void TBaseLikelihoods::reset(){
	set(Probability(1.0));
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

void TBaseCounts::fillFrequencies(TBaseData & freq){
	double tot = _data[BAM::A] + _data[BAM::C] + _data[BAM::G] + _data[BAM::T];
	freq[BAM::A] = _data[BAM::A] / tot;
	freq[BAM::C] = _data[BAM::C] / tot;
	freq[BAM::G] = _data[BAM::G] / tot;
	freq[BAM::T] = _data[BAM::T] / tot;
};

void TBaseCounts::fillCumulativeFrequencies(TBaseProbabilities & freq){
	double tot = sum();
	freq[BAM::A] = _data[BAM::A] / tot;
	freq[BAM::C] = freq[BAM::A] + _data[BAM::C] / tot;
	freq[BAM::G] = freq[BAM::C] + _data[BAM::G] / tot;
	freq[BAM::T] = 1.0;
};

void TBaseCounts::downsample(const uint32_t & max, TRandomGenerator & RandomGenerator){
	TBaseProbabilities probs;
	TBaseCounts newCounts;

	for(uint32_t i=0; i<max; ++i){
		fillCumulativeFrequencies(probs);
		uint8_t b = RandomGenerator.pickOne(probs);
		++newCounts[b];
		--_data[b];
	}

	//set counts
	*this = newCounts;
};

//--------------------------------------------------------------------
// TGenotypeLikelihoods
//--------------------------------------------------------------------
TGenotypeLikelihoods::TGenotypeLikelihoods(){
	reset();
};

void TGenotypeLikelihoods::fill(const std::vector<TBaseData> & bases){
	fill(bases, bases.size());
};

void TGenotypeLikelihoods::fill(const std::vector<TBaseData> & bases, const size_t size){
	//allows for vector to be longer than what is to be used
	//do in log if depth is high
	if(bases.size() > 50){
		//initialize
		set(0.0);

		//add to log genotype data
		for(size_t i=0; i<size; ++i){
			_data[BAM::AA] += log(bases[i][BAM::A]);
			_data[BAM::AC] += log(0.5*bases[i][BAM::A] + 0.5*bases[i][BAM::C]);
			_data[BAM::AG] += log(0.5*bases[i][BAM::A] + 0.5*bases[i][BAM::G]);
			_data[BAM::AT] += log(0.5*bases[i][BAM::A] + 0.5*bases[i][BAM::T]);
			_data[BAM::CC] += log(bases[i][BAM::C]);
			_data[BAM::CG] += log(0.5*bases[i][BAM::C] + 0.5*bases[i][BAM::G]);
			_data[BAM::CT] += log(0.5*bases[i][BAM::C] + 0.5*bases[i][BAM::T]);
			_data[BAM::GG] += log(bases[i][BAM::G]);
			_data[BAM::GT] += log(0.5*bases[i][BAM::G] + 0.5*bases[i][BAM::T]);
			_data[BAM::TT] += log(bases[i][BAM::T]);
		}

		//standardize and de-log
		double max = *std::max_element(_data.begin(), _data.end());
		_data[BAM::AA] = exp(_data[BAM::AA] - max);
		_data[BAM::AC] = exp(_data[BAM::AC] - max);
		_data[BAM::AG] = exp(_data[BAM::AG] - max);
		_data[BAM::AT] = exp(_data[BAM::AT] - max);
		_data[BAM::CC] = exp(_data[BAM::CC] - max);
		_data[BAM::CG] = exp(_data[BAM::CG] - max);
		_data[BAM::CT] = exp(_data[BAM::CT] - max);
		_data[BAM::GG] = exp(_data[BAM::GG] - max);
		_data[BAM::GT] = exp(_data[BAM::GT] - max);
		_data[BAM::TT] = exp(_data[BAM::TT] - max);
	} else { //on natural scale
		//initialize
		set(1.0);

		for(size_t i=0; i<size; ++i){
			_data[BAM::AA] *= bases[i][BAM::A];
			_data[BAM::AC] *= 0.5*bases[i][BAM::A] + 0.5*bases[i][BAM::C];
			_data[BAM::AG] *= 0.5*bases[i][BAM::A] + 0.5*bases[i][BAM::G];
			_data[BAM::AT] *= 0.5*bases[i][BAM::A] + 0.5*bases[i][BAM::T];
			_data[BAM::CC] *= bases[i][BAM::C];
			_data[BAM::CG] *= 0.5*bases[i][BAM::C] + 0.5*bases[i][BAM::G];
			_data[BAM::CT] *= 0.5*bases[i][BAM::C] + 0.5*bases[i][BAM::T];
			_data[BAM::GG] *= bases[i][BAM::G];
			_data[BAM::GT] *= 0.5*bases[i][BAM::G] + 0.5*bases[i][BAM::T];
			_data[BAM::TT] *= bases[i][BAM::T];
		}
	}
};

void TGenotypeLikelihoods::addNames(std::vector<std::string> & vec) const{
	for(uint16_t g = BAM::AA; g < BAM::NN; g++){
		vec.push_back( "P(D|" + (std::string) Genotype(static_cast<BAM::GenotypeEnum>(g)) + ")");
	}
};

//--------------------------------------------------------------------
// TGenotypeProbabilities
//--------------------------------------------------------------------
void TGenotypeProbabilities::fill(const TGenotypeLikelihoods & likelihoods, const TGenotypeData & prior){
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

double TGenotypeProbabilities::probHomozygous(){
	return _data[BAM::AA] + _data[BAM::CC] + _data[BAM::GG] + _data[BAM::TT];
};

double TGenotypeProbabilities::probHeterozygous(){
	return 1.0 - _data[BAM::AA] - _data[BAM::CC] - _data[BAM::GG] - _data[BAM::TT];
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
		//initialize
		_data[BAM::AA] = 0.0; _data[BAM::CC] = 0.0; _data[BAM::GG] = 0.0; _data[BAM::TT] = 0.0;


		//add to log genotype data
		for(size_t i=0; i<size; ++i){
			_data[BAM::AA] += log(bases[i][BAM::A]);
			_data[BAM::CC] += log(bases[i][BAM::C]);
			_data[BAM::GG] += log(bases[i][BAM::G]);
			_data[BAM::TT] += log(bases[i][BAM::T]);
		}

		//find max
		double max = _data[BAM::AA];
		if(_data[BAM::CC] > max) max = _data[BAM::CC];
		if(_data[BAM::GG] > max) max = _data[BAM::GG];
		if(_data[BAM::TT] > max) max = _data[BAM::TT];

		//standardize and de-log
		_data[BAM::AA] = exp(_data[BAM::AA] - max);
		_data[BAM::CC] = exp(_data[BAM::CC] - max);
		_data[BAM::GG] = exp(_data[BAM::GG] - max);
		_data[BAM::TT] = exp(_data[BAM::TT] - max);
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

}; // end namespace


