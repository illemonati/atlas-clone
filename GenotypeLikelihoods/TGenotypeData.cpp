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
// TData_base
//--------------------------------------------------------------------



//--------------------------------------------------------------------
// TBaseData (can also be used as haploid genotype likelihoods)
//--------------------------------------------------------------------
TBaseData::TBaseData(){
	reset();
};

TBaseData::TBaseData(const double val){
	set(val);
};

TBaseData::TBaseData(const BAM::Base trueBase, const double error){
	set(trueBase, error);
};

void TBaseData::operator=(const TBaseData & other){
	_data[BAM::A] = other[BAM::A];
	_data[BAM::C] = other[BAM::C];
	_data[BAM::G] = other[BAM::G];
	_data[BAM::T] = other[BAM::T];
};

void TBaseData::operator+=(const TBaseData & other){
	_data[BAM::A] += other._data[BAM::A];
	_data[BAM::C] += other._data[BAM::C];
	_data[BAM::G] += other._data[BAM::G];
	_data[BAM::T] += other._data[BAM::T];
};

void TBaseData::operator*=(const TBaseData & other){
	_data[BAM::A] *= other._data[BAM::A];
	_data[BAM::C] *= other._data[BAM::C];
	_data[BAM::G] *= other._data[BAM::G];
	_data[BAM::T] *= other._data[BAM::T];
};

void TBaseData::set(const double val){
	_data[BAM::A] = val;
	_data[BAM::G] = val;
	_data[BAM::C] = val;
	_data[BAM::T] = val;
};

void TBaseData::set(const BAM::Base trueBase, const double error){
	set(error / 3.0);
	_data[ (BAM::BaseEnum) trueBase] = 1.0 - error;
};

void TBaseData::reset(){
	set(1.0);
};

void TBaseData::add(const TBaseData & other){
	_data[BAM::A] += other._data[BAM::A];
	_data[BAM::C] += other._data[BAM::C];
	_data[BAM::G] += other._data[BAM::G];
	_data[BAM::T] += other._data[BAM::T];
};

void TBaseData::add(const Base base, const double value){
	_data[base] += value;
};

double TBaseData::sum() const{
	return _data[BAM::A] + _data[BAM::C] + _data[BAM::G] + _data[BAM::T];
};

double TBaseData::weightedSum(const TBaseData & weights) const{
	return   _data[BAM::A] * weights[BAM::A]
		   + _data[BAM::C] * weights[BAM::C]
		   + _data[BAM::G] * weights[BAM::G]
	       + _data[BAM::T] * weights[BAM::T];
};

void TBaseData::normalize(){
	double tot = sum();
	_data[BAM::A] /= tot;
	_data[BAM::C] /= tot;
	_data[BAM::G] /= tot;
	_data[BAM::T] /= tot;
};

std::ostream& operator<<(std::ostream& os, const TBaseData & baseData){
	os << "A: " << baseData[BAM::A] << ", C: " << baseData[BAM::C] << ", G: " << baseData[BAM::G] << ", T: " << baseData[BAM::T];
	return os;
};

//--------------------------------------------------------------------
// TBaseCounts
//--------------------------------------------------------------------
TBaseCounts::TBaseCounts(){
	reset();
};

void TBaseCounts::reset(){
	_data[BAM::A] = 0;
	_data[BAM::C] = 0;
	_data[BAM::G] = 0;
	_data[BAM::T] = 0;
	_data[BAM::N] = 0;
};

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

void TBaseCounts::fillCumulativeFrequencies(TBaseData & freq){
	double tot = _data[BAM::A] + _data[BAM::C] + _data[BAM::G] + _data[BAM::T];
	freq[BAM::A] = _data[BAM::A] / tot;
	freq[BAM::C] = freq[BAM::A] + _data[BAM::C] / tot;
	freq[BAM::G] = freq[BAM::C] + _data[BAM::G] / tot;
	freq[BAM::T] = 1.0;
};

void TBaseCounts::downsample(const uint32_t & max, TRandomGenerator & RandomGenerator){
	TBaseData probs;
	TBaseCounts newCounts;

	for(uint32_t i=0; i<max; ++i){
		fillCumulativeFrequencies(probs);
		uint8_t b = RandomGenerator.pickOne(probs);
		++newCounts[b];
		--_data[b];
	}

	//set counts
	_data[BAM::A] = newCounts[BAM::A];
	_data[BAM::C] = newCounts[BAM::C];
	_data[BAM::G] = newCounts[BAM::G];
	_data[BAM::T] = newCounts[BAM::T];
	_data[N] = 0;
};

//--------------------------------------------------------------------
// TGenotypeData
//--------------------------------------------------------------------

void TGenotypeData::operator=(const TGenotypeData & other){
	_data[BAM::AA] = other[BAM::AA];
	_data[BAM::AC] = other[BAM::AC];
	_data[BAM::AG] = other[BAM::AG];
	_data[BAM::AT] = other[BAM::AT];
	_data[BAM::CC] = other[BAM::CC];
	_data[BAM::CG] = other[BAM::CG];
	_data[BAM::CT] = other[BAM::CT];
	_data[BAM::GG] = other[BAM::GG];
	_data[BAM::GT] = other[BAM::GT];
	_data[BAM::TT] = other[BAM::TT];
};

void TGenotypeData::set(const double val){
	_data[BAM::AA] = val;
	_data[BAM::AC] = val;
	_data[BAM::AG] = val;
	_data[BAM::AT] = val;
	_data[BAM::CC] = val;
	_data[BAM::CG] = val;
	_data[BAM::CT] = val;
	_data[BAM::GG] = val;
	_data[BAM::GT] = val;
	_data[BAM::TT] = val;
};

void TGenotypeData::reset(){
	set(1.0);
};

void TGenotypeData::add(const TGenotypeData & other){
	_data[BAM::AA] += other[BAM::AA];
	_data[BAM::AC] += other[BAM::AC];
	_data[BAM::AG] += other[BAM::AG];
	_data[BAM::AT] += other[BAM::AT];
	_data[BAM::CC] += other[BAM::CC];
	_data[BAM::CG] += other[BAM::CG];
	_data[BAM::CT] += other[BAM::CT];
	_data[BAM::GG] += other[BAM::GG];
	_data[BAM::GT] += other[BAM::GT];
	_data[BAM::TT] += other[BAM::TT];
};

double TGenotypeData::weightedSum(const TGenotypeData & weights){
	return _data[BAM::AA] * weights[BAM::AA]
			+ _data[BAM::AC] * weights[BAM::AC]
		   	+ _data[BAM::AG] * weights[BAM::AG]
			+ _data[BAM::AT] * weights[BAM::AT]
			+ _data[BAM::CC] * weights[BAM::CC]
			+ _data[BAM::CG] * weights[BAM::CG]
			+ _data[BAM::CT] * weights[BAM::CT]
			+ _data[BAM::GG] * weights[BAM::GG]
			+ _data[BAM::GT] * weights[BAM::GT]
			+ _data[BAM::TT] * weights[BAM::TT];
};

void TGenotypeData::normalize(){
	double theSum = sum();

	_data[BAM::AA] = _data[BAM::AA] / theSum;
	_data[BAM::AC] = _data[BAM::AC] / theSum;
	_data[BAM::AG] = _data[BAM::AG] / theSum;
	_data[BAM::AT] = _data[BAM::AT] / theSum;
	_data[BAM::CC] = _data[BAM::CC] / theSum;
	_data[BAM::CG] = _data[BAM::CG] / theSum;
	_data[BAM::CT] = _data[BAM::CT] / theSum;
	_data[BAM::GG] = _data[BAM::GG] / theSum;
	_data[BAM::GT] = _data[BAM::GT] / theSum;
	_data[BAM::TT] = _data[BAM::TT] / theSum;
};

void TGenotypeData::addNames(std::vector<std::string> & vec) const{
	for(uint16_t g = BAM::AA; g < BAM::NN; g++){
		vec.push_back( (std::string) Genotype(static_cast<BAM::GenotypeEnum>(g)));
	}
};

void TGenotypeData::write(TOutputFile & out) const{
	out << _data[BAM::AA];
	out << _data[BAM::AC];
	out << _data[BAM::AG];
	out << _data[BAM::AT];
	out << _data[BAM::CC];
	out << _data[BAM::CG];
	out << _data[BAM::CT];
	out << _data[BAM::GG];
	out << _data[BAM::GT];
	out << _data[BAM::TT];
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
		double max = *std::max_element(&_data[BAM::AA], &_data[BAM::TT]);
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
// TGenotypeLikelihoodsHaploid
//--------------------------------------------------------------------
TGenotypeLikelihoodsHaploid::TGenotypeLikelihoodsHaploid(){
	reset();
};

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

double TGenotypeLikelihoodsHaploid::weightedSum(const TGenotypeData & weights){
	return _data[BAM::AA] * weights[BAM::AA]
			+ _data[BAM::CC] * weights[BAM::CC]
			+ _data[BAM::GG] * weights[BAM::GG]
			+ _data[BAM::TT] * weights[BAM::TT];
};

//--------------------------------------------------------------------
// TGenotypeProbabilities
//--------------------------------------------------------------------
TGenotypeProbabilities::TGenotypeProbabilities(){
	reset();
};

void TGenotypeProbabilities::reset(){
	set(0.1);
};

void TGenotypeProbabilities::fill(const TGenotypeData & likelihoods, const TGenotypeData & prior){
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

void TGenotypeProbabilities::addNames(std::vector<std::string> & vec) const{
	for(uint16_t g = BAM::AA; g < BAM::NN; g++){
		vec.push_back( "P(" + (std::string) Genotype(static_cast<BAM::GenotypeEnum>(g)) + "|D)");
	}
};

}; // end namespace


