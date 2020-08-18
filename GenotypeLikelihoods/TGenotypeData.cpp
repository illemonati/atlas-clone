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

TBaseData::TBaseData(const Base trueBase, const double error){
	set(trueBase, error);
};

void TBaseData::operator=(const TBaseData & other){
	_data[A] = other[A];
	_data[C] = other[C];
	_data[G] = other[G];
	_data[T] = other[T];
};

void TBaseData::operator+=(const TBaseData & other){
	_data[A] += other._data[A];
	_data[C] += other._data[C];
	_data[G] += other._data[G];
	_data[T] += other._data[T];
};

void TBaseData::operator*=(const TBaseData & other){
	_data[A] *= other._data[A];
	_data[C] *= other._data[C];
	_data[G] *= other._data[G];
	_data[T] *= other._data[T];
};

void TBaseData::set(const double val){
	_data[A] = val;
	_data[G] = val;
	_data[C] = val;
	_data[T] = val;
};

void TBaseData::set(const Base trueBase, const double error){
	set(error / 3.0);
	_data[trueBase] = 1.0 - error;
};

void TBaseData::reset(){
	set(1.0);
};

void TBaseData::add(const TBaseData & other){
	_data[A] += other._data[A];
	_data[C] += other._data[C];
	_data[G] += other._data[G];
	_data[T] += other._data[T];
};

void TBaseData::add(const Base base, const double value){
	_data[base] += value;
};

double TBaseData::sum() const{
	return _data[A] + _data[C] + _data[G] + _data[T];
};

double TBaseData::weightedSum(const TBaseData & weights) const{
	return   _data[A] * weights[A]
		   + _data[C] * weights[C]
		   + _data[G] * weights[G]
	       + _data[T] * weights[T];
};

void TBaseData::normalize(){
	double tot = sum();
	_data[A] /= tot;
	_data[C] /= tot;
	_data[G] /= tot;
	_data[T] /= tot;
};

std::ostream& operator<<(std::ostream& os, const TBaseData & baseData){
	os << "A: " << baseData[A] << ", C: " << baseData[C] << ", G: " << baseData[G] << ", T: " << baseData[T];
	return os;
};

//--------------------------------------------------------------------
// TBaseCounts
//--------------------------------------------------------------------
TBaseCounts::TBaseCounts(){
	reset();
};

void TBaseCounts::reset(){
	_data[A] = 0;
	_data[C] = 0;
	_data[G] = 0;
	_data[T] = 0;
	_data[N] = 0;
};

uint8_t TBaseCounts::numAlleles() const{
	uint8_t n = 0;
	if(_data[A] > 0) ++n;
	if(_data[C] > 0) ++n;
	if(_data[G] > 0) ++n;
	if(_data[T] > 0) ++n;
	return n;
};

void TBaseCounts::fillFrequencies(TBaseData & freq){
	double tot = _data[A] + _data[C] + _data[G] + _data[T];
	freq[A] = _data[A] / tot;
	freq[C] = _data[C] / tot;
	freq[G] = _data[G] / tot;
	freq[T] = _data[T] / tot;
};

void TBaseCounts::fillCumulativeFrequencies(TBaseData & freq){
	double tot = _data[A] + _data[C] + _data[G] + _data[T];
	freq[A] = _data[A] / tot;
	freq[C] = freq[A] + _data[C] / tot;
	freq[G] = freq[C] + _data[G] / tot;
	freq[T] = 1.0;
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
	_data[A] = newCounts[A];
	_data[C] = newCounts[C];
	_data[G] = newCounts[G];
	_data[T] = newCounts[T];
	_data[N] = 0;
};

//--------------------------------------------------------------------
// TGenotypeData
//--------------------------------------------------------------------
/*
void TGenotypeData::_copyFrom(const TGenotypeData & other){
	_data[AA] = other[AA];
	_data[AC] = other[AC];
	_data[AG] = other[AG];
	_data[AT] = other[AT];
	_data[CC] = other[CC];
	_data[CG] = other[CG];
	_data[CT] = other[CT];
	_data[GG] = other[GG];
	_data[GT] = other[GT];
	_data[TT] = other[TT];
};
*/

void TGenotypeData::operator=(const TGenotypeData & other){
	_data[AA] = other[AA];
	_data[AC] = other[AC];
	_data[AG] = other[AG];
	_data[AT] = other[AT];
	_data[CC] = other[CC];
	_data[CG] = other[CG];
	_data[CT] = other[CT];
	_data[GG] = other[GG];
	_data[GT] = other[GT];
	_data[TT] = other[TT];
};

void TGenotypeData::set(const double val){
	_data[AA] = val;
	_data[AC] = val;
	_data[AG] = val;
	_data[AT] = val;
	_data[CC] = val;
	_data[CG] = val;
	_data[CT] = val;
	_data[GG] = val;
	_data[GT] = val;
	_data[TT] = val;
};

void TGenotypeData::reset(){
	set(1.0);
};

void TGenotypeData::add(const TGenotypeData & other){
	_data[AA] += other[AA];
	_data[AC] += other[AC];
	_data[AG] += other[AG];
	_data[AT] += other[AT];
	_data[CC] += other[CC];
	_data[CG] += other[CG];
	_data[CT] += other[CT];
	_data[GG] += other[GG];
	_data[GT] += other[GT];
	_data[TT] += other[TT];
};

double TGenotypeData::weightedSum(const TGenotypeData & weights){
	return _data[AA] * weights[AA]
			+ _data[AC] * weights[AC]
		   	+ _data[AG] * weights[AG]
			+ _data[AT] * weights[AT]
			+ _data[CC] * weights[CC]
			+ _data[CG] * weights[CG]
			+ _data[CT] * weights[CT]
			+ _data[GG] * weights[GG]
			+ _data[GT] * weights[GT]
			+ _data[TT] * weights[TT];
};

void TGenotypeData::normalize(){
	double theSum = sum();

	_data[AA] = _data[AA] / theSum;
	_data[AC] = _data[AC] / theSum;
	_data[AG] = _data[AG] / theSum;
	_data[AT] = _data[AT] / theSum;
	_data[CC] = _data[CC] / theSum;
	_data[CG] = _data[CG] / theSum;
	_data[CT] = _data[CT] / theSum;
	_data[GG] = _data[GG] / theSum;
	_data[GT] = _data[GT] / theSum;
	_data[TT] = _data[TT] / theSum;
};

void TGenotypeData::addNames(std::vector<std::string> & vec, const TGenotypeMap & genoMap) const{
	for(uint16_t g=0; g<genoMap.numGenotypes; ++g){
		vec.push_back(genoMap.getGenotypeString(g));
	}
};

void TGenotypeData::write(TOutputFile & out) const{
	out << _data[AA];
	out << _data[AC];
	out << _data[AG];
	out << _data[AT];
	out << _data[CC];
	out << _data[CG];
	out << _data[CT];
	out << _data[GG];
	out << _data[GT];
	out << _data[TT];
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
			_data[AA] += log(bases[i][A]);
			_data[AC] += log(0.5*bases[i][A] + 0.5*bases[i][C]);
			_data[AG] += log(0.5*bases[i][A] + 0.5*bases[i][G]);
			_data[AT] += log(0.5*bases[i][A] + 0.5*bases[i][T]);
			_data[CC] += log(bases[i][C]);
			_data[CG] += log(0.5*bases[i][C] + 0.5*bases[i][G]);
			_data[CT] += log(0.5*bases[i][C] + 0.5*bases[i][T]);
			_data[GG] += log(bases[i][G]);
			_data[GT] += log(0.5*bases[i][G] + 0.5*bases[i][T]);
			_data[TT] += log(bases[i][T]);
		}

		//standardize and de-log
		double max = *std::max_element(&_data[AA], &_data[TT]);
		_data[AA] = exp(_data[AA] - max);
		_data[AC] = exp(_data[AC] - max);
		_data[AG] = exp(_data[AG] - max);
		_data[AT] = exp(_data[AT] - max);
		_data[CC] = exp(_data[CC] - max);
		_data[CG] = exp(_data[CG] - max);
		_data[CT] = exp(_data[CT] - max);
		_data[GG] = exp(_data[GG] - max);
		_data[GT] = exp(_data[GT] - max);
		_data[TT] = exp(_data[TT] - max);
	} else { //on natural scale
		//initialize
		set(1.0);

		for(size_t i=0; i<size; ++i){
			_data[AA] *= bases[i][A];
			_data[AC] *= 0.5*bases[i][A] + 0.5*bases[i][C];
			_data[AG] *= 0.5*bases[i][A] + 0.5*bases[i][G];
			_data[AT] *= 0.5*bases[i][A] + 0.5*bases[i][T];
			_data[CC] *= bases[i][C];
			_data[CG] *= 0.5*bases[i][C] + 0.5*bases[i][G];
			_data[CT] *= 0.5*bases[i][C] + 0.5*bases[i][T];
			_data[GG] *= bases[i][G];
			_data[GT] *= 0.5*bases[i][G] + 0.5*bases[i][T];
			_data[TT] *= bases[i][T];
		}
	}
};

void TGenotypeLikelihoods::addNames(std::vector<std::string> & vec, const TGenotypeMap & genoMap) const{
	for(uint16_t g=0; g<genoMap.numGenotypes; ++g){
		vec.push_back("P(D|" + genoMap.getGenotypeString(g) + ")");
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
	_data[AA] = 1.0; _data[CC] = 1.0; _data[GG] = 1.0; _data[TT] = 1.0;

	//initialize het to minimum
	_data[AC] = _MINLIKELIHOODVALUE; _data[AG] = _MINLIKELIHOODVALUE; _data[AT] = _MINLIKELIHOODVALUE;
	_data[CG] = _MINLIKELIHOODVALUE; _data[CT] = _MINLIKELIHOODVALUE; _data[GT] = _MINLIKELIHOODVALUE;
};

void TGenotypeLikelihoodsHaploid::fill(const std::vector<TBaseData> & bases, const size_t size){
	//allows for vector to be longer than what is to be used
	//initialize het to minimum
	_data[AC] = _MINLIKELIHOODVALUE; _data[AG] = _MINLIKELIHOODVALUE; _data[AT] = _MINLIKELIHOODVALUE;
	_data[CG] = _MINLIKELIHOODVALUE; _data[CT] = _MINLIKELIHOODVALUE; _data[GT] = _MINLIKELIHOODVALUE;

	//do in log if depth is high
	if(bases.size() > 50){
		//initialize
		_data[AA] = 0.0; _data[CC] = 0.0; _data[GG] = 0.0; _data[TT] = 0.0;


		//add to log genotype data
		for(size_t i=0; i<size; ++i){
			_data[AA] += log(bases[i][A]);
			_data[CC] += log(bases[i][C]);
			_data[GG] += log(bases[i][G]);
			_data[TT] += log(bases[i][T]);
		}

		//find max
		double max = _data[AA];
		if(_data[CC] > max) max = _data[CC];
		if(_data[GG] > max) max = _data[GG];
		if(_data[TT] > max) max = _data[TT];

		//standardize and de-log
		_data[AA] = exp(_data[AA] - max);
		_data[CC] = exp(_data[CC] - max);
		_data[GG] = exp(_data[GG] - max);
		_data[TT] = exp(_data[TT] - max);
	} else { //on natural scale
		//initialize
		_data[AA] = 1.0; _data[CC] = 1.0; _data[GG] = 1.0; _data[TT] = 1.0;

		for(size_t i=0; i<size; ++i){
			_data[AA] *= bases[i][A];
			_data[CC] *= bases[i][C];
			_data[GG] *= bases[i][G];
			_data[TT] *= bases[i][T];
		}
	}
};

double TGenotypeLikelihoodsHaploid::weightedSum(const TGenotypeData & weights){
	return _data[AA] * weights[AA]
			+ _data[CC] * weights[CC]
			+ _data[GG] * weights[GG]
			+ _data[TT] * weights[TT];
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
	_data[AA] = likelihoods[AA] * prior[AA];
	_data[AC] = likelihoods[AC] * prior[AC];
	_data[AG] = likelihoods[AG] * prior[AG];
	_data[AT] = likelihoods[AT] * prior[AT];
	_data[CC] = likelihoods[CC] * prior[CC];
	_data[CG] = likelihoods[CG] * prior[CG];
	_data[CT] = likelihoods[CT] * prior[CT];
	_data[GG] = likelihoods[GG] * prior[GG];
	_data[GT] = likelihoods[GT] * prior[GT];
	_data[TT] = likelihoods[TT] * prior[TT];

	normalize();
};

double TGenotypeProbabilities::probHomozygous(){
	return _data[AA] + _data[CC] + _data[GG] + _data[TT];
};

double TGenotypeProbabilities::probHeterozygous(){
	return 1.0 - _data[AA] - _data[CC] - _data[GG] - _data[TT];
};

void TGenotypeProbabilities::addNames(std::vector<std::string> & vec, const TGenotypeMap & genoMap) const{
	for(uint16_t g=0; g<genoMap.numGenotypes; ++g){
		vec.push_back("P(" + genoMap.getGenotypeString(g) + "|D)");
	}
};

}; // end namespace


