/*
 * TGenotypeData.cpp
 *
 *  Created on: May 15, 2020
 *      Author: phaentu
 */

#include "TGenotypeData.h"

namespace GenotypeLikelihoods{

//Note: for speed, all loops are unrolled

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
	_data[A] = other.at(A);
	_data[C] = other.at(C);
	_data[G] = other.at(G);
	_data[T] = other.at(T);
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

double TBaseData::min() const{
	return *std::min_element(_data.begin(), _data.end());
};

double TBaseData::max() const{
	return *std::max_element(_data.begin(), _data.end());
};

double TBaseData::weightedSum(const TBaseData & weights) const{
	return   _data[A] * weights.at(A)
		   + _data[C] * weights.at(C)
		   + _data[G] * weights.at(G)
	       + _data[T] * weights.at(T);
};

void TBaseData::normalize(){
	double tot = sum();
	_data[A] /= tot;
	_data[C] /= tot;
	_data[G] /= tot;
	_data[T] /= tot;
};

std::ostream& operator<<(std::ostream& os, const TBaseData & baseData){
	os << "A: " << baseData.at(A) << ", C: " << baseData.at(C) << ", G: " << baseData.at(G) << ", T: " << baseData.at(T);
	return os;
};

//--------------------------------------------------------------------
// TBaseCounts
//--------------------------------------------------------------------
TBaseCounts::TBaseCounts(){
	reset();
};

void TBaseCounts::reset(){
	_counts[A] = 0;
	_counts[C] = 0;
	_counts[G] = 0;
	_counts[T] = 0;
	_counts[N] = 0;
};

uint32_t TBaseCounts::size() const{
	return _counts[A] + _counts[C] + _counts[G] + _counts[T] + _counts[N];
};

uint8_t TBaseCounts::numAlleles() const{
	uint8_t n = 0;
	if(_counts[A] > 0) ++n;
	if(_counts[C] > 0) ++n;
	if(_counts[G] > 0) ++n;
	if(_counts[T] > 0) ++n;
	return n;
};

void TBaseCounts::fillFrequencies(TBaseData & freq){
	double tot = size() - _counts[N];
	freq[A] = _counts[A] / tot;
	freq[C] = _counts[C] / tot;
	freq[G] = _counts[G] / tot;
	freq[T] = _counts[T] / tot;
};

void TBaseCounts::fillCumulativeFrequencies(TBaseData & freq){
	double tot = size() - _counts[N];
	freq[A] = _counts[A] / tot;
	freq[C] = freq[A] + _counts[C] / tot;
	freq[G] = freq[C] + _counts[G] / tot;
	freq[T] = 1.0;
};

void TBaseCounts::downsample(const uint32_t & max, TRandomGenerator & RandomGenerator){
	TBaseData probs;
	TBaseCounts newCounts;

	for(uint32_t i=0; i<max; ++i){
		fillCumulativeFrequencies(probs);
		uint8_t b = RandomGenerator.pickOne(4, probs.data());
		++newCounts[b];
		--_counts[b];
	}

	//set counts
	_counts[A] = newCounts[A];
	_counts[C] = newCounts[C];
	_counts[G] = newCounts[G];
	_counts[T] = newCounts[T];
	_counts[N] = 0;
};

//--------------------------------------------------------------------
// TGenotypeData
//--------------------------------------------------------------------

void TGenotypeData::_copyFrom(const TGenotypeData & other){
	_data[AA] = other.at(AA);
	_data[AC] = other.at(AC);
	_data[AG] = other.at(AG);
	_data[AT] = other.at(AT);
	_data[CC] = other.at(CC);
	_data[CG] = other.at(CG);
	_data[CT] = other.at(CT);
	_data[GG] = other.at(GG);
	_data[GT] = other.at(GT);
	_data[TT] = other.at(TT);
};

void TGenotypeData::operator=(const TGenotypeData & other){
	_copyFrom(other);
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
	_data[AA] += other.at(AA);
	_data[AC] += other.at(AC);
	_data[AG] += other.at(AG);
	_data[AT] += other.at(AT);
	_data[CC] += other.at(CC);
	_data[CG] += other.at(CG);
	_data[CT] += other.at(CT);
	_data[GG] += other.at(GG);
	_data[GT] += other.at(GT);
	_data[TT] += other.at(TT);
};

double TGenotypeData::sum() const{
	return _data[AA] + _data[AC] + _data[AG] + _data[AT] + _data[CC] + _data[CG] + _data[CT] + _data[GG] + _data[GT] + _data[TT];
};

double TGenotypeData::min() const{
	return *std::min_element(_data.begin(), _data.end());
};

double TGenotypeData::max() const{
	return *std::max_element(_data.begin(), _data.end());
};

double TGenotypeData::weightedSum(const TGenotypeData & weights){
	return _data[AA] * weights.at(AA)
			+ _data[AC] * weights.at(AC)
		   	+ _data[AG] * weights.at(AG)
			+ _data[AT] * weights.at(AT)
			+ _data[CC] * weights.at(CC)
			+ _data[CG] * weights.at(CG)
			+ _data[CT] * weights.at(CT)
			+ _data[GG] * weights.at(GG)
			+ _data[GT] * weights.at(GT)
			+ _data[TT] * weights.at(TT);
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

void TGenotypeLikelihoods::operator=(const TGenotypeLikelihoods & other){
	_copyFrom(other);
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
			_data[AA] += log(bases[i].at(A));
			_data[AC] += log(0.5*bases[i].at(A) + 0.5*bases[i].at(C));
			_data[AG] += log(0.5*bases[i].at(A) + 0.5*bases[i].at(G));
			_data[AT] += log(0.5*bases[i].at(A) + 0.5*bases[i].at(T));
			_data[CC] += log(bases[i].at(C));
			_data[CG] += log(0.5*bases[i].at(C) + 0.5*bases[i].at(G));
			_data[CT] += log(0.5*bases[i].at(C) + 0.5*bases[i].at(T));
			_data[GG] += log(bases[i].at(G));
			_data[GT] += log(0.5*bases[i].at(G) + 0.5*bases[i].at(T));
			_data[TT] += log(bases[i].at(T));
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
			_data[AA] *= bases[i].at(A);
			_data[AC] *= 0.5*bases[i].at(A) + 0.5*bases[i].at(C);
			_data[AG] *= 0.5*bases[i].at(A) + 0.5*bases[i].at(G);
			_data[AT] *= 0.5*bases[i].at(A) + 0.5*bases[i].at(T);
			_data[CC] *= bases[i].at(C);
			_data[CG] *= 0.5*bases[i].at(C) + 0.5*bases[i].at(G);
			_data[CT] *= 0.5*bases[i].at(C) + 0.5*bases[i].at(T);
			_data[GG] *= bases[i].at(G);
			_data[GT] *= 0.5*bases[i].at(G) + 0.5*bases[i].at(T);
			_data[TT] *= bases[i].at(T);
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
			_data[AA] += log(bases[i].at(A));
			_data[CC] += log(bases[i].at(C));
			_data[GG] += log(bases[i].at(G));
			_data[TT] += log(bases[i].at(T));
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
			_data[AA] *= bases[i].at(A);
			_data[CC] *= bases[i].at(C);
			_data[GG] *= bases[i].at(G);
			_data[TT] *= bases[i].at(T);
		}
	}
};

double TGenotypeLikelihoodsHaploid::weightedSum(const TGenotypeData & weights){
	return _data[AA] * weights.at(AA)
			+ _data[CC] * weights.at(CC)
			+ _data[GG] * weights.at(GG)
			+ _data[TT] * weights.at(TT);
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
	_data[AA] = likelihoods.at(AA) * prior.at(AA);
	_data[AC] = likelihoods.at(AC) * prior.at(AC);
	_data[AG] = likelihoods.at(AG) * prior.at(AG);
	_data[AT] = likelihoods.at(AT) * prior.at(AT);
	_data[CC] = likelihoods.at(CC) * prior.at(CC);
	_data[CG] = likelihoods.at(CG) * prior.at(CG);
	_data[CT] = likelihoods.at(CT) * prior.at(CT);
	_data[GG] = likelihoods.at(GG) * prior.at(GG);
	_data[GT] = likelihoods.at(GT) * prior.at(GT);
	_data[TT] = likelihoods.at(TT) * prior.at(TT);

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


