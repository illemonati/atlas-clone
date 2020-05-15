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
// TBaseLikelihoods (can also be used as haploid genotype likelihoods)
//--------------------------------------------------------------------
TBaseLikelihoods::TBaseLikelihoods(){
	reset();
};

void TBaseLikelihoods::operator=(const TBaseLikelihoods & other){
	likelihoods[A] = other.at(A);
	likelihoods[C] = other.at(C);
	likelihoods[G] = other.at(G);
	likelihoods[T] = other.at(T);
	likelihoods[N] = other.at(N);
};

void TBaseLikelihoods::reset(){
	likelihoods[A] = 1.0;
	likelihoods[C] = 1.0;
	likelihoods[G] = 1.0;
	likelihoods[T] = 1.0;
	likelihoods[N] = 1.0;
};

//--------------------------------------------------------------------
// TGenotypeData
//--------------------------------------------------------------------

void TGenotypeData::_copyFrom(const TGenotypeData & other){
	data[AA] = other.at(AA);
	data[AC] = other.at(AC);
	data[AG] = other.at(AG);
	data[AT] = other.at(AT);
	data[CC] = other.at(CC);
	data[CG] = other.at(CG);
	data[CT] = other.at(CT);
	data[GG] = other.at(GG);
	data[GT] = other.at(GT);
	data[TT] = other.at(TT);
};

void TGenotypeData::operator=(const TGenotypeData & other){
	_copyFrom(other);
};

void TGenotypeData::set(const double val){
	data[AA] = val;
	data[AC] = val;
	data[AG] = val;
	data[AT] = val;
	data[CC] = val;
	data[CG] = val;
	data[CT] = val;
	data[GG] = val;
	data[GT] = val;
	data[TT] = val;
};

void TGenotypeData::reset(){
	set(1.0);
};

void TGenotypeData::add(const TGenotypeData & other){
	data[AA] += other.at(AA);
	data[AC] += other.at(AC);
	data[AG] += other.at(AG);
	data[AT] += other.at(AT);
	data[CC] += other.at(CC);
	data[CG] += other.at(CG);
	data[CT] += other.at(CT);
	data[GG] += other.at(GG);
	data[GT] += other.at(GT);
	data[TT] += other.at(TT);
};

double TGenotypeData::sum(){
	return data[AA] + data[AC] + data[AG] + data[AT] + data[CC] + data[CG] + data[CT] + data[GG] + data[GT] + data[TT];
};

double TGenotypeData::weightedSum(const TGenotypeData & weights){
	return data[AA] * weights.at(AA)
			+ data[AC] * weights.at(AC)
		   	+ data[AG] * weights.at(AG)
			+ data[AT] * weights.at(AT)
			+ data[CC] * weights.at(CC)
			+ data[CG] * weights.at(CG)
			+ data[CT] * weights.at(CT)
			+ data[GG] * weights.at(GG)
			+ data[GT] * weights.at(GT)
			+ data[TT] * weights.at(TT);
};

void TGenotypeData::normalize(){
	double theSum = sum();

	data[AA] = data[AA] / theSum;
	data[AC] = data[AC] / theSum;
	data[AG] = data[AG] / theSum;
	data[AT] = data[AT] / theSum;
	data[CC] = data[CC] / theSum;
	data[CG] = data[CG] / theSum;
	data[CT] = data[CT] / theSum;
	data[GG] = data[GG] / theSum;
	data[GT] = data[GT] / theSum;
	data[TT] = data[TT] / theSum;
};

void TGenotypeData::write(TOutputFileZipped & out) const{
	out << data[AA];
	out << data[AC];
	out << data[AG];
	out << data[AT];
	out << data[CC];
	out << data[CG];
	out << data[CT];
	out << data[GG];
	out << data[GT];
	out << data[TT];
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

void TGenotypeLikelihoods::fill(const std::vector<TBaseLikelihoods> & bases){
	fill(bases, bases.size());
};

void TGenotypeLikelihoods::fill(const std::vector<TBaseLikelihoods> & bases, const size_t size){
	//allows for vector to be longer than what is to be used
	//do in log if depth is high
	if(bases.size() > 50){
		//initialize
		set(0.0);

		//add to log genotype data
		for(size_t i=0; i<size; ++i){
			data[AA] += log(bases[i].at(A));
			data[AC] += log(0.5*bases[i].at(A) + 0.5*bases[i].at(C));
			data[AG] += log(0.5*bases[i].at(A) + 0.5*bases[i].at(G));
			data[AT] += log(0.5*bases[i].at(A) + 0.5*bases[i].at(T));
			data[CC] += log(bases[i].at(C));
			data[CG] += log(0.5*bases[i].at(C) + 0.5*bases[i].at(G));
			data[CT] += log(0.5*bases[i].at(C) + 0.5*bases[i].at(T));
			data[GG] += log(bases[i].at(G));
			data[GT] += log(0.5*bases[i].at(G) + 0.5*bases[i].at(T));
			data[TT] += log(bases[i].at(T));
		}

		//standardize and de-log
		double max = *std::max_element(&data[AA], &data[TT]);
		data[AA] = exp(data[AA] - max);
		data[AC] = exp(data[AC] - max);
		data[AG] = exp(data[AG] - max);
		data[AT] = exp(data[AT] - max);
		data[CC] = exp(data[CC] - max);
		data[CG] = exp(data[CG] - max);
		data[CT] = exp(data[CT] - max);
		data[GG] = exp(data[GG] - max);
		data[GT] = exp(data[GT] - max);
		data[TT] = exp(data[TT] - max);
	} else { //on natural scale
		//initialize
		set(1.0);

		for(size_t i=0; i<size; ++i){
			data[AA] *= bases[i].at(A);
			data[AC] *= 0.5*bases[i].at(A) + 0.5*bases[i].at(C);
			data[AG] *= 0.5*bases[i].at(A) + 0.5*bases[i].at(G);
			data[AT] *= 0.5*bases[i].at(A) + 0.5*bases[i].at(T);
			data[CC] *= bases[i].at(C);
			data[CG] *= 0.5*bases[i].at(C) + 0.5*bases[i].at(G);
			data[CT] *= 0.5*bases[i].at(C) + 0.5*bases[i].at(T);
			data[GG] *= bases[i].at(G);
			data[GT] *= 0.5*bases[i].at(G) + 0.5*bases[i].at(T);
			data[TT] *= bases[i].at(T);
		}
	}
};

//--------------------------------------------------------------------
// TGenotypePosteriorProbabilities
//--------------------------------------------------------------------
void TGenotypePosteriorProbabilities::reset(){
	set(0.1);
};

void TGenotypePosteriorProbabilities::fill(const TGenotypeData & likelihoods, const TGenotypeData & prior){
	//calculate normalized genotype probabilities according to Bayes rule
	data[AA] = likelihoods.at(AA) * prior.at(AA);
	data[AC] = likelihoods.at(AC) * prior.at(AC);
	data[AG] = likelihoods.at(AG) * prior.at(AG);
	data[AT] = likelihoods.at(AT) * prior.at(AT);
	data[CC] = likelihoods.at(CC) * prior.at(CC);
	data[CG] = likelihoods.at(CG) * prior.at(CG);
	data[CT] = likelihoods.at(CT) * prior.at(CT);
	data[GG] = likelihoods.at(GG) * prior.at(GG);
	data[GT] = likelihoods.at(GT) * prior.at(GT);
	data[TT] = likelihoods.at(TT) * prior.at(TT);

	normalize();
};




}; // end namespace


