/*
 * TBase.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TBASE_H_
#define TBASE_H_

#include "TGenotypeMap.h"

//---------------------------------------------------------------
//TEmissionProbabilities
//---------------------------------------------------------------
class TEmissionProbabilities{
public:
	double emission[10];

	TEmissionProbabilities(){
		for(int i=0; i<10; ++i) emission[i]=0.0;
	};
	void set(Genotype geno, double val){ emission[geno] = val; }
	double& get(Genotype geno){ return emission[geno]; };
	double& get(int geno){ return emission[geno]; };
	void print(){
		std::cout << "Emissions:";
		for(int i=0; i<10; ++i){
			std::cout << " " << emission[i];
		}
		std::cout << std::endl;
	};
};

//---------------------------------------------------------------
//TBase
//---------------------------------------------------------------

class TBase{
public:
	Base base;
	double errorRate;
	int distFrom5Prime; //zero based!
	int distFrom3Prime; //zero based!
	double PMD_CT, PMD_GA;
	uint16_t readGroup;
	BaseContext context;
	bool aligned;  //whether or not base is aligned to ref. Insertions and clipped bases are not aligned
	uint16_t alignedPos;
	bool isSecondMate; //false for single-end data as well as the first read of paired-end data. true for the second mate of paired-end data.
	bool isReverseStrand; //TODO: group three bool into one uint8_t to save two bytes.
	uint16_t mappingQuality;
	uint16_t fragmentLength;

	TBase(){
		base = N;
		errorRate = -1.0;
		distFrom5Prime = -1;
		distFrom3Prime = -1;
		PMD_CT = -1.0;
		PMD_GA = -1.0;
		readGroup = -1;
		context = cNN;
		aligned = false;
		alignedPos = 0;
		isSecondMate = false;
		isReverseStrand = false;
		mappingQuality = 0;
		fragmentLength = 0;
	};

	TBase(Base & Base, double & ErrorRate, int & PosInRead, int & DistFrom5Prime, int & DistFrom3Prime, double & thisPMD_CT, double & thisPMD_GA,  BaseContext & Context, int & ReadGroup, bool & Aligned, uint16_t & AlignedPos, bool & IsSecond, bool & IsReverseStrand, uint8_t & MappingQuality, uint16_t & FragmentLength){
		base = Base;
		errorRate = ErrorRate;
		distFrom5Prime = DistFrom5Prime;
		distFrom3Prime = DistFrom3Prime;
		PMD_CT = thisPMD_CT;
		PMD_GA = thisPMD_GA;
		readGroup = ReadGroup;
		context = Context;
		aligned = Aligned;
		alignedPos = AlignedPos;
		isSecondMate = IsSecond;
		isReverseStrand = IsReverseStrand;
		mappingQuality = MappingQuality;
		fragmentLength = FragmentLength;
	};

	~TBase(){};
	void addToEmissionProb(double genotypeLikelihoods[10]);
	void addToEmissionProbLog(double genotypeLikelihoods[10]);
	Base getBaseAsEnum(){ return base;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(base, weight); };

};

#endif /* TBASE_H_ */
