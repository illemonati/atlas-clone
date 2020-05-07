/*
 * TBase.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TBASE_H_
#define TBASE_H_

#include "TGenotypeMap.h"
#include <bitset>

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
//data container with minimal footprint, also used in recal
class TBaseData{
private:
	//flags: isReverseStrand, isSecondMate, isAligned
	std::bitset<3> flags; //initialized as 0,0,0

public:
	Base base;
	uint8_t qualityAsPhredInt; //Note: original quality as in BAM file, but transformed to phredInt
	uint16_t distFrom5Prime; //zero based!
	uint16_t distFrom3Prime; //zero based!	Do we need it if we also store fragment length?
	uint16_t readGroup;
	uint16_t fragmentLength;
	uint8_t mappingQuality;
	BaseContext context;


	TBaseData(){
		base = N;
		qualityAsPhredInt = 0;
		distFrom5Prime = -1;
		distFrom3Prime = -1;
		readGroup = -1;
		context = cNN;
		mappingQuality = 0;
		fragmentLength = 0;
	};

	//set and get flags
	bool isReverseStrand(){ return flags[0]; };
	bool isSecondMate(){ return flags[1]; };
	bool isAligned(){ return flags[2]; };

	void setReverseStrand(const bool status){ flags[0] = status; };
	void setSecondMate(const bool status){ flags[1] = status; };
	void setAligned(const bool status){ flags[2] = status; };

};

class TBase{
public:
	TBaseData data;

	//tmp variables. Remove to fuse TBase with TBaseData?
	uint16_t alignedPos; //takes value -1 when base is not aligned
	double errorRate;
	double PMD_CT, PMD_GA;

	TBase(){
		errorRate = -1.0;
		PMD_CT = -1.0;
		PMD_GA = -1.0;
		alignedPos = 0;
	};
	~TBase(){};

	bool isReverseStrand(){ return data.isReverseStrand(); };
	bool isSecondMate(){ return data.isSecondMate(); };
	bool isAligned(){ return data.isAligned(); };

	void setReverseStrand(const bool status){ data.setReverseStrand(status); };
	void setSecondMate(const bool status){ data.setSecondMate(status); };
	void setAligned(const bool status){ data.setAligned(status); };

	void addToEmissionProb(double genotypeLikelihoods[10]);
	void addToEmissionProbLog(double genotypeLikelihoods[10]);
	Base getBaseAsEnum(){ return data.base;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(data.base, weight); };

};

#endif /* TBASE_H_ */
