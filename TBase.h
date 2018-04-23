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
//	int quality;
	int posInRead; //zero based!
	int posInReadRev; //zero based!
	double PMD_CT, PMD_GA;
	int readGroup;
	BaseContext context;
	TEmissionProbabilities emissionProbabilities;

	TBase(){
		base = N;
		errorRate = -1.0;
//		quality = -1;
		posInRead = -1;
		posInReadRev = -1;
		PMD_CT = -1.0;
		PMD_GA = -1.0;
		readGroup = -1;
		context = cNN;
	}

	TBase(Base & Base, double & ErrorRate, int & PosInRead, int & PosInReadRev, double & thisPMD_CT, double & thisPMD_GA,  BaseContext & Context, int & ReadGroup){
		base = Base;
		errorRate = ErrorRate;
//		quality = Quality;
		posInRead = PosInRead;
		posInReadRev = PosInReadRev;
		PMD_CT = thisPMD_CT;
		PMD_GA = thisPMD_GA;
		readGroup = ReadGroup;
		context = Context;
	};

	virtual ~TBase(){};

	void fillEmissionProbabilitiesCore(double thisErrorRate);

	double getEmissionProbability(Genotype genotype){
		return emissionProbabilities.get(genotype);
	};
	double getEmissionProbability(int genotype){
		return emissionProbabilities.get(genotype);
	};

	void printEmissionProbs(){
		emissionProbabilities.print();
	}
//	virtual char getBaseAsChar(){ return '?'; };
	Base getBaseAsEnum(){ return N;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(base, weight); };

};
/*
//---------------------------------------------------------------
class TBaseDiploidA:public TBase{
public:
	TBaseDiploidA(int & Quality, int & PosInRead, int & PosInReadRev, double & thisPMD_CT, double & thisPMD_GA,  BaseContext & Context, int & ReadGroup):TBase(Quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup){};
	char getBase(){ return 'A'; };
	Base getBaseAsEnum(){ return A;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(A, weight); };
	void fillEmissionProbabilitiesCore(double thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidC:public TBase{
public:
	TBaseDiploidC(int & Quality, int & PosInRead, int & PosInReadRev, double & thisPMD_CT, double & thisPMD_GA, BaseContext & Context, int & ReadGroup):TBase(Quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup){};
	char getBase(){ return 'C'; };
	Base getBaseAsEnum(){ return C;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(C, weight); };
	void fillEmissionProbabilitiesCore(double thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidG:public TBase{
public:
	TBaseDiploidG(int & Quality, int & PosInRead, int & PosInReadRev, double & thisPMD_CT, double & thisPMD_GA, BaseContext & Context, int & ReadGroup):TBase(Quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup){};
	char getBase(){ return 'G'; };
	Base getBaseAsEnum(){ return G;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(G, weight); };
	void fillEmissionProbabilitiesCore(double thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidT:public TBase{
public:
	TBaseDiploidT(int & Quality, int & PosInRead, int & PosInReadRev, double & thisPMD_CT, double & thisPMD_GA,  BaseContext & Context, int & ReadGroup):TBase(Quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup){};
	char getBase(){ return 'T'; };
	Base getBaseAsEnum(){ return T;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(T, weight); };
	void fillEmissionProbabilitiesCore(double thisErrorRate);
};

*/
#endif /* TBASE_H_ */
