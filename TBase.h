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
class TEmissionProbabilitiesDiploid{
public:
	double emission[10];

	TEmissionProbabilitiesDiploid(){
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
/*class TBase{
public:
	int quality;
	int posInRead; //zero based!
	int posInReadRev; //zero based!
	double PMD_CT, PMD_GA;
	int readGroup;
	BaseContext context;

	TBase(int & Quality, int & PosInRead, int & PosInReadRev, double & thisPMD_CT, double & thisPMD_GA, BaseContext & Context, int & ReadGroup){
		quality = Quality;
		//errorRate = pow(10.0, (double) quality / -10.0);
		//transformedLogError = -log(1.0 / errorRate - 1.0);
		posInRead = PosInRead;
		posInReadRev = PosInReadRev;
		PMD_CT = thisPMD_CT;
		PMD_GA = thisPMD_GA;
		readGroup = ReadGroup;
		context = Context;
	};

	virtual ~TBase(){};

	//void fillEmissionProbabilities(TPMD & pmdObject);
	virtual void fillEmissionProbabilitiesCore(double thisErrorRate){
		throw "Function 'fillEmissionProbabilitiesCore' Not implemented for base class TBase!";
	};

	virtual void printEmissionProbs(){
		throw "Function 'printEmissionProbs' Not implemented for base class TBase!";
	};

	virtual char getBase(){ return '?'; };
	virtual Base getBaseAsEnum(){ return N;};
	virtual void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){};
	virtual double getEmissionProbability(int genotype){
		throw "Function 'getEmissionProbability' Not implemented for base class TBase!";
	};
};*/

class TBase{
public:
	int quality;
	int posInRead; //zero based!
	int posInReadRev; //zero based!
	double PMD_CT, PMD_GA;
	int readGroup;
	BaseContext context;
	TEmissionProbabilitiesDiploid emissionProbabilities;

	TBase(int & Quality, int & PosInRead, int & PosInReadRev, double & thisPMD_CT, double & thisPMD_GA,  BaseContext & Context, int & ReadGroup){
		quality = Quality;
		posInRead = PosInRead;
		posInReadRev = PosInReadRev;
		PMD_CT = thisPMD_CT;
		PMD_GA = thisPMD_GA;
		readGroup = ReadGroup;
		context = Context;
	};

	virtual ~TBase(){};

	virtual void fillEmissionProbabilitiesCore(double thisErrorRate){
		throw "Function 'fillEmissionProbabilitiesCore' Not implemented for base class TBase!";
	};

	double getEmissionProbability(Genotype genotype){
		return emissionProbabilities.get(genotype);
	};
	double getEmissionProbability(int genotype){
		return emissionProbabilities.get(genotype);
	};

	void printEmissionProbs(){
		emissionProbabilities.print();
	}
	virtual char getBase(){ return '?'; };
	virtual Base getBaseAsEnum(){ return N;};
	virtual void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){};

};

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


#endif /* TBASE_H_ */
