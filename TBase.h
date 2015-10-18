/*
 * TBase.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TBASE_H_
#define TBASE_H_

#include "TGenotypeMap.h"
#include "TPostMortemDamage.h"

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

class TEmissionProbabilitiesHaploid{
public:
	double emission[4];

	TEmissionProbabilitiesHaploid(){
		for(int i=0; i<4; ++i) emission[i]=0.0;
	};
	void set(Base geno, double val){ emission[geno] = val; }
	double& get(Base geno){ return emission[geno]; };
	double& get(int geno){ return emission[geno]; };
};

//---------------------------------------------------------------
//TBase
//---------------------------------------------------------------
class TBase{
public:
	int quality;
	//double errorRate;
	//double transformedLogError;
	int posInRead; //zero based!
	int pos5, pos3; //is distance and starts at 1 for position = 0
	int readGroup;
	BaseContext context;

	TBase(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup){
		quality = (int) Quality - 33;
		//errorRate = pow(10.0, (double) quality / -10.0);
		//transformedLogError = -log(1.0 / errorRate - 1.0);
		posInRead = PosInRead;
		pos5 = Pos5;
		pos3 = Pos3;
		readGroup = ReadGroup;
		context = Context;
	};

	virtual ~TBase(){};

	//void fillEmissionProbabilities(TPMD & pmdObject);
	virtual void fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate){
		throw "Function 'fillEmissionProbabilitiesCore' Not implemented for base class TBase!";
	};
	virtual char getBase(){ return '?'; };
	virtual Base getBaseAsEnum(){ return N;};
	virtual void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){};
	virtual double getEmissionProbability(int genotype){
		throw "Function 'getEmissionProbability' Not implemented for base class TBase!";
	};
};

class TBaseDiploid:public TBase{
public:
	TEmissionProbabilitiesDiploid emissionProbabilities;

	TBaseDiploid(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBase(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};

	virtual ~TBaseDiploid(){};

	double getEmissionProbability(Genotype genotype){
		return emissionProbabilities.get(genotype);
	};
	double getEmissionProbability(int genotype){
		return emissionProbabilities.get(genotype);
	};

};

class TBaseHaploid:public TBase{
public:
	TEmissionProbabilitiesHaploid emissionProbabilities;

	TBaseHaploid(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBase(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	virtual ~TBaseHaploid(){};

	double getEmissionProbability(Base genotype){
		return emissionProbabilities.get(genotype);
	};
	double getEmissionProbability(int genotype){
		return emissionProbabilities.get(genotype);
	};
};
//---------------------------------------------------------------
class TBaseDiploidA:public TBaseDiploid{
public:
	TBaseDiploidA(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseDiploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'A'; };
	Base getBaseAsEnum(){ return A;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(A, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate);
};
class TBaseHaploidA:public TBaseHaploid{
public:
	TBaseHaploidA(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseHaploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'A'; };
	Base getBaseAsEnum(){ return A;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(A, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidC:public TBaseDiploid{
public:
	TBaseDiploidC(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseDiploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'C'; };
	Base getBaseAsEnum(){ return C;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(C, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate);
};
class TBaseHaploidC:public TBaseHaploid{
public:
	TBaseHaploidC(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseHaploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'C'; };
	Base getBaseAsEnum(){ return C;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(C, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidG:public TBaseDiploid{
public:
	TBaseDiploidG(char & Quality, int & PosInRead, int & Pos5, int Pos3, BaseContext & Context, int & ReadGroup):TBaseDiploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'G'; };
	Base getBaseAsEnum(){ return G;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(G, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate);
};
class TBaseHaploidG:public TBaseHaploid{
public:
	TBaseHaploidG(char & Quality, int & PosInRead, int & Pos5, int Pos3, BaseContext & Context, int & ReadGroup):TBaseHaploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'G'; };
	Base getBaseAsEnum(){ return G;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(G, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidT:public TBaseDiploid{
public:
	TBaseDiploidT(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseDiploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'T'; };
	Base getBaseAsEnum(){ return T;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(T, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate);
};
class TBaseHaploidT:public TBaseHaploid{
public:
	TBaseHaploidT(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseHaploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'T'; };
	Base getBaseAsEnum(){ return T;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(T, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, double thisErrorRate);
};


#endif /* TBASE_H_ */
