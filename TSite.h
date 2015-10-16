/*
 * TBase.h
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#ifndef TSITE_H_
#define TSITE_H_

#include "stringFunctions.h"
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>
#include "TParameters.h"

enum Base {A=0, C, G, T, N};
enum Genotype {AA=0, AC, AG, AT, CC, CG, CT, GG, GT, TT};
enum PMDType {pmdCT=0, pmdGA};
enum BaseContext {cAA=0, cAC, cAG, cAT, cAN, cCA, cCC, cCG, cCT, cCN, cGA, cGC, cGG, cGT, cGN, cTA, cTC, cTG, cTT, cTN, cNA, cNC, cNG, cNT, cNN}; //N means "nothing", i.e. end of read or del

//---------------------------------------------------------------
//GenotypeMap
//---------------------------------------------------------------
//genotype map for enum type
struct GenotypeMap{
	Genotype** genotypeMap; //mapping base numbering to genotype enum
	BaseContext** contextMap; //mapping dinucleotide context to context enum

	GenotypeMap();
	~GenotypeMap(){
		for(int i=0; i<4; ++i){
			delete[] genotypeMap[i];
		}
		for(int i=0; i<5; ++i){
			delete[] contextMap[i];
		}
		delete[] genotypeMap;
		delete[] contextMap;
	};

	Base getBase(char & base){
		if(base == 'A') return A;
		if(base == 'C') return C;
		if(base == 'G') return G;
		if(base == 'T') return T;
		return N;
	};
	Genotype getGenotype(Base first, Base second){
		return genotypeMap[first][second];
	};
	Genotype getGenotype(int first, int second){
		return genotypeMap[first][second];
	};
	std::string getGenotypeString(int num);
	BaseContext getContext(Base first, Base second){
		return contextMap[first][second];
	};
	BaseContext getContext(int first, int second){
		return contextMap[first][second];
	};
	BaseContext getContext(char first, char second){
		return contextMap[getBase(first)][getBase(second)];
	};
};

//---------------------------------------------------------------
//TBaseFrequencies
//---------------------------------------------------------------
class TBaseFrequencies{
public:
	double freq[4];

	TBaseFrequencies(){
		for(int i = 0; i < 4; ++i) freq[i] = 0.0;
	};
	void add(Base B, double & weight){
		freq[B] += weight;
	}
	void normalize(){
		double sum = 0.0;
		for(int i = 0; i < 4; ++i) sum += freq[i];
		for(int i = 0; i < 4; ++i) freq[i] /= sum;
	};
	void clear(){
		for(int i = 0; i < 4; ++i) freq[i] = 0.0;
	};
	void print(){
		std::cout << "freq(A) = " << freq[0] << ", freq(C) = " << freq[1] << ", freq(G) = " << freq[2] << ", freq(T) = " << freq[3] << std::endl;
	};
	double& operator[](int pos){
		return freq[pos];
	}
};


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
//TPMD
//---------------------------------------------------------------
//Note: Base class is to be used when there is no PMD!
class TPMDFunction{
public:
	TPMDFunction(){};
	virtual ~TPMDFunction(){};
	virtual double getProb(double & pos){
		return 0.0;
	};
	virtual std::string getString(){ return "P(pmd|pos) = 0.0"; };
};

class TPMDSkoglund:public TPMDFunction{
private:
	double lambda, c;

public:
	TPMDSkoglund(double & Lambda, double & C){
		lambda = Lambda; c = C;
	};
	~TPMDSkoglund(){};
	double getProb(double & pos){
		return c + pow(1.0 - lambda, (double) pos - 1.0) * lambda;
	};
	std::string getString(){ return "P(pmd|pos) = p * (1 - p)^pos + c = " + toString(lambda) + " * (1 - " + toString(lambda) + ")^pos + " + toString(c); };
};

class TPMDVeeramah:public TPMDFunction{
private:
	double a,b,c;

public:
	TPMDVeeramah(double & A, double & B, double & C){
		a = A; b = B; c = C;
	}
	double getProb(double & pos){
		return a * exp(-pos * b) + c;
	};
	std::string getString(){ return "P(pmd|pos) = a * exp(- pos * b) + c = " + toString(a) + "* exp(- pos * " + toString(b) + ") + " + toString(c); };
};

class TPMD{
private:
	TPMDFunction* myFunctions[2];
	bool functionsInitialized[2];

public:
	TPMD();
	~TPMD(){
		if(functionsInitialized[pmdCT]) delete myFunctions[pmdCT];
		if(functionsInitialized[pmdGA]) delete myFunctions[pmdGA]; };
	void initializeFunction(std::string & pmdString, PMDType type);
	double getProb(double pos, PMDType type){ return myFunctions[type]->getProb(pos); };
	double getProbCT(double pos){ return myFunctions[pmdCT]->getProb(pos); };
	double getProbGA(double pos){ return myFunctions[pmdGA]->getProb(pos); };
	std::string getFunctionString(PMDType type){ return myFunctions[type]->getString(); };
};



//---------------------------------------------------------------
//TBase
//---------------------------------------------------------------
class TBase{
public:
	int quality;
	double errorRate;
	double transformedLogError;
	int posInRead; //zero based!
	int pos5, pos3; //is distance and starts at 1 for position = 0
	int readGroup;
	BaseContext context;

	TBase(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup){
		quality = (int) Quality - 33;
		errorRate = pow(10.0, (double) quality / -10.0);
		transformedLogError = -log(1.0 / errorRate - 1.0);
		posInRead = PosInRead;
		pos5 = Pos5;
		pos3 = Pos3;
		readGroup = ReadGroup;
		context = Context;
	};

	virtual ~TBase(){};

	void update(char & Quality, int & PosInRead, int & Pos5, int & Pos3, TPMD & pmdObject){
		quality = (int) Quality - 33;
		errorRate = pow(10.0, (double) quality / -10.0);
		posInRead = PosInRead;
		pos5 = Pos5;
		pos3 = Pos3;
		fillEmissionProbabilities(pmdObject);
	};

	void fillEmissionProbabilities(TPMD & pmdObject);
	virtual void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate){
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
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
class TBaseHaploidA:public TBaseHaploid{
public:
	TBaseHaploidA(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseHaploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'A'; };
	Base getBaseAsEnum(){ return A;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(A, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidC:public TBaseDiploid{
public:
	TBaseDiploidC(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseDiploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'C'; };
	Base getBaseAsEnum(){ return C;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(C, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
class TBaseHaploidC:public TBaseHaploid{
public:
	TBaseHaploidC(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseHaploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'C'; };
	Base getBaseAsEnum(){ return C;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(C, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidG:public TBaseDiploid{
public:
	TBaseDiploidG(char & Quality, int & PosInRead, int & Pos5, int Pos3, BaseContext & Context, int & ReadGroup):TBaseDiploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'G'; };
	Base getBaseAsEnum(){ return G;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(G, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
class TBaseHaploidG:public TBaseHaploid{
public:
	TBaseHaploidG(char & Quality, int & PosInRead, int & Pos5, int Pos3, BaseContext & Context, int & ReadGroup):TBaseHaploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'G'; };
	Base getBaseAsEnum(){ return G;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(G, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidT:public TBaseDiploid{
public:
	TBaseDiploidT(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseDiploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'T'; };
	Base getBaseAsEnum(){ return T;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(T, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
class TBaseHaploidT:public TBaseHaploid{
public:
	TBaseHaploidT(char & Quality, int & PosInRead, int & Pos5, int & Pos3, BaseContext & Context, int & ReadGroup):TBaseHaploid(Quality, PosInRead, Pos5, Pos3, Context, ReadGroup){};
	char getBase(){ return 'T'; };
	Base getBaseAsEnum(){ return T;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(T, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
//---------------------------------------------------------------
//TSite
//---------------------------------------------------------------
class TSite{
public:
	bool hasData;
	std::vector<TBase*> bases;
	int numGenotypes;
	double* emissionProbabilities;
	double* P_g; //P(g|d, theta, pi), see equation (3)
	char referenceBase; //optional

	TSite(){
		hasData = false;
		numGenotypes = 0;
		emissionProbabilities = NULL;
		P_g = NULL;
		referenceBase = 'N';
	};
	virtual ~TSite(){ clear(); };

	void clear();

	virtual void add(char & base, char & quality, int PosInRead, int pos5, int pos3, BaseContext & Context, int & ReadGroup){throw "Function 'add' Not implemented for base class TSite!"; };
	void setRefBase(char & Base){referenceBase = Base; };
	char getRefBase(){return referenceBase;};
	Base getRefBaseAsEnum(){
		if(referenceBase == 'A') return A;
		if(referenceBase == 'C') return C;
		if(referenceBase == 'G') return G;
		if(referenceBase == 'T') return T;
		return N;
	};
	void addToBaseFrequencies(TBaseFrequencies & frequencies);
	void calcEmissionProbabilities(TPMD & pmdObject);
	void callMLEGenotype(TPMD & pmdObject, GenotypeMap & genoMap, std::ofstream & out);
	void calculateP_g(double* genotypeProbabilities);
	double calculateWeightedSumOfEmissionProbs(double* weights);
	std::string getBases();
	std::string getEmissionProbs();
	double calculateLogLikelihood(double* genotypeProbabilities);
};

class TSiteDiploid:public TSite{
public:

	TSiteDiploid(){
		hasData = false;
		numGenotypes = 10;
		emissionProbabilities = new double[numGenotypes];
		P_g = new double[numGenotypes];
	};
	~TSiteDiploid(){
		delete[] emissionProbabilities;
		delete[] P_g;
	};
	void add(char & base, char & quality, int PosInRead, int pos5, int pos3, BaseContext & Context, int & ReadGroup);
};

class TSiteHaploid:public TSite{
public:

	TSiteHaploid(){
		hasData = false;
		numGenotypes = 4;
		emissionProbabilities = new double[numGenotypes];
		P_g = new double[numGenotypes];
	}
	~TSiteHaploid(){
		delete[] emissionProbabilities;
		delete[] P_g;
	};
	void add(char & base, char & quality, int PosInRead, int pos5, int pos3, BaseContext & Context, int & ReadGroup);
};


#endif /* TSITE_H_ */
