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

enum Base {A, C, G, T};
enum Genotype {AA, AC, AG, AT, CC, CG, CT, GG, GT, TT};

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
//TPMD
//---------------------------------------------------------------
//Note: Base class is to be used when there is no PMD!
class TPMDFunction{
public:
	TPMDFunction(){};
	virtual ~TPMDFunction(){};
	virtual double getProb(double pos){
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
	virtual std::string getString(){ return "P(pmd|pos) = p * (1 - p)^pos + c = " + toString(lambda) + " * (1 - " + toString(lambda) + ")^pos + " + toString(c); };
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
	virtual std::string getString(){ return "P(pmd|pos) = a * exp(- pos * b) + c = " + toString(a) + "* exp(- pos * " + toString(b) + ") + " + toString(c); };
};

class TPMD{
private:
	TPMDFunction* myFunction;
	bool functionInitialized;

public:
	TPMD(){ myFunction = NULL; functionInitialized = false; };
	~TPMD(){ if(functionInitialized) delete myFunction; };
	void initializeFunction(std::string & pmdString);
	double getProb(double pos){ return myFunction->getProb(pos); };
	std::string getFunctionString(){ return myFunction->getString(); };
};

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
};
//---------------------------------------------------------------
//TBase
//---------------------------------------------------------------
class TBase{
public:
	double errorRate;
	int pos5, pos3;
	TEmissionProbabilities emissionProbabilities;

	TBase(double & ErrorRate, int & Pos5, int & Pos3){
		errorRate = ErrorRate;
		pos5 = Pos5;
		pos3 = Pos3;
	};

	virtual ~TBase(){};

	void update(double & ErrorRate, int & Pos5, int & Pos3, TPMD* pmdCT, TPMD* pmdGA){
		errorRate = ErrorRate;
		pos5 = Pos5;
		pos3 = Pos3;
		fillEmissionProbabilities(pmdCT, pmdGA);
	};

	virtual void fillEmissionProbabilities(TPMD* pmdCT, TPMD* pmdGA){
		throw "Not implemented for base class!";
	};

	double getEmissionProbability(Genotype genotype){
		return emissionProbabilities.get(genotype);
	};
	double getEmissionProbability(int genotype){
		return emissionProbabilities.get(genotype);
	};

	virtual char getBase(){ return '?'; };
	virtual void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){};
};

//---------------------------------------------------------------
class TBaseA:public TBase{
public:
	TBaseA(double & ErrorRate, int & Pos5, int & Pos3, TPMD* pmdCT, TPMD* pmdGA):TBase(ErrorRate, Pos5, Pos3){
		fillEmissionProbabilities(pmdCT, pmdGA);
	};
	char getBase(){ return 'A'; };
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(A, weight); };
	void fillEmissionProbabilities(TPMD* pmdCT, TPMD* pmdGA);
};

//---------------------------------------------------------------
class TBaseC:public TBase{
public:
	TBaseC(double & ErrorRate, int & Pos5, int & Pos3, TPMD* pmdCT, TPMD* pmdGA):TBase(ErrorRate, Pos5, Pos3){
		fillEmissionProbabilities(pmdCT, pmdGA);
	};
	char getBase(){ return 'C'; };
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(C, weight); };
	void fillEmissionProbabilities(TPMD* pmdCT, TPMD* pmdGA);
};

//---------------------------------------------------------------
class TBaseG:public TBase{
public:
	TBaseG(double & ErrorRate, int & Pos5, int Pos3, TPMD* pmdCT, TPMD* pmdGA):TBase(ErrorRate, Pos5, Pos3){
		fillEmissionProbabilities(pmdCT, pmdGA);
	};
	char getBase(){ return 'G'; };
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(G, weight); };
	void fillEmissionProbabilities(TPMD* pmdCT, TPMD* pmdGA);
};

//---------------------------------------------------------------
class TBaseT:public TBase{
public:
	TBaseT(double & ErrorRate, int & Pos5, int & Pos3, TPMD* pmdCT, TPMD* pmdGA):TBase(ErrorRate, Pos5, Pos3){
		fillEmissionProbabilities(pmdCT, pmdGA);
	};
	char getBase(){ return 'T'; };
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(T, weight); };
	void fillEmissionProbabilities(TPMD* pmdCT, TPMD* pmdGA);
};


//---------------------------------------------------------------
//TSite
//---------------------------------------------------------------
class TSite{
public:
	std::vector<TBase*> bases;
	double emissionProbabilities[10];
	double P_g[10]; //P(g|d, theta, pi), see equation (3)
	bool hasData;

	TSite(){
		hasData = false;
	};
	~TSite(){
		clear();
	};

	void clear(){
		for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it)
			delete *it;
		bases.clear();
		hasData = false;
	};

	void add(char & base, char & quality, int pos5, int pos3, TPMD* pmdCT, TPMD* pmdGA);
	void addToBaseFrequencies(TBaseFrequencies & frequencies);
	void calcEmissionProbabilities();
	void calculateP_g(double* genotypeProbabilities);
	std::string getBases();
	std::string getEmissionProbs();
	double calculateLogLikelihood(double* genotypeProbabilities);
};

#endif /* TSITE_H_ */
