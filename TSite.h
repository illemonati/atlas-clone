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

//---------------------------------------------------------------
//GenotypeMap
//---------------------------------------------------------------
//genotype map for enum type
struct GenotypeMap{
	Genotype** genotypeMap; //mapping base numbering to genotype numbering

	GenotypeMap();
	~GenotypeMap(){
		for(int i=0; i<4; ++i)
			delete[] genotypeMap[i];
		delete[] genotypeMap;
	};
	Genotype getGenotype(Base first, Base second){
		return genotypeMap[first][second];
	};
	Genotype getGenotype(int first, int second){
		return genotypeMap[first][second];
	};
	std::string getGenotypeString(int num);
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
//Recalibration
//---------------------------------------------------------------
class TRecalibration{
public:
	bool doRecalibration;
	double a,b,c;

	TRecalibration(){
		doRecalibration = false;
		a = 0.0;
		b = 0.0;
		c = 1.0;
	};
	TRecalibration(const double & paramA, const double & paramB, const double & paramC){
		doRecalibration = true;
		a = paramA;
		b = paramB;
		c = paramC;
	};
	TRecalibration(std::string recalString);
	~TRecalibration(){};
	void set(const double & paramA, const double & paramB, const double & paramC){
		a = paramA;
		b = paramB;
		c = paramC;
	};
	double recalibrate(double & error);
	std::string getFunctionString();
};

//---------------------------------------------------------------
//RecalibrationEM
//---------------------------------------------------------------
class TBase;
class TRecalibrationEM{
public:
	int numParams;
	double* params;
	double* newParams; //used during EM
	arma::mat Jacobian;
	arma::vec F;
	arma::mat JxF;
	double maxF; //largest change during Newton-Ralphson
	long numSitesAdded;
	double logLikelihood;

	TRecalibrationEM(TParameters* arguments, TLog* logfile);
	~TRecalibrationEM(){
		delete[] params;
		delete[] newParams;
	};
	TRecalibrationEM(TLog* logfile);
	void setParams(double* Params);
	double calcEta(TBase* base);
	double calcEta(TBase* base, double* theseParams);
	double calcEpsilon(const double & eta);
	double calcEpsilon(TBase* base);
	double calcEpsilon(TBase* base, double* theseParams);
	void initEMStep();
	void initNetwonRalphsonStep();
	void saveParams();
	void addSiteToJacobianAndF(std::vector<TBase*> & bases, TBaseFrequencies* freqs);
	void runNewtonRalphson();
	void writeHeader(std::ofstream & out);
	void writeParams(std::ofstream & out);
	void resetLikelihood();
	void addSiteToLikelihood(std::vector<TBase*> & bases, TBaseFrequencies* freqs);
};

//---------------------------------------------------------------
//RecalibrationBQSR
//---------------------------------------------------------------
//covariates to take into account:
// - read base (A, G, C, T)
// -
class TRecalibrationBQSR{
public:
	//Table -> list of all cells


	TRecalibrationBQSR(TParameters* arguments, TLog* logfile);
	~TRecalibrationBQSR(){};

	//void addSites(TSite* site);

};

class TrecalibrationBQSR_cell{
public:
	double curEpsilon;
	bool estimationConverged;

	TrecalibrationBQSR_cell(double Error);
	virtual ~TrecalibrationBQSR_cell();

	virtual void addBase(TBase* base, char & RefBase){ throw "addBase not defined for base class 'TrecalibrationBQSR_cell'!"; };
	virtual bool estimateEpsilon();
};


class TrecalibrationBQSR_cellC:public TrecalibrationBQSR_cell{
public:
	long N_1, N_2;
	double D;
	TrecalibrationBQSR_cellC(double Error, int pos, TPMD* PmdObject);
	~TrecalibrationBQSR_cellC();

	void addBase(TBase* base, char & RefBase);
	virtual bool estimateEpsilon();
};

class TrecalibrationBQSR_cellT:public TrecalibrationBQSR_cellC{
public:
	long N_3;

	TrecalibrationBQSR_cellT(double Error, int pos, TPMD* PmdObject);
	~TrecalibrationBQSR_cellT();

	void addBase(TBase* base, char & RefBase);
	bool estimateEpsilon();
};

class TrecalibrationBQSR_cellA:public TrecalibrationBQSR_cell{
public:
	double firstDerivative, secondDerivative;
	TPMD* pmdObject;
	double convergenceThreshold;


	TrecalibrationBQSR_cellA(double Error, TPMD* PmdObject, double ConvergenceThreshold);
	~TrecalibrationBQSR_cellA();

	virtual double getD(TBase* base, char & RefBase);
	void addBase(TBase* base, char & RefBase);
	virtual bool estimateEpsilon();
};

class TrecalibrationBQSR_cellG:public TrecalibrationBQSR_cellA{
public:
	TrecalibrationBQSR_cellG(double Error, TPMD* PmdObject, double ConvergenceThreshold);
	~TrecalibrationBQSR_cellG();

	double getD(TBase* base, char & RefBase);
};



//---------------------------------------------------------------
//TBase
//---------------------------------------------------------------
class TBase{
public:
	double logError;
	double errorRate;
	double transformedLogError;
	int posInRead, pos5, pos3;
	int readGroup;

	TBase(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, int & ReadGroup){
		logError = LogErrorRate;
		errorRate = pow(10.0, logError);
		transformedLogError = -log(1.0 / errorRate - 1.0);
		posInRead = PosInRead;
		pos5 = Pos5;
		pos3 = Pos3;
		readGroup = ReadGroup;
	};

	virtual ~TBase(){};

	void update(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, TPMD & pmdObject){
		logError = LogErrorRate;
		errorRate = pow(10.0, logError);
		posInRead = PosInRead;
		pos5 = Pos5;
		pos3 = Pos3;
		fillEmissionProbabilities(pmdObject);
	};
	void fillEmissionProbabilities(TPMD & pmdObject);
	void fillEmissionProbabilitiesRecalibratedError(TPMD & pmdObject, TRecalibration & recal);
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

	TBaseDiploid(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, int & ReadGroup):TBase(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};

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

	TBaseHaploid(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, int & ReadGroup):TBase(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};
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
	TBaseDiploidA(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, int & ReadGroup):TBaseDiploid(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};
	char getBase(){ return 'A'; };
	Base getBaseAsEnum(){ return A;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(A, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
class TBaseHaploidA:public TBaseHaploid{
public:
	TBaseHaploidA(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, int & ReadGroup):TBaseHaploid(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};
	char getBase(){ return 'A'; };
	Base getBaseAsEnum(){ return A;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(A, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidC:public TBaseDiploid{
public:
	TBaseDiploidC(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, int & ReadGroup):TBaseDiploid(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};
	char getBase(){ return 'C'; };
	Base getBaseAsEnum(){ return C;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(C, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
class TBaseHaploidC:public TBaseHaploid{
public:
	TBaseHaploidC(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, int & ReadGroup):TBaseHaploid(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};
	char getBase(){ return 'C'; };
	Base getBaseAsEnum(){ return C;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(C, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidG:public TBaseDiploid{
public:
	TBaseDiploidG(double & LogErrorRate, int & PosInRead, int & Pos5, int Pos3, int & ReadGroup):TBaseDiploid(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};
	char getBase(){ return 'G'; };
	Base getBaseAsEnum(){ return G;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(G, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
class TBaseHaploidG:public TBaseHaploid{
public:
	TBaseHaploidG(double & LogErrorRate, int & PosInRead, int & Pos5, int Pos3, int & ReadGroup):TBaseHaploid(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};
	char getBase(){ return 'G'; };
	Base getBaseAsEnum(){ return G;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(G, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
//---------------------------------------------------------------
class TBaseDiploidT:public TBaseDiploid{
public:
	TBaseDiploidT(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, int & ReadGroup):TBaseDiploid(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};
	char getBase(){ return 'T'; };
	Base getBaseAsEnum(){ return T;};
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(T, weight); };
	void fillEmissionProbabilitiesCore(TPMD & pmdObject, const double & thisErrorRate);
};
class TBaseHaploidT:public TBaseHaploid{
public:
	TBaseHaploidT(double & LogErrorRate, int & PosInRead, int & Pos5, int & Pos3, int & ReadGroup):TBaseHaploid(LogErrorRate, PosInRead, Pos5, Pos3, ReadGroup){};
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

	double qualityToLogError(char & quality);
	virtual void add(char & base, char & quality, int PosInRead, int pos5, int pos3, int & ReadGroup){throw "Function 'add' Not implemented for base class TSite!"; };
	void setRefBase(char & Base){referenceBase = Base; };
	void addToBaseFrequencies(TBaseFrequencies & frequencies);
	void calcEmissionProbabilities(TPMD & pmdObject);
	void calcEmissionProbabilitiesScaledError(TPMD & pmdObject, TRecalibration & recal);
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
	void add(char & base, char & quality, int PosInRead, int pos5, int pos3, int & ReadGroup);
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
	void add(char & base, char & quality, int PosInRead, int pos5, int pos3, int & ReadGroup);
};


#endif /* TSITE_H_ */
