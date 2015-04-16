/*
 * loci.h
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#ifndef LOCI_H_
#define LOCI_H_

#include <math.h>
#include "stringFunctions.h"
#include "TLog.h"
#include "TParameters.h"
#include "bamtools/api/BamReader.h"
#include "bamtools/api/SamSequenceDictionary.h"
#include <vector>
#include <armadillo>

enum Base {A, C, G, T};
enum Genotype {AA, AC, AG, AT, CC, CG, CT, GG, GT, TT};


class TBaseFrequencies{
public:
	double freq[4];

	TBaseFrequencies(){
		for(int i = 0; i < 4; ++i) freq[i]=0.0;
	};
	void add(Base B, double & weight){
		freq[B] += weight;
	}
	void normalize(){
		double sum = 0.0;
		for(int i = 0; i < 4; ++i) sum += freq[i];
		for(int i = 0; i < 4; ++i) freq[i] /= sum;
	};
	void print(){
		std::cout << "freq(A) = " << freq[0] << ", freq(C) = " << freq[1] << ", freq(G) = " << freq[2] << ", freq(T) = " << freq[3] << std::endl;
	};
	double& operator[](int pos){
		return freq[pos];
	}
};

struct Constants{
	//post mortem damage
	double pdmC;
	double pdmLambda;

	//EM parameters
	int numIterations;
	double maxEpsilon;
	int NewtonRalphsonNumIterations;
	double NewtonRalphsonMaxF;

	//genotype map for enum type
	Genotype** genotypeMap; //mapping base numbering to genotype numbering

	Constants(){
		pdmC = 0.0;
		pdmLambda = 0.0;
		numIterations = -1;
		maxEpsilon = 0.0;
		NewtonRalphsonNumIterations = -1;
		NewtonRalphsonMaxF = 0.0;

		//set up genotype map
		genotypeMap = new Genotype*[4];
		for(int i=0; i<4; ++i)
			genotypeMap[i] = new Genotype[4];

		int geno = 0;
		for(int i=0; i<4; ++i){
			for(int j=i; j<4; ++j){
				genotypeMap[i][j] = static_cast<Genotype>(geno);
				genotypeMap[j][i] = genotypeMap[i][j];
				++geno;
			}
		}
	};

	~Constants(){
		for(int i=0; i<4; ++i) delete[] genotypeMap[i];
		delete[] genotypeMap;
	};

	Genotype getGenotype(Base first, Base second){
		return genotypeMap[first][second];
	};
	Genotype getGenotype(int first, int second){
		return genotypeMap[first][second];
	};
};

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

	void update(double & ErrorRate, int & Pos5, int & Pos3, Constants & constants){
		errorRate = ErrorRate;
		pos5 = Pos5;
		pos3 = Pos3;
		fillEmissionProbabilities(constants);
	};

	virtual void fillEmissionProbabilities(Constants & constants){
		throw "Not implemented for base class!";
	};

	double probPDM(int & pos, Constants & constants);

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
	TBaseA(double & ErrorRate, int & Pos5, int & Pos3, Constants & constants):TBase(ErrorRate, Pos5, Pos3){
		fillEmissionProbabilities(constants);
	};
	char getBase(){ return 'A'; };
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(A, weight); };
	void fillEmissionProbabilities(Constants & constants);
};

//---------------------------------------------------------------
class TBaseC:public TBase{
public:
	TBaseC(double & ErrorRate, int & Pos5, int & Pos3, Constants & constants):TBase(ErrorRate, Pos5, Pos3){
		fillEmissionProbabilities(constants);
	};
	char getBase(){ return 'C'; };
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(C, weight); };
	void fillEmissionProbabilities(Constants & constants);
};

//---------------------------------------------------------------
class TBaseG:public TBase{
public:
	TBaseG(double & ErrorRate, int & Pos5, int Pos3, Constants & constants):TBase(ErrorRate, Pos5, Pos3){
		fillEmissionProbabilities(constants);
	};
	char getBase(){ return 'G'; };
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(G, weight); };
	void fillEmissionProbabilities(Constants & constants);
};

//---------------------------------------------------------------
class TBaseT:public TBase{
public:
	TBaseT(double & ErrorRate, int & Pos5, int & Pos3, Constants & constants):TBase(ErrorRate, Pos5, Pos3){
		fillEmissionProbabilities(constants);
	};
	char getBase(){ return 'T'; };
	void addToBaseFrequencies(TBaseFrequencies & frequencies, double & weight){ frequencies.add(T, weight); };
	void fillEmissionProbabilities(Constants & constants);
};

//---------------------------------------------------------------
//TSite
//---------------------------------------------------------------
class TSite{
public:
	std::vector<TBase*> bases;
	double emissionProbabilities[10];
	double P_g[10]; //P(g|d, theta, pi), see equation (3)

	TSite(){};
	~TSite(){
		clear();
	};

	void clear(){
		for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it)
			delete *it;
		bases.clear();
	};

	void add(char & base, char & quality, int pos5, int pos3, Constants & constants);
	void addToBaseFrequencies(TBaseFrequencies & frequencies);
	void calcEmissionProbabilities();
	void calculateP_g(double* genotypeProbabilities);
	std::string getBases();
	std::string getEmissionProbs();

};

//---------------------------------------------------------------
//TWindow
//---------------------------------------------------------------
class TWindow{
public:
	long start;
	long end; //end NOT included in window!
	long length;
	TSite* sites;
	bool sitesInitialized;
	double coverage, fractionSitesNoData, fractionsitesCoverageAtLeastTwo;;
	TBaseFrequencies baseFreq;
	double theta;
	Constants* sharedConstants;
	TWindow(Constants* SharedConstants);
	TWindow(long Start, long End, Constants* SharedConstants);
	~TWindow(){
		if(sitesInitialized) delete[] sites;
	};
	void initSites(long newLength);
	void clear();
	void move(long Start, long End);
	bool addFromRead(BamTools::BamAlignment & bamAlignement);
	void runEM();
	void writeEMResults(std::ofstream & out);
	void printPileup();
	void calcCoverage();
};

//---------------------------------------------------------------
//TGenome
//---------------------------------------------------------------
class TGenome{
public:
	//std::vector<TChromosome> chromosomes;
	std::string filename;
	TLog* logfile;
	BamTools::BamReader bamReader;
 	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	BamTools::BamAlignment bamAlignement;
 	BamTools::SamSequenceIterator chrIterator;
 	int chrNumber;
 	long chrLength;
 	long curStart, curEnd;
 	long windowSize;
 	double maxMissing;
 	Constants sharedConstants;
 	long oldPos;
 	TWindow* curWindow;
 	TWindow* nextWindow;
 	std::string outputName;

	TGenome(TLog* Logfile, TParameters & params);
	~TGenome(){
		delete curWindow;
		delete nextWindow;
	};
	bool iterateChromosome();
	bool iterateWindow();
	bool readData();
	void estimateTheta();
};




#endif /* LOCI_H_ */
