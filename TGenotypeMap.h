/*
 * TGenotypeMap.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TGENOTYPEMAP_H_
#define TGENOTYPEMAP_H_

#include "stringFunctions.h"
#include <math.h>

enum Base : uint8_t {A=0, C, G, T, N};
enum Genotype : uint8_t {AA=0, AC, AG, AT, CC, CG, CT, GG, GT, TT};
enum BaseContext : uint8_t {cAA=0, cAC, cAG, cAT, cCA, cCC, cCG, cCT, cGA, cGC, cGG, cGT, cTA, cTC, cTG, cTT, cNA, cNC, cNG, cNT, cAN, cCN, cGN, cTN, cNN}; //N means unknwon base or "nothing", i.e. end of read or del

//---------------------------------------------------------------
//TGenotypeMap
//---------------------------------------------------------------
//genotype map for enum type
class TGenotypeMap{
public:
	Genotype** genotypeMap; //mapping base numbering to genotype enum
	BaseContext** contextMap; //mapping dinucleotide context to context enum
	Base** genotypeToBase; //mapping genotypes to bases
	Base** alleleicCombinationToBase; //mapping the allelic combination to the two bases
	Genotype** alleleicCombinationToGenotypes; //mapping the alleleic combinations to the homozygous, heterozygosu and other homozygous genotype
	char** alleleicCombinationToBaseChar; //mapping the allelic combination to the two bases

	char* baseToChar;
	Base* baseToFlippedBase;
	int numGenotypes;
	int numContexts;
	int numContextsNotN;

    TGenotypeMap();
    ~TGenotypeMap();

	//TODO: also make an array to speed up?
    Base getBase(const char base);
    Base getBaseOnlyCapitals(const char base);
    char getBaseAsChar(Base base);
    char getBaseAsChar(int base);
    Base flipBase(char & base);
    Genotype getGenotype(Base first, Base second);
    Genotype getGenotype(int first, int second);
    Genotype getGenotype(char first, char second);
    std::string getGenotypeString(int num);
    std::string getGenotypeString(Genotype geno);
    std::string getGenotypeStringKnownAlleles(int num, char ref, char alt);
    std::pair<Base,Base> getBasesOfGenotype(int num);
    int getNumContext();
    BaseContext getContext(Base first, Base second);
    BaseContext getContext(int first, int second);
    BaseContext getContext(char first, char second);
    BaseContext getContextReverseRead(char first, char second);
    std::string getContextString(int num);
};

//---------------------------------------------------------------
//TBaseFrequencies
//---------------------------------------------------------------
class TBaseFrequencies{
public:
	double freq[4];
	bool wasNormalized;

    TBaseFrequencies();

    void add(Base B, double & weight);
    void addNoRef(Base B, double weight);
    void normalize();
    void setEqualBaseFreq();
    void clear();
    void print();
    double& operator[](int pos);
};

/*
//---------------------------------------------------------------
//TQualityMap
//---------------------------------------------------------------
class TQualityMap{
public:
	//IMPORTANT NOMENCLATURE
	//error is error rate between 0 and 1
	//phred is phred-scaled error as phred = -10 * log10(error)
	//phredInt is (int) phred
	//quality is phredInt + 33
	double* phredIntToErrorMap;
	double* qualityToErrorMap;
	double* phredIntToLogErrorMap;
	int* illuminaQualityBins;
	double minError;
	int maxPhred, sizeQual;
	double maxError;

	TQualityMap(){
		//only up to phred = 255, else always return 255
		maxPhred = 255;
		maxError = 0.9999999999;
		sizeQual = maxPhred + 33;
		phredIntToErrorMap = new double[maxPhred + 1];
		phredIntToLogErrorMap = new double[maxPhred + 1];
		qualityToErrorMap = new double[sizeQual];

		//initialize quality <= 0
		for(int i=0; i<34; ++i)
			qualityToErrorMap[i] = maxError;
		phredIntToErrorMap[0] = maxError;
		phredIntToLogErrorMap[0] = 0.0;

//		phredIntToErrorMap[33] = 0.9999999999;

		//and now others
		double tmp = - log(10) / 10.0;
		for(int i=1; i<maxPhred; ++i){
			phredIntToErrorMap[i] = phredToError(i);
			phredIntToLogErrorMap[i] = (double) i * tmp;
			qualityToErrorMap[i+33] = phredIntToErrorMap[i];
		}

		minError = phredToError(maxPhred);

		//Create map of illumina quality bins
		illuminaQualityBins = new int[sizeQual];
		for(int i=0; i<35; ++i)
			illuminaQualityBins[i] = 33;
		for(int i=35; i<43; ++i)
			illuminaQualityBins[i] = 39;
		for(int i=43; i<53; ++i)
			illuminaQualityBins[i] = 48;
		for(int i=53; i<58; ++i)
			illuminaQualityBins[i] = 55;
		for(int i=58; i<63; ++i)
			illuminaQualityBins[i] = 60;
		for(int i=63; i<68; ++i)
			illuminaQualityBins[i] = 66;
		for(int i=68; i<72; ++i)
			illuminaQualityBins[i] = 70;
		for(int i=72; i<sizeQual; ++i)
			illuminaQualityBins[i] = 73;
	};

	~TQualityMap(){
		delete[] phredIntToErrorMap;
		delete[] qualityToErrorMap;
		delete[] illuminaQualityBins;
	};

	double phredIntToError(int phredInt){
		if(phredInt >= maxPhred)
			return minError;
		else return phredIntToErrorMap[phredInt];
	};

	inline double phredToError(double phred){
		return pow(10.0, -phred/10.0);
	};

	double qualityToError(int qual){
		return phredIntToError(qual - 33);
	};

	int phredIntToQuality(int phredInt){
		return phredInt + 33;
	};

	inline int errorToPhredInt(const double & errorRate){
		return round(errorToPhred(errorRate));
	};

	inline double errorToPhred(const double & errorRate){
		if(errorRate < minError)
			return maxPhred;
		else
			return -10.0 * log10(errorRate);
	};

	inline int errorToQuality(const double & errorRate){
		return round(-10.0 * log10(errorRate)) + 33;
	};

	double& operator[](int phred){
		//Note: no check on range!
		return phredIntToErrorMap[phred];
	};

	double& operator[](uint8_t & phred){
		//Note: no check on range!
		return phredIntToErrorMap[phred];
	};
};


*/

#endif /* TGENOTYPEMAP_H_ */
