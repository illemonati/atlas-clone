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

enum Base {A=0, C, G, T, N};
enum Genotype {AA=0, AC, AG, AT, CC, CG, CT, GG, GT, TT};
enum BaseContext {cAA=0, cAC, cAG, cAT, cCA, cCC, cCG, cCT, cGA, cGC, cGG, cGT, cTA, cTC, cTG, cTT, cNA, cNC, cNG, cNT}; //N means "nothing", i.e. end of read or del

//---------------------------------------------------------------
//GenotypeMap
//---------------------------------------------------------------
//genotype map for enum type
class TGenotypeMap{
public:
	Genotype** genotypeMap; //mapping base numbering to genotype enum
	BaseContext** contextMap; //mapping dinucleotide context to context enum
	Base** genotypeToBase; //mapping genotypes to bases
	char* baseToChar;
	int numContexts;

	TGenotypeMap(){
		//create genotype map
		genotypeMap = new Genotype*[4];
		for(int i=0; i<4; ++i)
			genotypeMap[i] = new Genotype[4];

		//fill genotype map
		int geno = 0;
		for(int i=0; i<4; ++i){
			for(int j=i; j<4; ++j){
				genotypeMap[i][j] = static_cast<Genotype>(geno);
				genotypeMap[j][i] = genotypeMap[i][j];
				++geno;
			}
		}

		//create and fill genotypeToBase
		genotypeToBase = new Base*[10];
		for(int i=0; i<10; ++i){
			genotypeToBase[i] = new Base[2];
		}
		genotypeToBase[0][0] = A; genotypeToBase[0][1] = A;
		genotypeToBase[1][0] = A; genotypeToBase[1][1] = C;
		genotypeToBase[2][0] = A; genotypeToBase[2][1] = G;
		genotypeToBase[3][0] = A; genotypeToBase[3][1] = T;
		genotypeToBase[4][0] = C; genotypeToBase[4][1] = C;
		genotypeToBase[5][0] = C; genotypeToBase[5][1] = G;
		genotypeToBase[6][0] = C; genotypeToBase[6][1] = T;
		genotypeToBase[7][0] = G; genotypeToBase[7][1] = G;
		genotypeToBase[8][0] = G; genotypeToBase[8][1] = T;
		genotypeToBase[9][0] = T; genotypeToBase[9][1] = T;

		//create and fill context map
		numContexts = 20;
		contextMap = new BaseContext*[5];
		int context = 0;
		for(int i=0; i<5; ++i){
			contextMap[i] = new BaseContext[4];
			for(int j=0; j<4; ++j){
				contextMap[i][j] = static_cast<BaseContext>(context);
				++context;
			}
		}

		//fill base to char map
		baseToChar = new char[5];
		baseToChar[A] = 'A';
		baseToChar[C] = 'C';
		baseToChar[G] = 'G';
		baseToChar[T] = 'T';
		baseToChar[N] = 'N';
	};

	~TGenotypeMap(){
		for(int i=0; i<4; ++i){
			delete[] genotypeMap[i];
		}
		for(int i=0; i<5; ++i){
			delete[] contextMap[i];
		}
		delete[] genotypeMap;
		for(int i=0; i<10; ++i)
			delete[] genotypeToBase[i];
		delete[] genotypeToBase;
		delete[] contextMap;
		delete[] baseToChar;
	};

	Base getBase(char & base){
		if(base == 'A') return A;
		if(base == 'C') return C;
		if(base == 'G') return G;
		if(base == 'T') return T;
		if(base == 'a') return A;
		if(base == 'c') return C;
		if(base == 'g') return G;
		if(base == 't') return T;
		return N;
	};

	Base getBaseOnlyCapitals(char & base){
		if(base == 'A') return A;
		if(base == 'C') return C;
		if(base == 'G') return G;
		if(base == 'T') return T;
		return N;
	};

	char getBaseAsChar(Base base){
		if(base == A) return 'A';
		if(base == C) return 'C';
		if(base == G) return 'G';
		if(base == T) return 'T';
		return 'N';
	};

	char getBaseAsChar(int base){
		if(base == A) return 'A';
		if(base == C) return 'C';
		if(base == G) return 'G';
		if(base == T) return 'T';
		return 'N';
	};

	Base flipBase(char & base){
		if(base == 'A') return T;
		if(base == 'C') return G;
		if(base == 'G') return C;
		if(base == 'T') return A;
		if(base == 'a') return T;
		if(base == 'c') return G;
		if(base == 'g') return C;
		if(base == 't') return A;
		return N;
	};

	Genotype getGenotype(Base first, Base second){
		if(first < second) return genotypeMap[first][second];
		else return genotypeMap[second][first];
	};

	Genotype getGenotype(int first, int second){
		if(first < second) return genotypeMap[first][second];
		else return genotypeMap[second][first];
	};

	Genotype getGenotype(char first, char second){
		Base Bfirst = getBase(first);
		Base Bsecond = getBase(second);
		if(Bfirst < Bsecond) return genotypeMap[Bfirst][Bsecond];
		else return genotypeMap[Bsecond][Bfirst];
	};

	std::string getGenotypeString(int num){
		if(num==0) return "AA";
		if(num==1) return "AC";
		if(num==2) return "AG";
		if(num==3) return "AT";
		if(num==4) return "CC";
		if(num==5) return "CG";
		if(num==6) return "CT";
		if(num==7) return "GG";
		if(num==8) return "GT";
		if(num==9) return "TT";
		throw "GenotypeMap: Unknown genotype with number " + toString(num) + "!";
	};

	std::pair<Base,Base> getBasesOfGenotype(int num){
		if(num==0) return std::pair<Base,Base>(A,A);
		if(num==1) return std::pair<Base,Base>(A,C);
		if(num==2) return std::pair<Base,Base>(A,G);
		if(num==3) return std::pair<Base,Base>(A,T);
		if(num==4) return std::pair<Base,Base>(C,C);
		if(num==5) return std::pair<Base,Base>(C,G);
		if(num==6) return std::pair<Base,Base>(C,T);
		if(num==7) return std::pair<Base,Base>(G,G);
		if(num==8) return std::pair<Base,Base>(G,T);
		if(num==9) return std::pair<Base,Base>(T,T);
		throw "GenotypeMap: Unknown genotype with number " + toString(num) + "!";
	}

	int getNumContext(){
		return 20;
	};

	BaseContext getContext(Base first, Base second){
		if(second == N) throw "Context not defined with second base = N!";
		return contextMap[first][second];
	};

	BaseContext getContext(int first, int second){
		if(second > 3) throw "Context not defined with second base = N!";
		return contextMap[first][second];
	};

	BaseContext getContext(char first, char second){
		return getContext(getBase(first), getBase(second));
	};

	BaseContext getContextReverseRead(char first, char second){
		return getContext(flipBase(first), flipBase(second));
	};

	std::string getContextString(int num){
		if(num==0) return "AA";
		if(num==1) return "AC";
		if(num==2) return "AG";
		if(num==3) return "AT";
		if(num==4) return "CA";
		if(num==5) return "CC";
		if(num==6) return "CG";
		if(num==7) return "CT";
		if(num==8) return "GA";
		if(num==9) return "GC";
		if(num==10) return "GG";
		if(num==11) return "GT";
		if(num==12) return "TA";
		if(num==13) return "TC";
		if(num==14) return "TG";
		if(num==15) return "TT";
		if(num==16) return "-A";
		if(num==17) return "-C";
		if(num==18) return "-G";
		if(num==19) return "-T";
		throw "GenotypeMap: Unknown text with number " + toString(num) + "!";
	};
};

//---------------------------------------------------------------
//TBaseFrequencies
//---------------------------------------------------------------
class TBaseFrequencies{
public:
	double freq[4];
	bool wasNormalized;

	TBaseFrequencies(){
		for(int i = 0; i < 4; ++i) freq[i] = 0.0;
		wasNormalized = false;
	};
	void add(Base B, double & weight){
		freq[B] += weight;
	};
	void addNoRef(Base B, double weight){
		freq[B] += weight;
	};
	void normalize(){
		if(!wasNormalized){
			double sum = 0.0;
			for(int i = 0; i < 4; ++i) sum += freq[i];
			sum += 4.0;
			for(int i = 0; i < 4; ++i) freq[i] = (freq[i] + 1.0) / sum;
			wasNormalized = true;
		}
	};
	void clear(){
		for(int i = 0; i < 4; ++i) freq[i] = 0.0;
		wasNormalized = false;
	};
	void print(){
		std::cout << "freq(A) = " << freq[0] << ", freq(C) = " << freq[1] << ", freq(G) = " << freq[2] << ", freq(T) = " << freq[3] << std::endl;
	};
	double& operator[](int pos){
		return freq[pos];
	};
};

//---------------------------------------------------------------
//TQualityMap
//---------------------------------------------------------------
class TQualityMap{
public:
	double* phredToErrorMap;
	double* charToErrorMap;
	double min;
	int size;

	TQualityMap(){
		//only up to phred = 255, else always return 256
		size = 256;
		phredToErrorMap = new double[size];
		charToErrorMap = new double[size + 33];

		//initialize quality <= 0
		for(int i=0; i<33; ++i)
			charToErrorMap[i] = 1.0;
		phredToErrorMap[0] = 1.0;

		//and now others
		for(int i=0; i<256; ++i){
			phredToErrorMap[i] = pow(10.0, (double) -i/10.0);
			charToErrorMap[i+33] = phredToErrorMap[i];
		}
		min = pow(10.0, (double) -256/10.0);
	};

	~TQualityMap(){
		delete[] phredToErrorMap;
		delete[] charToErrorMap;
	};

	double phredToError(int phred){
		if(phred>255)
			return min;
		else return phredToErrorMap[phred];
	};

	inline int errorToPhred(const double & errorRate){
		return round(-10.0 * log10(errorRate));
	};

	inline int errorToQuality(const double & errorRate){
		return round(-10.0 * log10(errorRate)) + 33;
	};

	double& operator[](int phred){
		//Note: no check on range!
		return phredToErrorMap[phred];
	};
};

#endif /* TGENOTYPEMAP_H_ */
