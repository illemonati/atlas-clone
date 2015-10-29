/*
 * TGenotypeMap.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TGENOTYPEMAP_H_
#define TGENOTYPEMAP_H_

#include "stringFunctions.h"

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

		//create and fill context map
		contextMap = new BaseContext*[5];
		int context = 0;
		for(int i=0; i<5; ++i){
			contextMap[i] = new BaseContext[4];
			for(int j=0; j<4; ++j){
				contextMap[i][j] = static_cast<BaseContext>(context);
				++context;
			}
		}
	};

	~TGenotypeMap(){
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

	Genotype getGenotype(Base first, Base second){
		if(first < second) return genotypeMap[first][second];
		else return genotypeMap[second][first];
	};

	Genotype getGenotype(int first, int second){
		if(first < second) return genotypeMap[first][second];
		else return genotypeMap[second][first];
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

	TBaseFrequencies(){
		for(int i = 0; i < 4; ++i) freq[i] = 0.0;
	};
	void add(Base B, double & weight){
		freq[B] += weight;
	}
	void normalize(){
		double sum = 0.0;
		for(int i = 0; i < 4; ++i) sum += freq[i];
		sum += 4.0;
		for(int i = 0; i < 4; ++i) freq[i] = (freq[i] + 1.0) / sum;
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


#endif /* TGENOTYPEMAP_H_ */
