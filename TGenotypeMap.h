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
	BaseContext** contextMapFlipped; //mapping dinucleotide context to context enum
	BaseContext flippedContext[25] = {cTT, cTG, cTC, cTA, cGT, cGG, cGC, cGA, cCT, cCG, cCC, cCA, cAT, cAG, cAC, cAA, cNT, cNG, cNC, cNA, cTN, cGN, cCN, cAN, cNN}; //N means unknwon base or "nothing", i.e. end of read or del
	Base** genotypeToBase; //mapping genotypes to bases
	Base** alleleicCombinationToBase; //mapping the allelic combination to the two bases
	Genotype** alleleicCombinationToGenotypes; //mapping the alleleic combinations to the homozygous, heterozygosu and other homozygous genotype
	char** alleleicCombinationToBaseChar; //mapping the allelic combination to the two bases

	char* baseToChar;
	Base* baseToFlippedBase;
	int numGenotypes;
	int numContexts;
	int numContextsNotN;

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
		numGenotypes = 10;
		genotypeToBase = new Base*[numGenotypes];
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
		numContexts = 25;
		numContextsNotN = 20;
		contextMap = new BaseContext*[5];

		//now fill regular context
		int context = 0;
		for(int i=0; i<5; ++i){
			contextMap[i] = new BaseContext[5];
			for(int j=0; j<4; ++j){
				contextMap[i][j] = static_cast<BaseContext>(context);
				++context;
			}
		}

		//Now add those that should not occur, but sometimes do in bam files
		//Note that these should never occur in our data processing as they imply the base is N
		contextMap[0][4] = cAN;
		contextMap[1][4] = cCN;
		contextMap[2][4] = cGN;
		contextMap[3][4] = cTN;
		contextMap[4][4] = cNN;

		//Same for FLIPPED context
		contextMapFlipped = new BaseContext*[5];
		context = 0;
		for(int i=0; i<5; ++i){
			contextMapFlipped[i] = new BaseContext[5];
			for(int j=0; j<4; ++j){
				contextMapFlipped[i][j] = flippedContext[context];
				++context;
			}
		}
		contextMapFlipped[0][4] = cTN;
		contextMapFlipped[1][4] = cGN;
		contextMapFlipped[2][4] = cCN;
		contextMapFlipped[3][4] = cAN;
		contextMapFlipped[4][4] = cNN;


		//fill base to char map
		baseToChar = new char[5];
		baseToChar[A] = 'A';
		baseToChar[C] = 'C';
		baseToChar[G] = 'G';
		baseToChar[T] = 'T';
		baseToChar[N] = 'N';

		//fill baseToFlippedBase map
		baseToFlippedBase = new Base[5];
		baseToFlippedBase[A] = T;
		baseToFlippedBase[C] = G;
		baseToFlippedBase[G] = C;
		baseToFlippedBase[T] = A;
		baseToFlippedBase[N] = N;

		//fill alleleicCombinationToBase
		alleleicCombinationToBase = new Base*[6];
		alleleicCombinationToBase[0] = new Base[2]; alleleicCombinationToBase[0][0] = A; alleleicCombinationToBase[0][1] = C;
		alleleicCombinationToBase[1] = new Base[2]; alleleicCombinationToBase[1][0] = A; alleleicCombinationToBase[1][1] = G;
		alleleicCombinationToBase[2] = new Base[2]; alleleicCombinationToBase[2][0] = A; alleleicCombinationToBase[2][1] = T;
		alleleicCombinationToBase[3] = new Base[2]; alleleicCombinationToBase[3][0] = C; alleleicCombinationToBase[3][1] = G;
		alleleicCombinationToBase[4] = new Base[2]; alleleicCombinationToBase[4][0] = C; alleleicCombinationToBase[4][1] = T;
		alleleicCombinationToBase[5] = new Base[2]; alleleicCombinationToBase[5][0] = G; alleleicCombinationToBase[5][1] = T;

		//fill alleleicCombinationToBaseChar
		alleleicCombinationToBaseChar = new char*[6];
		alleleicCombinationToBaseChar[0] = new char[2]; alleleicCombinationToBaseChar[0][0] = 'A'; alleleicCombinationToBaseChar[0][1] = 'C';
		alleleicCombinationToBaseChar[1] = new char[2]; alleleicCombinationToBaseChar[1][0] = 'A'; alleleicCombinationToBaseChar[1][1] = 'G';
		alleleicCombinationToBaseChar[2] = new char[2]; alleleicCombinationToBaseChar[2][0] = 'A'; alleleicCombinationToBaseChar[2][1] = 'T';
		alleleicCombinationToBaseChar[3] = new char[2]; alleleicCombinationToBaseChar[3][0] = 'C'; alleleicCombinationToBaseChar[3][1] = 'G';
		alleleicCombinationToBaseChar[4] = new char[2]; alleleicCombinationToBaseChar[4][0] = 'C'; alleleicCombinationToBaseChar[4][1] = 'T';
		alleleicCombinationToBaseChar[5] = new char[2]; alleleicCombinationToBaseChar[5][0] = 'G'; alleleicCombinationToBaseChar[5][1] = 'T';

		//fill alleleicCombinationToGenotypes (hom, het, hom2)
		alleleicCombinationToGenotypes = new Genotype*[6];
		alleleicCombinationToGenotypes[0] = new Genotype[3]; alleleicCombinationToGenotypes[0][0] = AA; alleleicCombinationToGenotypes[0][1] = AC; alleleicCombinationToGenotypes[0][2] = CC;
		alleleicCombinationToGenotypes[1] = new Genotype[3]; alleleicCombinationToGenotypes[1][0] = AA; alleleicCombinationToGenotypes[1][1] = AG; alleleicCombinationToGenotypes[1][2] = GG;
		alleleicCombinationToGenotypes[2] = new Genotype[3]; alleleicCombinationToGenotypes[2][0] = AA; alleleicCombinationToGenotypes[2][1] = AT; alleleicCombinationToGenotypes[2][2] = TT;
		alleleicCombinationToGenotypes[3] = new Genotype[3]; alleleicCombinationToGenotypes[3][0] = CC; alleleicCombinationToGenotypes[3][1] = CG; alleleicCombinationToGenotypes[3][2] = GG;
		alleleicCombinationToGenotypes[4] = new Genotype[3]; alleleicCombinationToGenotypes[4][0] = CC; alleleicCombinationToGenotypes[4][1] = CT; alleleicCombinationToGenotypes[4][2] = TT;
		alleleicCombinationToGenotypes[5] = new Genotype[3]; alleleicCombinationToGenotypes[5][0] = GG; alleleicCombinationToGenotypes[5][1] = GT; alleleicCombinationToGenotypes[5][2] = TT;
	};

	~TGenotypeMap(){
		for(int i=0; i<4; ++i){
			delete[] genotypeMap[i];
		}
		for(int i=0; i<5; ++i){
			delete[] contextMap[i];
			delete[] contextMapFlipped[i];
		}
		delete[] genotypeMap;
		for(int i=0; i<10; ++i)
			delete[] genotypeToBase[i];
		delete[] genotypeToBase;
		delete[] contextMap;
		delete[] contextMapFlipped;
		delete[] baseToChar;
		delete[] baseToFlippedBase;

		for(int i=0; i<6; ++i){
			delete[] alleleicCombinationToBase[i];
			delete[] alleleicCombinationToBaseChar[i];
			delete[] alleleicCombinationToGenotypes[i];
		}
		delete[] alleleicCombinationToBase;
		delete[] alleleicCombinationToBaseChar;
		delete[] alleleicCombinationToGenotypes;
	};

	Base getBase(const char base){
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

	Base getBaseOnlyCapitals(const char base){
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

	std::string getGenotypeString(Genotype geno){
		if(geno==0) return "AA";
		if(geno==1) return "AC";
		if(geno==2) return "AG";
		if(geno==3) return "AT";
		if(geno==4) return "CC";
		if(geno==5) return "CG";
		if(geno==6) return "CT";
		if(geno==7) return "GG";
		if(geno==8) return "GT";
		if(geno==9) return "TT";
		throw "GenotypeMap: Unknown genotype with number " + toString((int) geno) + "!";
	};

	std::string getGenotypeStringKnownAlleles(int num, char ref, char alt){
		std::string genotype = "";
		if(num == 0){
			genotype += ref;
			genotype += ref;
			return genotype;
		}
		else if(num == 2){
			genotype += alt;
			genotype += alt;
			return genotype;
		}
		else if(num == 1){
			Base refBase = getBase(ref);
			Base altBase = getBase(alt);
			if(refBase > altBase){
				genotype += alt;
				genotype += ref;
				return genotype;
			}
			else{
				genotype += ref;
				genotype += alt;
				return genotype;
			}
		}
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
	void setEqualBaseFreq(){
		for(int i = 0; i < 4; ++i) freq[i] = 0.25;
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
