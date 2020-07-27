/*
 * TGenotypeMap.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TGENOTYPEMAP_H_
#define TGENOTYPEMAP_H_

#include "stringFunctions.h"
#include "TBase.h"
#include <math.h>

enum Genotype : uint8_t {AA=0, AC, AG, AT, CC, CG, CT, GG, GT, TT, NN};
enum BaseContext : uint8_t {cAA=0, cAC, cAG, cAT, cCA, cCC, cCG, cCT, cGA, cGC, cGG, cGT, cTA, cTC, cTG, cTT, cNA, cNC, cNG, cNT, cAN, cCN, cGN, cTN, cNN}; //N means unknown base or "nothing", i.e. end of read

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
	Genotype** alleleicCombinationToGenotypes; //mapping the allelic combinations to the homozygous, heterozygosu and other homozygous genotype
	char** alleleicCombinationToBaseChar; //mapping the allelic combination to the two bases
	uint16_t** allelicCombinationsMatchingBase; //mapping all allelic combinations that contain a specific base

	char* baseToChar;
	Base* baseToFlippedBase;
	size_t numGenotypes;
	size_t numContexts;
	size_t numContextsNotN;

	TGenotypeMap(){
		//create genotype map
		genotypeMap = new Genotype*[4];
		for(size_t i=0; i<4; ++i)
			genotypeMap[i] = new Genotype[4];

		//fill genotype map
		int geno = 0;
		for(size_t i=0; i<4; ++i){
			for(size_t j=i; j<4; ++j){
				genotypeMap[i][j] = static_cast<Genotype>(geno);
				genotypeMap[j][i] = genotypeMap[i][j];
				++geno;
			}
		}

		//create and fill genotypeToBase
		numGenotypes = 10;
		genotypeToBase = new Base*[numGenotypes];
		for(size_t i=0; i<10; ++i){
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

		//fill allelicCombinationsMatchingBase: gives the three allelelic combinations containing a specific base
		allelicCombinationsMatchingBase = new uint16_t*[4];
		allelicCombinationsMatchingBase[A] = new uint16_t[3]; allelicCombinationsMatchingBase[A][0] = 0; allelicCombinationsMatchingBase[A][1] = 1; allelicCombinationsMatchingBase[A][2] = 2;
		allelicCombinationsMatchingBase[C] = new uint16_t[3]; allelicCombinationsMatchingBase[C][0] = 0; allelicCombinationsMatchingBase[C][1] = 3; allelicCombinationsMatchingBase[C][2] = 4;
		allelicCombinationsMatchingBase[G] = new uint16_t[3]; allelicCombinationsMatchingBase[G][0] = 1; allelicCombinationsMatchingBase[G][1] = 3; allelicCombinationsMatchingBase[G][2] = 5;
		allelicCombinationsMatchingBase[T] = new uint16_t[3]; allelicCombinationsMatchingBase[T][0] = 2; allelicCombinationsMatchingBase[T][1] = 4; allelicCombinationsMatchingBase[T][2] = 5;
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

		for(int i=0; i<4; ++i){
			delete[] allelicCombinationsMatchingBase[i];
		}
		delete[] allelicCombinationsMatchingBase;
	};

	bool isValidBase(const char base){
		if(base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'N' || base == 'a' || base == 'c' || base == 'g' || base == 't' || base == 'n')
			return true;
		return false;
	};

	Base toBase(const char base) const{
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

	Base getBaseOnlyCapitals(const char base) const{
		if(base == 'A') return A;
		if(base == 'C') return C;
		if(base == 'G') return G;
		if(base == 'T') return T;
		return N;
	};

	char getBaseAsChar(const Base base) const{
		if(base == A) return 'A';
		if(base == C) return 'C';
		if(base == G) return 'G';
		if(base == T) return 'T';
		return 'N';
	};

	char getBaseAsChar(const uint8_t base) const{
		if(base == A) return 'A';
		if(base == C) return 'C';
		if(base == G) return 'G';
		if(base == T) return 'T';
		return 'N';
	};

	Base flipBase(const char base) const{
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

	uint16_t getNumGenotypes() const{
		return numGenotypes;
	};

	Genotype toGenotype(const Base first, const Base second) const{
		if(first < second) return genotypeMap[first][second];
		else return genotypeMap[second][first];
	};

	Genotype toGenotype(const uint8_t first, const uint8_t second) const{
		if(first < second) return genotypeMap[first][second];
		else return genotypeMap[second][first];
	};

	Genotype toGenotype(const char first, const char second) const{
		Base Bfirst = toBase(first);
		Base Bsecond = toBase(second);
		if(Bfirst < Bsecond) return genotypeMap[Bfirst][Bsecond];
		else return genotypeMap[Bsecond][Bfirst];
	};

	std::string getGenotypeString(const uint8_t num) const{
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

	std::string getGenotypeString(const Genotype geno) const{
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

	std::string getGenotypeStringKnownAlleles(const uint8_t num, const char ref, const char alt) const{
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
			Base refBase = toBase(ref);
			Base altBase = toBase(alt);
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

	std::pair<Base,Base> getBasesOfGenotype(const uint8_t num) const{
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

	size_t getNumContextsNotN() const{
		return numContextsNotN;
	};

	size_t getNumContexts() const{
		return numContexts;
	};

	BaseContext toContext(const Base first, const Base second) const{
		if(second == N) throw "Context not defined with second base = N!";
		return contextMap[first][second];
	};

	BaseContext toContext(const uint8_t first, const uint8_t second) const{
		if(second > 3) throw "Context not defined with second base = N!";
		return contextMap[first][second];
	};

	BaseContext toContext(const char first, const char second) const{
		return toContext(toBase(first), toBase(second));
	};

	BaseContext getContextReverseRead(const char first, const char second) const{
		return toContext(flipBase(first), flipBase(second));
	};

	std::string getContextString(const uint8_t context) const{
		return getContextString(static_cast<BaseContext>(context));
	};

	std::string getContextString(BaseContext context) const{
		if(context == cAA) return "AA";
		if(context == cAC) return "AC";
		if(context == cAG) return "AG";
		if(context == cAT) return "AT";
		if(context == cCA) return "CA";
		if(context == cCC) return "CC";
		if(context == cCG) return "CG";
		if(context == cCT) return "CT";
		if(context == cGA) return "GA";
		if(context == cGC) return "GC";
		if(context == cGG) return "GG";
		if(context == cGT) return "GT";
		if(context == cTA) return "TA";
		if(context == cTC) return "TC";
		if(context == cTG) return "TG";
		if(context == cTT) return "TT";
		if(context == cNA) return "NA";
		if(context == cNC) return "NC";
		if(context == cNG) return "NG";
		if(context == cNT) return "NT";
		if(context == cNA) return "NA";
		if(context == cCN) return "NC";
		if(context == cGN) return "NG";
		if(context == cTN) return "NT";
		if(context == cAN) return "NA";
		else return "NN";
	};

	void addContextNames(std::vector<std::string> & vec){
		for(uint8_t c=0; c<numContexts; ++c){
			vec.push_back(getContextString(c));
		}
	};
};

#endif /* TGENOTYPEMAP_H_ */
