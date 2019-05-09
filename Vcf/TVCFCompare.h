/*
 * TCompareVCF.h
 *
 *  Created on: May 7, 2019
 *      Author: wegmannd
 */

#ifndef VCF_TVCFCOMPARE_H_
#define VCF_TVCFCOMPARE_H_

#include "TGenotypeMap.h"
#include "TLog.h"
#include "TParameters.h"
#include "TVcfFile.h"
#include "TFile.h"

class TGenotypeComparisonTable{
private:
	int** counts;
	int size;
	int missingIndex;
	int firstDiploidIndex;

public:
	TGenotypeComparisonTable();
	~TGenotypeComparisonTable();

	//add haploid genotypes
	void add(const Base b1, const Base b2);
	void addOtherMissing(const int sample, const Base b);
	void addFirstMissing(const Base b2);
	void addSecondMissing(const Base b1);

	//add diploid genotypes
	void add(const Genotype g1, const Genotype g2);
	void addOtherMissing(const int sample, const Genotype g);
	void addFirstMissing(const Genotype g2);
	void addSecondMissing(const Genotype g1);

	//add haploid / diploid combination of genotypes
	void add(const Genotype g1, const Base b2);
	void add(const Base b1, const Genotype g2);

	//write output
	void write(const std::string filename, TGenotypeMap & genoMap);
};

class TVCFComapreVCF{
private:
	int sampleIndex;
	std::vector<std::string> parsedChromosomes;
	TVcfFileSingleLine* vcfFile;
	bool vcfFileOpen;

public:

	TVCFComapreVCF(std::string & filename, std::string & sampleName, TLog* logfile);
	TVCFComapreVCF(TVCFComapreVCF&& other);
	TVCFComapreVCF& operator=(TVCFComapreVCF&& other);
	~TVCFComapreVCF();


	void next();
	void next(const int mindepth, const double minQual);
	void depthFilter(const int minDepth);
	void genotypeQualityFilter(const double minQual);

	bool eof(){ return vcfFile->eof; };
	bool isMissing(){ return vcfFile->sampleIsMissing(sampleIndex); };
	bool isDiploid(){ return vcfFile->sampleIsDiploid(sampleIndex); };
	Genotype genotype(TGenotypeMap & genoMap){ return vcfFile->sampleGenotype(sampleIndex, genoMap); };
	Base base(TGenotypeMap & genoMap){ return vcfFile->getFirstAlleleOfSample(sampleIndex, genoMap); };
	std::string chr(){ return vcfFile->chr(); };
	long position(){ return vcfFile->position(); };
	bool chrParsed(const std::string chr);
};

class TVCFCompare{
private:
	TLog* logfile;
	TGenotypeMap genoMap;

	std::vector<TVCFComapreVCF> vcfFiles;

	void addToOtherMissing(TGenotypeComparisonTable & counts, const int sample);

public:
	TVCFCompare(TLog* Logfile);
	~TVCFCompare(){};

	void compareVCFFiles(TParameters & parameters);
};




#endif /* VCF_TVCFCOMPARE_H_ */
