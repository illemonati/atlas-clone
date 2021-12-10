/*
 * TCompareVCF.h
 *
 *  Created on: May 7, 2019
 *      Author: wegmannd
 */

#ifndef VCF_TVCFCOMPARE_H_
#define VCF_TVCFCOMPARE_H_

#include "TLog.h"
#include "TParameters.h"
#include "TFile.h"
#include "TVcfFile.h"
#include "TTask.h"

namespace VCF{

using genometools::Base;
using genometools::Genotype;

using coretools::TLog;
using coretools::TParameters;


//--------------------------------------------------------------
// TGenotypeComparisonTable
//--------------------------------------------------------------
class TGenotypeComparisonTable{
private:
	static constexpr uint8_t size = 15;
	static constexpr uint8_t  missingIndex = 14;
	static constexpr uint8_t firstDiploidIndex = 4;

	std::array<std::array<uint32_t, size>, size> counts;

public:
	TGenotypeComparisonTable();
	~TGenotypeComparisonTable() = default;

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
	void write(const std::string filename);
};

//--------------------------------------------------------------
// TVCFComapreVCF
//--------------------------------------------------------------
class TVcfComapreVCF{
private:
	int sampleIndex;
	std::vector<std::string> parsedChromosomes;
	TVcfFileSingleLine* vcfFile;
	bool vcfFileOpen;

	int minDepth;
	double minQual;

public:
	TVcfComapreVCF();
	TVcfComapreVCF(std::string & filename, std::string & sampleName, TLog* logfile);
	TVcfComapreVCF(TVcfComapreVCF&& other);
	TVcfComapreVCF& operator=(TVcfComapreVCF&& other);
	~TVcfComapreVCF();


	void next();
	void setFilters(const int mindepth, const double minQual);

	bool eof(){ return vcfFile->eof; };
	bool isMissing(){ return vcfFile->sampleIsMissing(sampleIndex); };
	bool isDiploid(){ return vcfFile->sampleIsDiploid(sampleIndex); };
	Genotype genotype(){ return vcfFile->sampleGenotype(sampleIndex); };
	Base base(){ return vcfFile->getFirstAlleleOfSample(sampleIndex); };
	std::string chr(){ return vcfFile->chr(); };
	long position(){ return vcfFile->position(); };
	bool chrParsed(const std::string chr);
};

//--------------------------------------------------------------
// TVCFCompare
//--------------------------------------------------------------
class TVcfCompare{
private:
	TLog* logfile;
	std::vector<TVcfComapreVCF> vcfFiles;

	void addToOtherMissing(TGenotypeComparisonTable & counts, const int sample);

public:
	TVcfCompare(TLog* Logfile);
	~TVcfCompare(){};

	void compareVCFFiles(TParameters & parameters);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_VcfCompare:public coretools::TTask{
public:
	TTask_VcfCompare(){ _explanation = "Comparing genotype calls in two VCF files"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TVcfCompare compare(Logfile);
		compare.compareVCFFiles(Parameters);
	};
};


}; //end namespace

#endif /* VCF_TVCFCOMPARE_H_ */
