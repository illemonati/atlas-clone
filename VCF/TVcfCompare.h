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



//--------------------------------------------------------------
// Tgenometools::GenotypeComparisonTable
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
	void add(const genometools::Base b1, const genometools::Base b2);
	void addOtherMissing(const int sample, const genometools::Base b);
	void addFirstMissing(const genometools::Base b2);
	void addSecondMissing(const genometools::Base b1);

	//add diploid genotypes
	void add(const genometools::Genotype g1, const genometools::Genotype g2);
	void addOtherMissing(const int sample, const genometools::Genotype g);
	void addFirstMissing(const genometools::Genotype g2);
	void addSecondMissing(const genometools::Genotype g1);

	//add haploid / diploid combination of genotypes
	void add(const genometools::Genotype g1, const genometools::Base b2);
	void add(const genometools::Base b1, const genometools::Genotype g2);

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
	genometools::TVcfFileSingleLine* vcfFile;
	bool vcfFileOpen;

	int minDepth;
	double minQual;

public:
	TVcfComapreVCF();
	TVcfComapreVCF(std::string & filename, std::string & sampleName, coretools::TLog* logfile);
	TVcfComapreVCF(TVcfComapreVCF&& other);
	TVcfComapreVCF& operator=(TVcfComapreVCF&& other);
	~TVcfComapreVCF();


	void next();
	void setFilters(const int mindepth, const double minQual);

	bool eof(){ return vcfFile->eof; };
	bool isMissing(){ return vcfFile->sampleIsMissing(sampleIndex); };
	bool isDiploid(){ return vcfFile->sampleIsDiploid(sampleIndex); };
	genometools::Genotype genotype(){ return vcfFile->sampleGenotype(sampleIndex); };
	genometools::Base base(){ return vcfFile->getFirstAlleleOfSample(sampleIndex); };
	std::string chr(){ return vcfFile->chr(); };
	long position(){ return vcfFile->position(); };
	bool chrParsed(const std::string chr);
};

//--------------------------------------------------------------
// TVCFCompare
//--------------------------------------------------------------
class TVcfCompare{
private:
	coretools::TLog* logfile;
	std::vector<TVcfComapreVCF> vcfFiles;

	void addToOtherMissing(TGenotypeComparisonTable & counts, const int sample);

public:
	TVcfCompare(coretools::TLog* Logfile);
	~TVcfCompare(){};

	void compareVCFFiles(coretools::TParameters & parameters);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_VcfCompare : public coretools::TTask {
public:
	TTask_VcfCompare() { _explanation = "Comparing genotype calls in two VCF files"; };

	void run() {
		using namespace coretools::instances;
		TVcfCompare compare(&logfile());
		compare.compareVCFFiles(parameters());
	};
};

}; //end namespace

#endif /* VCF_TVCFCOMPARE_H_ */
