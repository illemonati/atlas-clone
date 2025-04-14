/*
 * TCompareVCF.h
 *
 *  Created on: May 7, 2019
 *      Author: wegmannd
 */

#ifndef VCF_TVCFCOMPARE_H_
#define VCF_TVCFCOMPARE_H_

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "genometools/VCF/TVcfFile.h"

namespace VCF{

//--------------------------------------------------------------
// Tgenometools::GenotypeComparisonTable
//--------------------------------------------------------------
class TGenotypeComparisonTable{
private:
	static constexpr uint8_t _size = 15;
	static constexpr uint8_t _firstDiploidIndex = 4;

	std::array<std::array<uint32_t, _size>, _size> _counts{};

public:
	//add haploid genotypes
	void add(const genometools::Base b1, const genometools::Base b2);
	void addOtherMissing(const int sample, const genometools::Base b);

	//add diploid genotypes
	void add(const genometools::Genotype g1, const genometools::Genotype g2);
	void addOtherMissing(const int sample, const genometools::Genotype g);

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
	int _sampleIndex = 0;
	std::vector<std::string> _parsedChromosomes;
	std::unique_ptr<genometools::TVcfFileSingleLine> _vcfFile;

	int _minDepth = 0;
	double _minQual = 0.;

public:
	TVcfComapreVCF(std::string_view filename, std::string_view sampleName);

	void next();
	void setFilters(const int mindepth, const double minQual);

	bool eof(){ return _vcfFile->eof; };
	bool isMissing(){ return _vcfFile->sampleIsMissing(_sampleIndex); };
	bool isDiploid(){ return _vcfFile->sampleIsDiploid(_sampleIndex); };
	genometools::Genotype genotype(){ return _vcfFile->sampleGenotype(_sampleIndex); };
	genometools::Base base(){ return _vcfFile->getFirstAlleleOfSample(_sampleIndex); };
	std::string chr(){ return _vcfFile->chr(); };
	long position(){ return _vcfFile->position(); };
	bool chrParsed(const std::string chr);
};

//--------------------------------------------------------------
// TVCFCompare
//--------------------------------------------------------------
class TVcfCompare{
private:
	std::vector<TVcfComapreVCF> _vcfFiles;
	bool _limitLines = false;
	long _lineLimit  = -1;
	std::string _outName;

	void addToOtherMissing(TGenotypeComparisonTable & counts, const int sample);

public:
	TVcfCompare();
	void run();
};

} // namespace VCF

#endif /* VCF_TVCFCOMPARE_H_ */
