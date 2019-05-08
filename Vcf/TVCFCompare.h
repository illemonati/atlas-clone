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

class TGenotypeComparisonTable{
private:
	int** counts;
	int size;
	int missingIndex;
	TGenotypeMap genoMap;

public:
	TGenotypeComparisonTable();

	//add genotypes
	void add(const Genotype g1, const Genotype g2);
	void addFirstMissing(const Genotype g2);
	void addSecondMissing(const Genotype g1);

	//add from bases
	void add(const char ind1_first, const char ind1_second, const char ind2_first, const char ind2_second);
	void addFirstMissing(const char ind2_first, const char ind2_second);
	void addSecondMissing(const char ind1_first, const char ind1_second);

	void add(TVcfFileSingleLine & vcfFile1, TVcfFileSingleLine & vcfFile2);
	void addFirstMissing(TVcfFileSingleLine & vcfFile2);
	void addSecondMissing(TVcfFileSingleLine & vcfFile1);
};


class TVCFCompare{
private:
	TLog* logfile;

	void openVCF(std::string & filename, TVcfFileSingleLine & vcfFile);

public:
	TVCFCompare(TLog* Logfile);
	~TVCFCompare();

	void compareVCFFiles(TParameters & parameters);
};




#endif /* VCF_TVCFCOMPARE_H_ */
