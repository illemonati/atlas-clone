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
	int size = 11;

public:
	TGenotypeComparisonTable();

	void add();
};


class TVCFCompare{
private:
	TGenotypeComparisonTable countTable;
	TLog* logfile;

	void openVCF(std::string & filename, TVcfFileSingleLine & vcfFile);

public:
	TVCFCompare(TLog* Logfile);
	~TVCFCompare();

	void compareVCFFiles(TParameters & parameters);
};




#endif /* VCF_TVCFCOMPARE_H_ */
