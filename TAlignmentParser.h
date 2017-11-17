/*
 * TAlignmentParser.h
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#ifndef TALIGNMENTPARSER_H_
#define TALIGNMENTPARSER_H_

#include "bamtools/api/BamReader.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamSequenceDictionary.h"
#include "bamtools/utils/bamtools_fasta.h"
#include "TGenotypeMap.h"
#include "TReadGroups.h"
#include "TLog.h"
#include "TRecalibration.h"
#include "TPostMortemDamage.h"

class TAlignmentParser{
private:
	//variables
	unsigned int maxSize;
	TGenotypeMap genoMap;
	TQualityMap qualityMap;
	TReadGroups* readGroupTable;
	TLog* logfile;
	bool _keepDuplicates;

	//data
	bool initialized;
	bool parsed;
	bool changed;
	int* qualityOriginal; //Note: quality is char as int: quality = (int) bam.quality
	int* qualityRecalibrated;
	double* errorRates;

	//tmp variables
	unsigned int i;
	int d, k, p;
	std::string tmpString;
	std::string tmpString2;

	//functions
	inline int toQual(const char & q);
	void parseBasesQualities();
	void setDistancesFromEnds();
	void fillContext();

public:
	//alignment: goal is to make this private!
	BamTools::BamAlignment bamAlignment;

	//details
	int length;
	int chrNumber;
	std::string readGroup;
	int readGroupId;

	int32_t position;
	bool isReverseStrand;
	bool passedFilters;
	bool recalibrated;

	//per base data
	Base* base;
	char* baseAsChar; //to be removed, if possible
	BaseContext* context;
	int* quality; //pointer to qualities ot be used
	bool* aligned; //whether or not base is aligned to ref. Insertions are not aligned
	int* alignedPos;
	int* distFrom3Prime;
	int* distFrom5Prime;

	//soft clipped data
	uint8_t softClippedEntry; //0 means start, 1 means end of read
	int* softClippedLength;
	char** softClippedBase;
	char** softClippedQuality;

	//construction
	TAlignmentParser();
	TAlignmentParser(TReadGroups* readGroupTable, unsigned int MaxSize, TLog* Logfile);
	~TAlignmentParser(){
		clear();
	}
	void init(TReadGroups* readGroupTable, unsigned int MaxSize, TLog* Logfile);
	void keepDuplicates(){_keepDuplicates = true;};
	void clear();

	//functions to read and parse
	bool readAlignment(BamTools::BamReader & bamReader);
	void parse();

	//functions to access data
	std::string& name(){return bamAlignment.Filename;};
	void recalibrate(TRecalibration & recalObject);
	void recalibrate(TRecalibration & recalObject, TPMD* pmdObjects, BamTools::Fasta & reference);
	void save(BamTools::BamWriter & bamWriter);
	void print();
};



#endif /* TALIGNMENTPARSER_H_ */
