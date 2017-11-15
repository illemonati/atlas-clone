/*
 * TAlignmentParser.h
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#ifndef TALIGNMENTPARSER_H_
#define TALIGNMENTPARSER_H_

#include "bamtools/api/BamReader.h"
#include "bamtools/api/SamSequenceDictionary.h"
#include "TGenotypeMap.h"
#include "TReadGroups.h"
#include "TLog.h"

class TAlignmentParser{
private:
	//variables
	unsigned int maxSize;
	TGenotypeMap genoMap;
	TReadGroups* readGroupTable;
	TLog* logfile;
	bool _keepDuplicates;
	bool initialized;

	//tmp variables
	unsigned int i;
	int d, k, p;

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
	bool parsed;

	//per base data
	Base* base;
	char* baseAsChar; //to be removed, if possible
	BaseContext* context;
	int* quality;
	bool* aligned; //whether or not base is aligned to ref. Insertions are not aligned
	int* alignedPos;
	int* distFrom3Prime;
	int* distFrom5Prime;

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
	void print();

	//functions to access data
	std::string& name(){return bamAlignment.Filename;};
	void recalibrate(TRecalibration & recalObject);
};



#endif /* TALIGNMENTPARSER_H_ */
