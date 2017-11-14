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

class TAlignmentParser{
private:
	int maxSize;
	TGenotypeMap genoMap;
	TReadGroups* readGroupTable;
	bool initialized;

	std::string name;

	//tmp variables
	unsigned int i;
	int d, k, p;

	inline int toQual(const char & q);
	void parseBasesQualities(BamTools::BamAlignment & bamAlignment);
	void setDistancesFromEnds(BamTools::BamAlignment & bamAlignment);
	void fillContext();

public:
	//details
	int length;
	std::string readGroup;
	int readGroupId;
	int32_t position;
	bool isReverseStrand;

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
	TAlignmentParser(TReadGroups* readGroupTable, int Size);
	~TAlignmentParser(){
		clear();
	}
	void init(TReadGroups* readGroupTable, int MaxSize);
	void clear();

	//functions
	void parse(BamTools::BamAlignment & bamAlignment);
	void print();
};



#endif /* TALIGNMENTPARSER_H_ */
