/*
 * TAlignment.h
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#ifndef TALIGNMENT_H_
#define TALIGNMENT_H_



class TAlignment{
private:
	//details
	bool empty;
	int length;
	int chrNumber;
	std::string readGroup;
	int readGroupId;
	int32_t position;
	bool isReverseStrand;
	bool isProperPair;
	bool passedFilters;
	bool recalibrated;

	//data
	unsigned int maxSize;
	bool initialized;
	bool parsed;
	bool changed;
	int* qualityOriginal; //Note: quality is char as int: quality = (int) bam.quality
	int* qualityRecalibrated;
	double* errorRates;

	Base* base;
	char* baseAsChar; //TODO: to be removed, if possible
	BaseContext* context;
	int* quality; //pointer to qualities to be used
	bool* aligned; //whether or not base is aligned to ref. Insertions are not aligned
	int* alignedPos;
	int* distFrom3Prime;
	int* distFrom5Prime;
	double* pmdCT;
	double* pmdGA;


	//reference
	bool hasReference;
	std::string referenceSequence;

	//soft clipped data
	uint8_t softClippedEntry; //0 means start, 1 means end of read
	int* softClippedLength;
	char** softClippedBase;
	char** softClippedQuality;

	//functions
	void initStorage();
	void freeStorage();
	void clear();

	void parseBasesQualities();
	void setDistancesFromEnds();
	void fillContext();


public:
	TAlignment(int size);

	void fill(BamTools::BamAlignment & bamAlignment, int ReadGroupId);
	void setFiltersPassed(bool passed);
};

#endif /* TALIGNMENT_H_ */
