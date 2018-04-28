/*
 * TAlignment.h
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#ifndef TALIGNMENT_H_
#define TALIGNMENT_H_

#include "stringFunctions.h"
#include "TBase.h"
#include "TPostMortemDamage.h"
#include "TRecalibration.h"
#include "bamtools/api/BamAlignment.h"
#include "bamtools/utils/bamtools_fasta.h"
#include "bamtools/api/BamWriter.h"

//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
//a buffer class to speed up adding the reference sequence to each read
//This class makes use of the fact that bam files are sorted, hence the buffer can always start at the current position

class TFastaBuffer{
private:
	BamTools::Fasta* reference;
	int bufferSize;
	std::string referenceSequence;

	int curChr;
	long curStart, curEnd;

	void moveTo(const int & chr, const int32_t & pos);

public:
	TFastaBuffer(BamTools::Fasta* Reference);
	void fill(const int & chr, const int32_t & start, const int32_t end, std::string & ref);
};

class TAlignmentParser;

//-----------------------------------------------------
//TAlignment
//-----------------------------------------------------

class TAlignment{
private:
	//details
	bool empty;
	std::string readGroup;
	bool recalibrated;

	//data
	BamTools::BamAlignment bamAlignment;
	unsigned int maxSize;
	int length;

	int* alignedPos;
	bool storageInitialized;
	bool parsed;
	bool changed;

	//per base data
	TBase* bases;
	bool* aligned; //whether or not base is aligned to ref. Insertions are not aligned
	//soft clipped data
	int* softClippedLength;
	char** softClippedBase;
	char** softClippedQuality;
	int* qualityOriginal; //Note: quality is char as int: quality = (int) bam.quality
	int* qualityRecalibrated;
	int* quality; //pointer to qualities to be used
	char* baseAsChar; //TODO: to be removed, if possible



	//per base data
/*	Base* base;
	BaseContext* context;
	double* errorRates;

	int* distFrom3Prime;
	int* distFrom5Prime;
	double* pmdCT;
	double* pmdGA;

	//soft clipped data
	int* softClippedLength;
	char** softClippedBase;
	char** softClippedQuality;
*/
	uint8_t softClippedEntry; //0 means start, 1 means end of read


	//reference
	bool hasReference;
	std::string referenceSequence;

	//functions
	void initStorage();
	void freeStorage();

	//functions to read and parse
	void setDistancesFromEnds();
	void parseBasesQualities(TGenotypeMap & genoMap, TQualityMap & qualityMap);
	void fillContext(TGenotypeMap & genoMap);
	void fillPmdProbabilities(TPMD* pmdObjects);

	//functions to modify data
	void filterForPrintingBaseQuality(std::string & qual, int & minQualForPrinting, int & maxQualForPrinting);
	void trimRead(int & trimmingLength3Prime, int & trimmingLength5Prime);

	//functions that access data
	std::string name(){return alignmentName;};



public:
	TAlignment();
	TAlignment(unsigned int MaxSize);
	~TAlignment(){
		freeStorage();
	}

	void fill(BamTools::BamAlignment & bamAlignment, int ReadGroupId);
	void setFiltersPassed(bool passed);
	void setReferenceAdded();
	int getPosition(){return position;};

	//functions to write / print alignment
	void save(BamTools::BamWriter & bamWriter, TGenotypeMap & genoMap, int & minQualForPrinting, int & maxQualForPrinting);
	void print(TGenotypeMap & genoMap);

	//accessed by alignmentParser

	void filterForBaseQuality(int & minQual, int & maxQual);
	void clear();
	void parse(TGenotypeMap & genoMap, TQualityMap & qualityMap);

	//accessed by TGenome
	int readGroupId;
	bool isReverseStrand;
	bool isProperPair;
	int mappingQuality;
	bool passedFilters;
	int chrNumber;
	int32_t position;
	std::string alignmentName;

	//TODO: move these functions to TGenome
	void recalibrate(TRecalibration & recalObject, TQualityMap & qualityMap);
	void recalibrate(TRecalibration & recalObject, TPMD* pmdObjects, TFastaBuffer* fastaBuffer, TQualityMap & qualityMap);
	void fillReadGroupInfo(int & readGroupID);
	void binQualityScores(TQualityMap & qualityMap);
	void updateOptionalSamField(std::string tag, float value);
	void downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator);

	void addToPMDTables(TPMDTables & pmdTables, TGenotypeMap & genoMap);
	double calculatePMDS(double & pi, TPMD* pmdObjects);
	void assessSoftClipping(int & S_left, int & middle, int & S_right);
	void addToQualityTable(TQualityTable & qualTable);

	friend class TAlignmentParser;
	friend class TWindow;

};

#endif /* TALIGNMENT_H_ */
