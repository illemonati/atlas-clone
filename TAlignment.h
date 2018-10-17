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
	unsigned int maxSize;
	int length;
	int lastAlignedPos;

//	int* alignedPos;

	bool parsed;
	bool changed;

	//per base data
	TBase* bases;
//	bool* aligned; //whether or not base is aligned to ref. Insertions and clipped bases are not aligned
	int* softClippedLength;
	char** softClippedBase;
	char** softClippedQuality;
	int* qualityOriginal; //Note: quality is char as int: quality = (int) bam.quality
	int numInsertions;
	int numDeletions;
//	int* qualityRecalibrated;
//	int* quality; //pointer to qualities to be used
//	char* baseAsChar; //TODO: to be removed, if possible
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
	void setReadTrimming(int trim3Prime, int trim5Prime);

	//functions that access data
	std::string name(){return alignmentName;};



public:

	bool storageInitialized;
	BamTools::BamAlignment bamAlignment;


	TAlignment();
	TAlignment(unsigned int MaxSize);
	TAlignment(TAlignment & Alignment);

	~TAlignment(){
		freeStorage();
	}

	void fill(BamTools::BamAlignment & bamAlignment, int ReadGroupId);
	void setReferenceAdded();
	int getPosition(){return position;};
	int getLength(){return length;};

	//functions to write / print alignment
	void setToSingleEnd();
	void save(BamTools::BamWriter & bamWriter, TGenotypeMap & genoMap, int & minQualForPrinting, int & maxQualForPrinting, TQualityMap & qualMap);
	void print(TGenotypeMap & genoMap, TQualityMap & qualMap);

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
	int32_t lastPositionPlusOne;
	std::string alignmentName;

	//TODO: move these functions to TGenome
	void fillReadGroupInfo(int & readGroupID);
	void binQualityScores(TQualityMap & qualityMap);
	void updateOptionalSamField(std::string tag, float value);
	void updateOptionalSamField(std::string tag, std::string value);
	void downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator, TQualityMap & qualMap);

	void addToPMDTables(TPMDTables & pmdTables, TGenotypeMap & genoMap);
	void recalibrateWithPMD(TRecalibration* recalObject, TQualityMap & qualMap);
	double calculatePMDS(double & pi, TPMD* pmdObjects);
	void assessSoftClipping(int & S_left, int & middle, int & S_right);
	int measureOverlap();
	void addToQualityTable(TQualityTable & qualTable, TQualityMap & qualMap);

	friend class TAlignmentParser;
	friend class TWindow;

};

#endif /* TALIGNMENT_H_ */
