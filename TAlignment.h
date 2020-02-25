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
#include "TQualityMap.h"
#include "QualityTables.h"
#include "TContextStats.h"
#include "TSoftClipping.h"

class TAlignmentParser;

//-----------------------------------------------------
//TAlignment
//-----------------------------------------------------

class TAlignment{
private:
	//details
	bool empty;
	bool recalibrated;

	//data
	unsigned int maxSize;
	uint16_t length;
	uint16_t fragmentLength; //is insert size - deletions + insertions for paired-end; read length - deletions + insertions for single end

	//BamAlignment data
	std::string alignmentName;
	int32_t insertSize;
	uint32_t position;

	//booleans
	bool parsed;
	bool changed;

	//per base data
	TBase* bases;
	int* softClippedLength;
	char** softClippedBase;
	char** softClippedQuality;
	int numInsertions;
	int numDeletions;
	uint8_t softClippedEntry; //0 means start, 1 means end of read

	//reference
	bool hasReference;
	std::string referenceSequence;

	//functions
	void initStorage();
	void freeStorage();

	//functions to read and parse
	void setDistancesFromEnds();
	void copyDataToBase(TBase & base, const char baseAsChar, const char qualAsChar, TGenotypeMap & genoMap, TQualityMap & qualityMap);
	void parseBasesQualities(TGenotypeMap & genoMap, TQualityMap & qualityMap);
	void fillContext(TGenotypeMap & genoMap);
	void fillPmdProbabilities(TPMD* pmdObjects);

	//functions to modify data
	void filterForPrintingBaseQuality(std::string & bases, std::string & qual, int minQualForPrinting, int maxQualForPrinting);
	void trimRead(int & trimmingLength3Prime, int & trimmingLength5Prime);
	void setReadTrimming(int trim3Prime, int trim5Prime);

public:
	bool storageInitialized;
	BamTools::BamAlignment bamAlignment;

	TAlignment();
	TAlignment(unsigned int MaxSize);
	TAlignment(const TAlignment & Alignment);

	~TAlignment(){
		freeStorage();
	};

	void fill(BamTools::BamAlignment & bamAlignment, const uint16_t ReadGroupId);
	void setReferenceAdded();
	uint32_t getPosition(){ return position; };
	uint16_t getParsedLength(){ return length; };
	uint32_t getBamAlignmentLength(){ return bamAlignment.Length; };
	std::string name(){ return alignmentName; };
	int32_t getInsertSize(){ return insertSize; };

	//functions to write / print alignment
	void setToSingleEnd();
	void setIsProperPair(const bool & ok);
	void save(BamTools::BamWriter & bamWriter, TGenotypeMap & genoMap, int minQualForPrinting, int maxQualForPrinting, TQualityMap & qualMap);
	void print(TGenotypeMap & genoMap, TQualityMap & qualMap);
	void setAlignmentHasChanged(){ changed = true; };

	//accessed by alignmentParser
	void filterForBaseQualityAsPhredInt(int & minQual, int & maxQual);
	void filterForContext(std::map<BaseContext,int> ignoreTheseContexts);
	void clear();
	void parse(TGenotypeMap & genoMap, TQualityMap & qualityMap);

	//accessed by TGenome
	uint16_t readGroupId;
	bool isReverseStrand;
	bool isPaired;
	bool isProperPair;
	uint16_t mappingQuality;
	bool isSecondMate;
	int32_t matePosition;
	uint32_t chrNumber;
	int32_t lastPositionPlusOne;
	int32_t lastAlignedPositionWithRespectToRef;
	int32_t lastAlignedPos;

	//TODO: move these functions to TGenome
	void fillReadGroupInfo(int & readGroupID);
	void binQualityScores(TQualityMap & qualityMap);
	void updateOptionalSamField(std::string tag, float value);
	void updateOptionalSamField(std::string tag, std::string value);
	void downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator, TQualityMap & qualMap);

	void addToPMDTables(TPMDTables & pmdTables, TGenotypeMap & genoMap);
	void recalibrateWithPMD(TRecalibration* recalObject, TQualityMap & qualMap);
	double calculatePMDS(double & pi, TPMD* pmdObjects);
	void removeSoftClippedBases(TSoftClippingData & softClippingData);
	int measureOverlap();
	int getUsableLength(const int minQual, const int maxQual);
	void addToQualityTable(TQualityTable & qualTable, TQualityMap & qualMap);
	void addToContextStats(TContextStats & contextStats);

	friend class TAlignmentParser;
	friend class TWindow;

};

#endif /* TALIGNMENT_H_ */
