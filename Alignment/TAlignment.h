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
#include "IOTools/IOAbstractClasses/BamAlignment.h"
#include "IOTools/IOAbstractClasses/Fasta.h"
#include "IOTools/IOAbstractClasses/BamWriter.h"
#include "IOTools/IOAbstractClasses/Cigar.h"
#include "IOTools/IOTool.h"
#include "TQualityMap.h"
#include "QualityTables.h"

class TAlignmentParser;
class TWindow;

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
	int length;

	//BamAlignment data
	std::string alignmentName;
	int32_t insertSize;
	int32_t position;

	//booleans
	bool parsed;
	bool changed;

	//per base data
	TBase* bases;
	int* softClippedLength;
	char** softClippedBase;
	char** softClippedQuality;
	int* qualityOriginal; //Note: quality is char as int: quality = (int) bam.quality
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
	void parseBasesQualities(TGenotypeMap & genoMap, TQualityMap & qualityMap);
	void fillContext(TGenotypeMap & genoMap);
	void fillPmdProbabilities(TPMD* pmdObjects);

	//functions to modify data
	void filterForPrintingBaseQuality(std::string & qual, int & minQualForPrinting, int & maxQualForPrinting);
	void trimRead(int & trimmingLength3Prime, int & trimmingLength5Prime);
	void setReadTrimming(int trim3Prime, int trim5Prime);

public:
	bool storageInitialized;
    BamAlignment* bamAlignment;

	TAlignment();
	TAlignment(unsigned int MaxSize);
	TAlignment(const TAlignment & Alignment);

	~TAlignment(){
		freeStorage();
	};

    void fill(BamAlignment & bamAlignment, int ReadGroupId);
	void setReferenceAdded();
    int getPosition(){ return position; }
    int getParsedLength(){ return length; }
    int getBamAlignmentLength(){ return bamAlignment->GetLength(); }
    std::string name(){ return alignmentName; }
    int32_t getInsertSize(){ return insertSize; }

	//functions to write / print alignment
	void setToSingleEnd();
	void setIsProperPair(const bool & ok);
    void save(BamWriter & bamWriter, TGenotypeMap & genoMap, int & minQualForPrinting, int & maxQualForPrinting, TQualityMap & qualMap);
	void print(TGenotypeMap & genoMap, TQualityMap & qualMap);
	void setAlignmentHasChanged(){ changed = true; };

	//accessed by alignmentParser
	void filterForBaseQuality(int & minQual, int & maxQual);
	void clear();
	void parse(TGenotypeMap & genoMap, TQualityMap & qualityMap);

	//accessed by TGenome
	int readGroupId;
	bool isReverseStrand;
	bool isPaired;
	bool isProperPair;
	int mappingQuality;
	bool passedFilters;
	bool isSecondMate;
	int32_t matePosition;
	int chrNumber;
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
	void assessSoftClipping(int & S_left, int & middle, int & S_right, std::string & S_string_left, std::string & S_string_middle, std::string & S_qualities_middle, std::string & S_string_right, TGenotypeMap & genoMap);
	void addToAlignment(BamAlignment* bamAlignment, int counter, std::string & queryBases, std::string & qualities, CigarOp& op);
	void removeSoftClippedBases(int & S_left, int & middle, int & S_right, std::string & S_string_left, std::string & S_string_middle, std::string & S_qualities_middle, std::string & S_string_right, TGenotypeMap & genoMap);
	int measureOverlap();
	void addToQualityTable(TQualityTable & qualTable, TQualityMap & qualMap);

	friend class TAlignmentParser;
	friend class TWindow;

};

#endif /* TALIGNMENT_H_ */
