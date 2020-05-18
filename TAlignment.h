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
#include "bamtools/api/BamAlignment.h"
#include "bamtools/utils/bamtools_fasta.h"
#include "bamtools/api/BamWriter.h"

#include "TGenotypeLikelihoodCalculator.h"
#include "TQualityMap.h"
#include "QualityTables.h"
#include "TContextStats.h"
#include "TSoftClipping.h"
#include "TFastaBuffer.h"

class TAlignmentParser;



//-----------------------------------------------------
//TAlignment
//-----------------------------------------------------

class TAlignment{
private:
	//details
	bool empty;

	//data
	uint16_t maxSize;
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
	//ToDo: turn into std storage where possible
	TBase* bases;
	uint16_t* alignedPosition;
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
	void _initStorage();
	void _freeStorage();

	//functions to read and parse
	void _setDistancesFromEnds();
	void _parseBasesQualities(TGenotypeMap & genoMap, TQualityMap & qualityMap);
	void _fillContext(TGenotypeMap & genoMap);
	void _fillPmdProbabilities(GenotypeLikelihoods::TPMDDoubleStrand* pmdObjects);

	//functions to modify data
	void _trimRead(int & trimmingLength3Prime, int & trimmingLength5Prime);
	void _setReadTrimming(int trim3Prime, int trim5Prime);

public:
	bool storageInitialized;
	//BamTools::BamAlignment bamAlignment;

	TAlignment();
	TAlignment(uint16_t MaxSize);
	TAlignment(const TAlignment & Alignment);

	~TAlignment(){
		_freeStorage();
	};

	void fill(BamTools::BamAlignment & bamAlignment, const uint16_t ReadGroupId);
	uint32_t getPosition(){ return position; };
	uint16_t getParsedLength(){ return length; };
	uint32_t getBamAlignmentLength(){ return bamAlignment.Length; };
	std::string name(){ return alignmentName; };
	int32_t getInsertSize(){ return insertSize; };

	//functions to write / print alignment
	void setIsProperPair(const bool & ok);
	void save(BamTools::BamWriter & bamWriter, const TGenotypeMap & genoMap, const TQualityMap & qualMap);
	void print(TGenotypeMap & genoMap, TQualityMap & qualMap);
	void setAlignmentHasChanged(){ changed = true; };

	//accessed by alignmentParser
	void filterForBaseQualityAsPhredInt(int & minQual, int & maxQual);
	void filterForContext(std::map<BaseContext,int> ignoreTheseContexts);
	void clear();
	void parse(TGenotypeMap & genoMap, TQualityMap & qualityMap, GenotypeLikelihoods::TSequencingErrorModels & seqErrorModels);

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

	void fillReadGroupInfo(uint16_t & ReadGroupID);
	void addReference(TFastaBuffer & fasta);
	void binQualityScores(TQualityMap & qualityMap);
	void updateOptionalSamField(std::string tag, float value);
	void updateOptionalSamField(std::string tag, std::string value);
	void downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator, TQualityMap & qualMap);

	void addToPMDTables(GenotypeLikelihoods::TPMDTables & pmdTables, TGenotypeMap & genoMap);
	void addSitesToQualityTransformTable(TQualityTransformTables & QTtables);
	void addSitesToQualityTransformTable(GenotypeLikelihoods::TSequencingErrorModels & otherSeqErrors, TQualityTransformTables & QTtables);
	void recalibrateWithPMD(GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator);
	double calculatePMDS(double & pi, GenotypeLikelihoods::TPMDDoubleStrand* pmdObjects);
	void removeSoftClippedBases(TSoftClippingData & softClippingData);
	int measureOverlap();
	int getUsableLength(const int minQual, const int maxQual);
	void addToQualityTable(TQualityTable & qualTable, TQualityMap & qualMap);
	void addToContextStats(TContextStats & contextStats);

	friend class TAlignmentParser;
	friend class TWindow;

};

#endif /* TALIGNMENT_H_ */
