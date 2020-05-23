/*
 * TAlignment.h
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#ifndef TALIGNMENT_H_
#define TALIGNMENT_H_

#include <BAMData.h>
#include <TBase.h>
#include "stringFunctions.h"
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
	//Alignment data
	std::string name;
	BAM::TSamFlags flags;
	uint32_t refID;
	uint32_t position;
	uint16_t mappingQuality;
	BAM::TCigar cigar;
	uint32_t mateRefID;
	int32_t matePosition;
	int32_t insertSize_TLEN;
	std::string sequence;
	std::string qualities;
	uint16_t readGroupID;

	int32_t lastAlignedPositionWithRespectToRef;
	int32_t lastAlignedPos;


	//details


	//data
	uint16_t maxSize;
	uint16_t length;
	uint16_t fragmentLength; //is insert size - deletions + insertions for paired-end; read length - deletions + insertions for single end


	//booleans
	bool storageInitialized;
	bool empty;
	bool parsed;
	bool changed;
	bool sequenceAndQualitiesChanged;

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
	void _parseBasesQualities(const TGenotypeMap & genoMap, const TQualityMap & qualityMap);
	void _setDistancesFromEnds();
	void _fillContext();

	void _updateSequenceAndQualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap);

	//functions to modify data
	void _trimRead(int & trimmingLength3Prime, int & trimmingLength5Prime);
	void _setReadTrimming(int trim3Prime, int trim5Prime);

public:


	TAlignment();
	TAlignment(uint16_t MaxSize);
	TAlignment(const TAlignment & Alignment);

	~TAlignment(){
		_freeStorage();
	};

	//clear, fill and parse
	void clear();
	void fill(const	std::string Name,
			  const BAM::TSamFlags Flags,
			  const uint32_t RefID,
			  const uint32_t Position,
			  const uint16_t MappingQuality,
			  const BAM::TCigar Cigar,
			  const uint32_t MateRefID,
			  const int32_t MatePosition,
			  const int32_t InsertSize_TLEN,
			  const std::string Sequence,
			  const std::string Qualities,
			  const uint16_t ReadGroupId);
	void parse(const TGenotypeMap & genoMap, const TQualityMap & qualityMap, const GenotypeLikelihoods::TSequencingErrorModels & seqErrorModels);


	//getters
	std::string getName() const{ return name; };
	uint32_t getPosition() const{ return position; };
	uint16_t getParsedLength() const{ return length; };
	int32_t getInsertSize() const{ return insertSize_TLEN; };
	std::string getSequence(const TGenotypeMap & genoMap, const TQualityMap & qualMap);
	std::string getQualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap);

	//functions to write / print alignment
	void setIsProperPair(const bool & ok);
	void save(BamTools::BamWriter & bamWriter, const TGenotypeMap & genoMap, const TQualityMap & qualMap);
	void print(TGenotypeMap & genoMap, TQualityMap & qualMap);
	void setAlignmentHasChanged(){ changed = true; };

	//accessed by alignmentParser
	void filterForBaseQualityAsPhredInt(int & minQual, int & maxQual);
	void filterForContext(std::map<BaseContext,int> ignoreTheseContexts);




	void fillReadGroupInfo(uint16_t & ReadGroupID);
	void addReference(TFastaBuffer & fasta);
	void binQualityScores(TQualityMap & qualityMap);

	/*
	 * TODO: discuss if needed. Is only used by PMDS and maybe we just do not add that info to the BAM file?
	void updateOptionalSamField(std::string tag, float value);
	void updateOptionalSamField(std::string tag, std::string value);
	*/

	void downsampleAlignment(double& fraction, TRandomGenerator& randomGenerator, TQualityMap & qualMap);

	void addToPMDTables(GenotypeLikelihoods::TPMDTables & pmdTables, TGenotypeMap & genoMap);
	void addSitesToQualityTransformTable(TQualityTransformTables & QTtables);
	void addSitesToQualityTransformTable(GenotypeLikelihoods::TSequencingErrorModels & otherSeqErrors, TQualityTransformTables & QTtables);
	void recalibrateWithPMD(GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator);
	double calculatePMDS(double & pi, GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator);
	void removeSoftClippedBases(TSoftClippingData & softClippingData);
	int measureOverlap();
	int getUsableLength(const int minQual, const int maxQual);
	void addToQualityTable(TQualityTable & qualTable, TQualityMap & qualMap);
	void addToContextStats(TContextStats & contextStats);

	friend class TAlignmentParser;
	friend class TWindow;
	friend class TBamFile;

};

#endif /* TALIGNMENT_H_ */
