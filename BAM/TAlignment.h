/*
 * TAlignment.h
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#ifndef TALIGNMENT_H_
#define TALIGNMENT_H_

#include <BamData.h>
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

namespace BAM{

//-----------------------------------------------------
//TAlignment
//-----------------------------------------------------

class TAlignment{
private:
	//Alignment data
	std::string _name;
	BAM::TSamFlags _flags;
	uint32_t _refID;
	uint32_t _position;
	uint16_t _mappingQuality;
	BAM::TCigar _cigar;
	uint32_t _mateRefID;
	int32_t _matePosition;
	int32_t _insertSize_TLEN;
	std::string _sequence;
	std::string _qualities;
	uint16_t _readGroupID;
	uint16_t _fragmentLength;

	int32_t _lastAlignedPositionWithRespectToRef;
	int32_t _lastAlignedPos;
	//uint16_t maxSize;

	//booleans
	bool _empty;
	bool _parsed;
	bool _changed;
	bool _sequenceAndQualitiesChanged;

	//per base data
	std::vector<TBase> _bases;
	std::vector<uint16_t> _alignedPosition;

	//reference
	bool _hasReference;
	std::string _referenceSequence;

	//functions to read and parse
	void _parseBasesQualities(const TGenotypeMap & genoMap, const TQualityMap & qualityMap);
	void _setDistancesFromEnds();
	void _fillContext();

	//functions to modify data
	void _updateSequenceAndQualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap);

public:
	TAlignment();
	//TAlignment(const TAlignment & Alignment); //TODO: check that copy constructor works automatically. It should!

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
	void parse(const TGenotypeMap & genoMap, const TQualityMap & qualityMap);
	void parse(const TGenotypeMap & genoMap, const TQualityMap & qualityMap, const GenotypeLikelihoods::TSequencingErrorModels & seqErrorModels);
	void addReference(TFastaBuffer & fasta);
	void setHasChanged(){ _changed = true; };
	void setReadGroup(const uint16_t readGroupId);

	//getters
	std::string name() const{ return _name; };
	uint32_t position() const{ return _position; };
	uint32_t refID() const{ return _refID; };
	uint32_t lastAlingedInternalPos() const{ return _lastAlignedPos; };
	uint32_t lastAlignedPositionWithRespectToRef() const{ return _lastAlignedPositionWithRespectToRef; };
	bool isAlignedAtInternalPos(const uint32_t internalPosition) const;
	char referenceAtInternalPos(const uint32_t internalPosition) const;
	uint32_t positionInRef(const uint32_t internalPosition) const;
	uint16_t readGroupId() const { return _readGroupID; };
	uint16_t parsedLength() const{ return _cigar.lengthSequenced(); };
	int32_t insertSize() const{ return _insertSize_TLEN; };
	TBase& base(const uint32_t internalPos){ return _bases[internalPos]; };
	std::string sequence(const TGenotypeMap & genoMap, const TQualityMap & qualMap);
	std::string qualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap);
	bool isReverseStrand() const{ return _flags.isReverseStrand(); };
	bool isPaired() const{ return _flags.isPaired(); };
	bool isProperPair() const{ return _flags.isProperPair(); };
	bool isParsed() const{ return _parsed; };

	//looping
	std::vector<TBase>::iterator begin(){ return _bases.begin(); };
	std::vector<TBase>::iterator end(){ return _bases.end(); };

	//filters and other functions to modify data
	void filterForBaseQuality(TQualityFilter & qualFilter);
	void filterForContext(const std::map<BaseContext,int> & ignoreTheseContexts, const TGenotypeMap & genoMap);
	void trimRead(const int & trimmingLength3Prime, const int & trimmingLength5Prime);
	void removeSoftClippedBases();
	void binQualityScores(TQualityMap & qualityMap);
	void recalibrateWithPMD(GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator);
	void setIsProperPair(const bool & ok);
	void downsampleAlignment(const double fraction, TRandomGenerator& randomGenerator);

	//debug functions
	void print(TGenotypeMap & genoMap, TQualityMap & qualMap);
};

}; //end namespace

#endif /* TALIGNMENT_H_ */
