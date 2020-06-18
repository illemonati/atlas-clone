/*
 * TAlignment.h
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#ifndef TALIGNMENT_H_
#define TALIGNMENT_H_

#include "TCigar.h"
#include "TSamFlags.h"
#include "TBase.h"
#include "stringFunctions.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TQualityMap.h"
#include "TFastaBuffer.h"
#include "TRandomGenerator.h"
#include "TGenomePosition.h"

namespace BAM{

//-----------------------------------------------------
//TAlignment
//-----------------------------------------------------
class TAlignment{
private:
	//Alignment data
	std::string _name;
	TSamFlags _flags;
	TGenomePosition _genomicPosition;
	uint16_t _mappingQuality;
	TCigar _cigar;
	TGenomePosition _mateGenomicPosition;
	int32_t _insertSize_TLEN;
	uint16_t _readGroupID;
	uint16_t _fragmentLength;

	TGenomePosition _lastAlignedPositionWithRespectToRef;
	int32_t _lastAlignedPos;

	//booleans
	bool _empty;
	bool _parsed;
	bool _changed;

	//sequence and qualities. Mutable so that they can be recreated from bases even for const TAlignment
	mutable std::string _sequence;
	mutable std::string _qualities;
	mutable bool _sequenceAndQualitiesChanged;

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
	void _updateSequenceAndQualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap) const;

public:
	TAlignment();
	//TAlignment(const TAlignment & Alignment); //TODO: check that copy constructor works automatically. It should!

	//clear, fill and parse
	void clear();
	void fill(const	std::string Name,
			  const TSamFlags Flags,
			  const uint32_t RefID,
			  const uint32_t Position,
			  const uint16_t MappingQuality,
			  const TCigar Cigar,
			  const uint32_t MateRefID,
			  const uint32_t MatePosition,
			  const int32_t InsertSize_TLEN,
			  const std::string Sequence,
			  const std::string Qualities,
			  const uint16_t ReadGroupId);
	void parse(const TGenotypeMap & genoMap, const TQualityMap & qualityMap);
	void parse(const TGenotypeMap & genoMap, const TQualityMap & qualityMap, const GenotypeLikelihoods::TSequencingErrorModels & seqErrorModels);

	//setters
	void addReference(TFastaBuffer & fasta);
	void setHasChanged(){ _changed = true; };
	void setName(const std::string Name){ _name = Name; };
	void setGenomicPosition(const uint32_t RefID, const uint32_t Position){ _genomicPosition.update(RefID, Position); };
	void setMappingQuality(const uint16_t Mappingquality){ _mappingQuality = Mappingquality; };
	void setMateGenomicPosition(const uint32_t RefID, const uint32_t Position){ _mateGenomicPosition.update(RefID, Position); };
	void setInsertSize(const int32_t InsertSize){ _insertSize_TLEN = InsertSize; };
	void setSequenceQualities(const TCigar Cigar, const std::string Sequence, const std::string Qualities);
	void setReadGroup(const uint16_t readGroupId);
	void setIsReverseStrand(const bool IsReverse){ _flags.setIsReverseStrand(IsReverse); };

	//getters
	std::string name() const{ return _name; };
	uint32_t position() const{ return _genomicPosition.position(); };
	uint32_t refID() const{ return _genomicPosition.refID(); };
	const BAM::TGenomePosition& genomicPosition() const{ return _genomicPosition; };
	uint32_t matePosition() const{ return _mateGenomicPosition.position(); };
	uint32_t mateRefID() const{ return _mateGenomicPosition.refID(); };
	const BAM::TGenomePosition& mateGenomicPosition() const{ return _mateGenomicPosition; };
	uint32_t lastAlingedInternalPos() const{ return _lastAlignedPos; };
	TGenomePosition lastAlignedPositionWithRespectToRef() const{ return _lastAlignedPositionWithRespectToRef; };
	bool isAlignedAtInternalPos(const uint32_t internalPosition) const;
	char referenceAtInternalPos(const uint32_t internalPosition) const;
	uint32_t positionInRef(const uint32_t internalPosition) const;
	uint16_t readGroupId() const { return _readGroupID; };
	uint16_t parsedLength() const{ return _cigar.lengthSequenced(); };
	int32_t insertSize() const{ return _insertSize_TLEN; };
	TBase& base(const uint32_t internalPos){ return _bases[internalPos]; };
	std::string sequence(const TGenotypeMap & genoMap, const TQualityMap & qualMap) const;
	std::string qualities(const TGenotypeMap & genoMap, const TQualityMap & qualMap) const;
	bool isReverseStrand() const{ return _flags.isReverseStrand(); };
	bool isPaired() const{ return _flags.isPaired(); };
	bool isProperPair() const{ return _flags.isProperPair(); };
	bool isParsed() const{ return _parsed; };

	//looping
	std::vector<TBase>::iterator begin(){ return _bases.begin(); };
	std::vector<TBase>::iterator end(){ return _bases.end(); };

	bool operator<(const TAlignment & other) const{
		return this->_genomicPosition < other._genomicPosition;
	};
	bool operator<(const BAM::TGenomePosition & other) const{
		return this->_genomicPosition < other;
	};

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
