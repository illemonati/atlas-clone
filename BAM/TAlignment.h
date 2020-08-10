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
class TAlignment:public TGenomePosition{
private:
	//Alignment data
	std::string _name;
	TSamFlags _flags;
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
	void _parseBasesQualities(const GenotypeLikelihoods::TGenotypeMap & genoMap, const TQualityMap & qualityMap);
	void _setDistancesFromEnds();
	void _fillContext();

	//functions to modify data
	void _updateSequenceAndQualities(const GenotypeLikelihoods::TGenotypeMap & genoMap, const TQualityMap & qualMap) const;

public:
	TAlignment();
	//TAlignment(const TAlignment & Alignment); //TODO: check that copy constructor works automatically. It should!

	//clear, fill and parse
	void clear();
	void fill(const	std::string & Name,
			  const TSamFlags & Flags,
			  const uint32_t & RefID,
			  const uint32_t & Position,
			  const uint16_t & MappingQuality,
			  const TCigar & Cigar,
			  const uint32_t & MateRefID,
			  const uint32_t & MatePosition,
			  const int32_t & InsertSize_TLEN,
			  const std::string & Sequence,
			  const std::string & Qualities,
			  const uint16_t & ReadGroupId);
	void parse(const GenotypeLikelihoods::TGenotypeMap & genoMap, const TQualityMap & qualityMap);
	void parse(const GenotypeLikelihoods::TGenotypeMap & genoMap, const TQualityMap & qualityMap, const GenotypeLikelihoods::TSequencingErrorModels & seqErrorModels);

	//setters
	void addReference(TFastaBuffer & fasta);
	void setHasChanged(){ _changed = true; };
	void setName(const std::string Name){ _name = Name; };
	void setMappingQuality(const uint16_t Mappingquality){ _mappingQuality = Mappingquality; };
	void setMateGenomicPosition(const uint32_t RefID, const uint32_t Position){ _mateGenomicPosition.move(RefID, Position); };
	void setInsertSize(const int32_t InsertSize){ _insertSize_TLEN = InsertSize; };
	void setSequenceQualities(const TCigar Cigar, const std::string Sequence, const std::string Qualities);
	void setSequenceQualitiesOnlyMatches(const std::string Sequence, const std::string Qualities);
	void setReadGroup(const uint16_t readGroupId);
	void setIsReverseStrand(const bool IsReverse){ _flags.setIsReverseStrand(IsReverse); };

	//getters: position
	uint32_t lastAlingedInternalPos() const{ return _lastAlignedPos; };
	TGenomePosition lastAlignedPositionWithRespectToRef() const{ return _lastAlignedPositionWithRespectToRef; };
	bool isAlignedAtInternalPos(const uint32_t internalPosition) const;
	char referenceAtInternalPos(const uint32_t internalPosition) const;
	TGenomePosition positionInRef(const uint32_t internalPosition) const;
	const BAM::TGenomePosition& mateGenomicPosition() const{ return _mateGenomicPosition; };
	uint32_t matePosition() const{ return _mateGenomicPosition.position(); };
	uint32_t mateRefID() const{ return _mateGenomicPosition.refID(); };

	std::string name() const{ return _name; };
	uint16_t readGroupId() const { return _readGroupID; };
	uint16_t parsedLength() const;
	int32_t insertSize() const{ return _insertSize_TLEN; };
	uint16_t mappingQuality() const{ return _mappingQuality; };
	uint16_t flags() const{ return _flags.asInt(); };
	const TCigar& cigar() const{ return _cigar; };

	TBase& operator[](const uint32_t internalPos){ return _bases[internalPos]; };
	const TBase& at(const uint32_t internalPos) const { return _bases[internalPos]; };

	std::string sequence(const GenotypeLikelihoods::TGenotypeMap & genoMap, const TQualityMap & qualMap) const;
	std::string qualities(const GenotypeLikelihoods::TGenotypeMap & genoMap, const TQualityMap & qualMap) const;
	bool isReverseStrand() const{ return _flags.isReverseStrand(); };
	bool isPaired() const{ return _flags.isPaired(); };
	bool isProperPair() const{ return _flags.isProperPair(); };
	bool isParsed() const{ return _parsed; };

	//looping
	std::vector<TBase>::iterator begin(){ return _bases.begin(); };
	std::vector<TBase>::iterator end(){ return _bases.end(); };

	//filters and other functions to modify data
	void filterForBaseQuality(TQualityFilter & qualFilter);
	void filterForContext(const std::map<GenotypeLikelihoods::BaseContext,int> & ignoreTheseContexts, const GenotypeLikelihoods::TGenotypeMap & genoMap);
	void trimRead(const int & trimmingLength3Prime, const int & trimmingLength5Prime);
	void removeSoftClippedBases();
	void binQualityScores(TQualityMap & qualityMap);
	void recalibrateWithPMD(const GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator);
	void setIsProperPair(const bool & ok);
	void downsampleAlignment(const double fraction, TRandomGenerator& randomGenerator);

	//debug functions
	void print(const GenotypeLikelihoods::TGenotypeMap & genoMap, const TQualityMap & qualMap);
};

}; //end namespace

#endif /* TALIGNMENT_H_ */
