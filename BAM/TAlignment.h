/*
 * TAlignment.h
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#ifndef TALIGNMENT_H_
#define TALIGNMENT_H_

#include <stdint.h>
#include <string>
#include <vector>

#include "GenotypeTypes.h"
#include "TCigar.h"
#include "TGenomePosition.h"
#include "TSamFlags.h"
#include "TSequencedBase.h"
#include "probability.h"
#include "devtools.h"

namespace BAM { class TBaseFilter; }
namespace BAM { class TFastaBuffer; }
namespace GenotypeLikelihoods { class TGenotypeLikelihoodCalculator; }
namespace GenotypeLikelihoods { namespace SequencingError { class TModels; } }
namespace coretools { class TRandomGenerator; }
namespace genometools { class PhredIntProbability; }

namespace BAM{

//-----------------------------------------------------
//TAlignment
//-----------------------------------------------------
class TAlignment:public genometools::TGenomePosition{
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
	std::vector<TSequencedBase> _bases;
	std::vector<int> _alignedPosition;

	//reference
	bool _hasReference;
	std::vector<genometools::Base> _referenceSequence;

	void _initialize();

	//functions to read and parse
	void _parseBasesQualities();
	void _parseBasesQualities(const std::vector<genometools::Base> & Sequence, const std::vector<genometools::PhredIntProbability> & Qualities);
	void _setQualitiesNoRecal();
	void _setDistancesFromEnds();
	void _fillContext();

	//functions to modify data
	void _updateSequenceAndQualities() const;

public:
	TAlignment(uint32_t RefID, uint32_t Position);
	TAlignment(const TGenomePosition & other);
	TAlignment();
	~TAlignment() = default;
	TAlignment(const TAlignment&) = default;
	TAlignment(TAlignment&&) = default;
	TAlignment& operator=(const TAlignment&) = default;
	TAlignment& operator=(TAlignment&&) = default;

	//clear, fill and parse
	void clear();
	void fill(const	std::string & Name,
			  const TSamFlags & Flags,
			  uint32_t RefID,
			  uint32_t Position,
			  uint16_t MappingQuality,
			  const TCigar & Cigar,
			  uint32_t MateRefID,
			  uint32_t MatePosition,
			  const int32_t & InsertSize_TLEN,
			  const std::string & Sequence,
			  const std::string & Qualities,
			  uint16_t ReadGroupId);
	void parse();
	void parse(const GenotypeLikelihoods::SequencingError::TModels & seqErrorModels);

	//setters
	void addReference(TFastaBuffer & fasta);
	void setSequenceAndQualitiesChanged(){ _sequenceAndQualitiesChanged = true; };
	void setName(const std::string Name){ _name = Name; };
	void setMappingQuality(const uint16_t Mappingquality){ _mappingQuality = Mappingquality; };
	void setMateGenomicPosition(const uint32_t RefID, const uint32_t Position){ _mateGenomicPosition.move(RefID, Position); };
	void setInsertSize(const int32_t InsertSize){ _insertSize_TLEN = InsertSize; };
	void setSequenceQualities(const TCigar & Cigar, const std::vector<genometools::Base> & Sequence, const std::vector<genometools::PhredIntProbability> & Quals);
	void setReadGroup(const uint16_t readGroupId);
	void setIsReverseStrand(const bool IsReverse){ _flags.setIsReverseStrand(IsReverse); };
    void setIsRead1(const bool IsRead1){ _flags.setIsRead1(IsRead1); };
    void setIsRead2(const bool IsRead2){ _flags.setIsRead2(IsRead2); };
    void setSamFlags(const BAM::TSamFlags Flags){_flags = Flags; };

	//getters: position
	uint32_t lastAlingedInternalPos() const{ return _lastAlignedPos; };
	TGenomePosition lastAlignedPositionWithRespectToRef() const{ return _lastAlignedPositionWithRespectToRef; };
	bool isAlignedAtInternalPos(const uint32_t internalPosition) const;
	genometools::Base referenceAtInternalPos(uint32_t internalPosition) const;
	TGenomePosition positionInRef(uint32_t internalPosition) const;
	const genometools::TGenomePosition& mateGenomicPosition() const{ return _mateGenomicPosition; };
	uint32_t matePosition() const{ return _mateGenomicPosition.position(); };
	uint32_t mateRefID() const{ return _mateGenomicPosition.refID(); };

	std::string name() const { return _name; };
	uint16_t readGroupId() const { return _readGroupID; };
	uint32_t parsedLength() const { return _parsed ? length() : 0; };
	uint32_t length() const { return _cigar.lengthRead(); };
	int32_t insertSize() const { return _insertSize_TLEN; };
	uint16_t mappingQuality() const { return _mappingQuality; };
	uint16_t flags() const { return _flags.asInt(); };
	const TCigar &cigar() const { return _cigar; };

	TSequencedBase& operator[](const uint32_t internalPos){ return _bases[internalPos]; };
	const TSequencedBase& operator[](const uint32_t internalPos) const { return _bases[internalPos]; };

	std::string sequence() const;
	std::string qualities() const;
	bool isEmpty() const { return _empty; };
	bool isParsed() const { return _parsed; };
	bool isReverseStrand() const { return _flags.isReverseStrand(); };
	bool isPaired() const { return _flags.isPaired(); };
	bool isProperPair() const { return _flags.isProperPair(); };
	bool isSecondMate() const { return _flags.isSecondMate(); };

	//looping
	std::vector<TSequencedBase>::iterator begin(){ return _bases.begin(); };
	std::vector<TSequencedBase>::iterator end(){ return _bases.end(); };

	//filters and other functions to modify data
	void filter(const TBaseFilter & Filter);
	void trimRead(int trimmingLength3Prime, int trimmingLength5Prime);
	void removeSoftClippedBases();
	void binQualityScoresIllumina();
	void recalibrateWithPMD(const GenotypeLikelihoods::TGenotypeLikelihoodCalculator & GLCalculator);
	void setIsProperPair(const bool & ok);
	void downsampleAlignment(const coretools::Probability & fraction, coretools::TRandomGenerator& randomGenerator);

	//debug functions
	void print();
};

}; //end namespace

#endif /* TALIGNMENT_H_ */
