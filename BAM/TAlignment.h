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
#include "stringFunctions.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TQualityMap.h"
#include "TFastaBuffer.h"
#include "TRandomGenerator.h"
#include "TGenomePosition.h"
#include "TSequencedBase.h"
#include "GenotypeTypes.h"
#include <vector>

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
	std::vector<TSequencedBase> _bases;
	std::vector<int> _alignedPosition;

	//reference
	bool _hasReference;
	std::vector<genometools::Base> _referenceSequence;

	void _initialize();

	//functions to read and parse
	template <typename BASETYPE, typename QUALTYPE>
	void _parseBasesQualities(const BASETYPE & Sequence, const QUALTYPE & Qualities){
		//initialize
		_bases.resize(_cigar.lengthRead());
		_alignedPosition.resize(_cigar.lengthRead());
		int d = 0; //index regarding data structures and inside read
		int p = 0; //index regarding reference position (!= d for soft clipping & indels)

		//loop over cigar operations
		for(auto& cigarIter : _cigar){
			switch ( cigarIter.type ) {

				// for 'M', '=' or 'X': just copy
				case ('M') :
				case ('=') :
				case ('X') :
					//soft-clipped bases on 5' are before bamAlignment.Position
					for(unsigned int i=0; i<cigarIter.length; ++i, ++d, ++p){
						_bases[d].base = _sequence[d];
						_bases[d].originalQuality_phredInt = _qualities[d];
						_bases[d].setAligned(true);
						_alignedPosition[d] = p;
					}
					break;

				//for 'S' - soft clip: copy by set aligned = false
				case ('S') :
					//add bases to softclipped entries
					for(unsigned int i=0; i<cigarIter.length; ++i, ++d){
						//soft-clipped bases on 5' are before bamAlignment.Position
						//need to initialize quality for quality filter and bases for context
						_bases[d].base = _sequence[d];
						_bases[d].originalQuality_phredInt = _qualities[d];
						_bases[d].setAligned(false);
						_alignedPosition[d] = -1;
					}
					break;

				//for 'I' - insertion: copy bases, but put aligned  = false
				case ('I')      :
					for(unsigned int i=0; i<cigarIter.length; ++i, ++d){
						_bases[d].base = _sequence[d];
						_bases[d].originalQuality_phredInt = _qualities[d];
						_bases[d].setAligned(false);
						_alignedPosition[d] = -1;
					}
					break;

				// for 'D' - deletion: just add to position
				case ('D') :
					p += cigarIter.length;
					break;

				// for 'N' - skipped region in reference: only advance reference position
				case ('N') :
					p += cigarIter.length;
					break;

				// for 'H' or 'P' - hard clip: do nothing as these bases are not present in SEQ
				case ('H') :
				case ('P') :
					break;

				// invalid CIGAR op-code
				default:
					throw (std::string) "CIGAR operation '" + cigarIter.type + "' not supported!";
			}
		}

		//update length and last aligned position
		_lastAlignedPositionWithRespectToRef = *this + (p - 1);
		_lastAlignedPos = p - 1; //why -1? -> same reason as above

		//then update distances from ends
		_setDistancesFromEnds();

		//fill context for each base
		_fillContext();

		//set mapping quality and whether read is first or second
		for(auto& b : _bases){
			b.readGroupID = _readGroupID;
			b.mappingQuality = _mappingQuality;
			b.fragmentLength = _fragmentLength;
			b.setSecondMate(_flags.isSecondMate());
			b.setReverseStrand(_flags.isReverseStrand());
		}

		_parsed = true;
		_sequenceAndQualitiesChanged = false;
	};

	void _setDistancesFromEnds();
	void _fillContext();

	//functions to modify data
	void _updateSequenceAndQualities() const;

public:
	TAlignment(const uint32_t& RefID, const uint32_t& Position);
	TAlignment(const TGenomePosition & other);
	TAlignment();

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
	void parse();
	void parse(const GenotypeLikelihoods::TSequencingErrorModels & seqErrorModels);

	//setters
	void addReference(TFastaBuffer & fasta);
	void setHasChanged(){ _changed = true; };
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
	BAM::Base referenceAtInternalPos(const uint32_t & internalPosition) const;
	TGenomePosition positionInRef(const uint32_t & internalPosition) const;
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

	TSequencedBase& operator[](const uint32_t internalPos){ return _bases[internalPos]; };
	const TSequencedBase& operator[](const uint32_t internalPos) const { return _bases[internalPos]; };

	std::string sequence() const;
	std::string qualities() const;
	bool isEmpty() const { return _empty; };
	bool isParsed() const{ return _parsed; };
	bool isReverseStrand() const{ return _flags.isReverseStrand(); };
	bool isPaired() const{ return _flags.isPaired(); };
	bool isProperPair() const{ return _flags.isProperPair(); };

	//looping
	std::vector<TSequencedBase>::iterator begin(){ return _bases.begin(); };
	std::vector<TSequencedBase>::iterator end(){ return _bases.end(); };

	//filters and other functions to modify data
	void filter(const TBaseFilter & Filter);
	void trimRead(const int & trimmingLength3Prime, const int & trimmingLength5Prime);
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
