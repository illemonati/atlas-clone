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

#include "coretools/Types/probability.h"
#include "genometools/GenomePositions/TGenomePosition.h"

#include "TCigar.h"
#include "TSamFlags.h"
#include "TSequencedBase.h"
#include "genometools/Genotypes/Base.h"

namespace GenotypeLikelihoods {
class TErrorModels;
}
namespace GenotypeLikelihoods {
namespace SequencingError {
class TModels;
}
} // namespace GenotypeLikelihoods
namespace genometools {
class TFastaReader;
}

namespace BAM {

coretools::PhredInt makeIllumina(coretools::PhredInt value) noexcept;
char makeIllumina(char Quality) noexcept;

//-----------------------------------------------------
// TAlignment
//-----------------------------------------------------
class TAlignment final : public genometools::TGenomePosition {
private:
	// Alignment data
	std::string _name;
	TSamFlags _flags;
	coretools::PhredInt _mappingQuality;
	TCigar _cigar;
	TGenomePosition _mateGenomicPosition;
	int32_t _insertSize_TLEN = 0;
	size_t _readGroupID      = 0;
	size_t _bamID            = 0;
	uint16_t _fragmentLength = 0;

	size_t _refSize         = 0;
	int32_t _lastAlignedPos = 0;

	// booleans
	bool _empty   = true;
	bool _parsed  = false;

	// sequence and qualities. Mutable so that they can be recreated from bases even for const TAlignment
	mutable std::string _sequence;
	mutable std::string _qualities;
	mutable bool _sequenceAndQualitiesChanged = false;

	// per base data
	std::vector<TSequencedBase> _bases;
	std::vector<int> _alignedPosition;

	// reference
	std::vector<genometools::Base> _referenceSequence;

	// functions to read and parse
	void _setCigar(const TCigar &Cigar);
	void _parseBasesQualities();
	void _setQualitiesNoRecal();
	void _setDistancesFromEnds();
	void _fillContext();

	// functions to modify data
	void _updateSequenceAndQualities() const;

public:
	TAlignment(uint32_t RefID, uint32_t Position) : TGenomePosition(RefID, Position) {}
	TAlignment(const TGenomePosition &other) : TGenomePosition(other){};
	TAlignment() = default;

	// clear, fill and parse
	void clear();
	void fill(const std::string &Name, const TSamFlags &Flags, uint32_t RefID, uint32_t Position,
			  uint16_t MappingQuality, const TCigar &Cigar, uint32_t MateRefID, uint32_t MatePosition,
			  const int32_t &InsertSize_TLEN, const std::string &Sequence, const std::string &Qualities,
			  size_t BamID, size_t ReadGroupId);
	void parse();
	void parse(const GenotypeLikelihoods::SequencingError::TModels &seqErrorModels);

	// setters
	void addReference(const genometools::TFastaReader &fasta);
	void setSequenceAndQualitiesChanged() { _sequenceAndQualitiesChanged = true; };
	void setName(std::string Name) { _name = std::move(Name); };
	void setMappingQuality(coretools::PhredInt Mappingquality) { _mappingQuality = Mappingquality; };
	void setMateGenomicPosition(const TGenomePosition & Position) {
		_mateGenomicPosition.move(Position);
	};
	void setInsertSize(int32_t InsertSize) { _insertSize_TLEN = InsertSize; };
	void setSequenceQualities(const TCigar &Cigar, const std::vector<genometools::Base> &Sequence,
							  const std::vector<coretools::PhredInt> &Quals);
	void setReadGroup(uint16_t readGroupId);
	void setIsReverseStrand(bool IsReverse) { _flags.setIsReverseStrand(IsReverse); };
	void setIsRead1(bool IsRead1) { _flags.setIsRead1(IsRead1); };
	void setIsRead2(bool IsRead2) { _flags.setIsRead2(IsRead2); };
	void setSamFlags(BAM::TSamFlags Flags) { _flags = std::move(Flags); };
	void setCigarForUnitTest(const TCigar &Cigar) {_cigar = Cigar;}

	// getters: position
	uint32_t lastAlingedInternalPos() const { return _lastAlignedPos; };
	uint32_t getLastInternalPos() const;
	TGenomePosition lastAlignedPositionWithRespectToRef() const { return *this + (_refSize - 1); };
	bool isAlignedAtInternalPos(size_t internalPosition) const;
	genometools::Base referenceAtInternalPos(size_t internalPosition) const;
	TGenomePosition positionInRef(size_t internalPosition) const;
	const genometools::TGenomePosition &mateGenomicPosition() const { return _mateGenomicPosition; };
	uint32_t matePosition() const { return _mateGenomicPosition.position(); };
	uint32_t mateRefID() const { return _mateGenomicPosition.refID(); };
	uint16_t fragmentLength() const { return _fragmentLength; };

	const std::string& name() const { return _name; };
	uint16_t readGroupId() const { return _readGroupID; };
	size_t parsedLength() const { return _alignedPosition.size(); };
	uint32_t length() const { return _cigar.lengthRead(); };
	int32_t insertSize() const { return _insertSize_TLEN; };
	coretools::PhredInt mappingQuality() const { return _mappingQuality; };
	uint16_t flags() const { return _flags.asInt(); };
	const TCigar &cigar() const { return _cigar; };

	TSequencedBase &operator[](size_t internalPos) noexcept {
		assert(internalPos < _bases.size());
		return _bases[internalPos];
	};
	const TSequencedBase &operator[](size_t internalPos) const noexcept {
		assert(internalPos < _bases.size());
		return _bases[internalPos];
	};

	std::string sequence() const;
	std::string qualities() const;
	bool isEmpty() const noexcept { return _empty; };
	bool isParsed() const noexcept { return _parsed; };
	bool isReverseStrand() const noexcept { return _flags.isReverseStrand(); };
	bool isPaired() const noexcept { return _flags.isPaired(); };
	bool isProperPair() const noexcept { return _flags.isProperPair(); };
	bool isSecondMate() const noexcept { return _flags.isSecondMate(); };
	Mate mate() const noexcept {return static_cast<Mate>(_flags.isSecondMate());}

	// looping
	std::vector<TSequencedBase>::iterator begin() noexcept { return _bases.begin(); };
	std::vector<TSequencedBase>::iterator end() noexcept { return _bases.end(); };
	std::vector<TSequencedBase>::reverse_iterator rbegin() noexcept { return _bases.rbegin(); };
	std::vector<TSequencedBase>::reverse_iterator rend() noexcept { return _bases.rend(); };
	std::vector<TSequencedBase>::const_iterator begin() const noexcept { return _bases.begin(); };
	std::vector<TSequencedBase>::const_iterator end() const noexcept { return _bases.end(); };
	std::vector<TSequencedBase>::const_reverse_iterator rbegin() const noexcept { return _bases.rbegin(); };
	std::vector<TSequencedBase>::const_reverse_iterator rend() const noexcept { return _bases.rend(); };
	
	size_t size() const noexcept { return _bases.size(); }

	// filters and other functions to modify data
	template<typename Filter> void filter(const Filter &F) {
		// set quality = 0 and base = N if outside quality filter
		for (auto &b : _bases) {
			if (!F.pass(b)) {
				b.base         = genometools::Base::N;
				b.recalQuality = coretools::PhredInt::highest();
			}
		}
	}
	void trimRead(int trimmingLength3Prime, int trimmingLength5Prime);
	void trimSoftClips();
	void trimSoftClips(size_t maxNumberOfSoftClippedBases);
	void binQualityScoresIllumina();
	void recalibrateWithPMD(const GenotypeLikelihoods::TErrorModels &GLCalculator);
	void setIsProperPair(bool ok);
	void downsampleAlignment(coretools::Probability fraction);
	void merge(uint16_t overlapLength, size_t &mappedBasesClipped);
};

}; // namespace BAM

#endif /* TALIGNMENT_H_ */
