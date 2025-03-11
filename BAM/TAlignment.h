/*
 * TAlignment.h
 *
 *  Created on: Apr 19, 2018
 *      Author: wegmannd
 */

#ifndef TALIGNMENT_H_
#define TALIGNMENT_H_

#include <cstdint>
#include <stdint.h>
#include <string>
#include <vector>

#include "coretools/Types/probability.h"
#include "genometools/GenomePositions/TGenomePosition.h"

#include "TCigar.h"
#include "TSamFlags.h"
#include "TSequencedData.h"
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
class TAlignment {
private:
	genometools::TGenomePosition _position;
	// Alignment data
	std::string _name;
	TSamFlags _flags;
	coretools::PhredInt _mappingQuality;
	TCigar _cigar;
	genometools::TGenomePosition _mateGenomicPosition;
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
	std::vector<TSequencedData> _data;
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
	TAlignment(uint32_t RefID, uint32_t Position) : _position(RefID, Position) {}
	TAlignment(const genometools::TGenomePosition &other) : _position(other){}
	TAlignment() = default;
	constexpr size_t refID() const noexcept {return _position.refID();}
	constexpr size_t position() const noexcept {return _position.position();}
	constexpr void move(const genometools::TGenomePosition& GPos) noexcept {_position.move(GPos);}
	constexpr void advanceOnRef(size_t Len) noexcept { _position += Len; }

	constexpr genometools::TGenomePosition from() const noexcept {return _position;}
	constexpr genometools::TGenomePosition to() const { return _position + _refSize; }

	// clear, fill and parse
	void fill(const std::string &Name, const TSamFlags &Flags, uint32_t RefID, uint32_t Position,
			  uint16_t MappingQuality, const TCigar &Cigar, uint32_t MateRefID, uint32_t MatePosition,
			  const int32_t &InsertSize_TLEN, const std::string &Sequence, const std::string &Qualities,
			  size_t BamID, size_t ReadGroupId);
	void parse();
	void parse(const GenotypeLikelihoods::SequencingError::TModels &seqErrorModels);

	// setters
	void addReference(const genometools::TFastaReader &fasta);
	void setSequenceAndQualitiesChanged() { _sequenceAndQualitiesChanged = true; }
	void setName(std::string Name) { _name = std::move(Name); }
	void setMappingQuality(coretools::PhredInt Mappingquality) { _mappingQuality = Mappingquality; }
	void setMateGenomicPosition(const genometools::TGenomePosition &Position) { _mateGenomicPosition.move(Position); }
	void setInsertSize(int32_t InsertSize) { _insertSize_TLEN = InsertSize; }
	void setSequenceQualities(const TCigar &Cigar, const std::vector<genometools::Base> &Sequence,
							  const std::vector<coretools::PhredInt> &Quals);
	void setReadGroup(uint16_t readGroupId);
	void setIsPaired(bool ok) { _flags.setIsPaired(ok); }
	void setIsProperPair(bool ok) { _flags.setIsProperPair(ok); }
	void setIsReverseStrand(bool IsReverse) { _flags.setIsReverseStrand(IsReverse); }
	void setIsSecondMate(bool ok) { _flags.setIsSecondMate(ok); }
	void setSamFlags(BAM::TSamFlags Flags) { _flags = std::move(Flags); }
	void setCigarForUnitTest(const TCigar &Cigar) {_cigar = Cigar;}

	// getters: position
	uint32_t lastAlingedInternalPos() const { return _lastAlignedPos; }
	uint32_t getLastInternalPos() const;
	bool isAlignedAtInternalPos(size_t internalPosition) const;
	genometools::Base referenceAtInternalPos(size_t internalPosition) const;
	genometools::TGenomePosition positionInRef(size_t internalPosition) const;
	const genometools::TGenomePosition &mateGenomicPosition() const { return _mateGenomicPosition; }
	uint32_t matePosition() const { return _mateGenomicPosition.position(); }
	uint32_t mateRefID() const { return _mateGenomicPosition.refID(); }
	uint16_t fragmentLength() const { return _fragmentLength; }

	const std::string& name() const { return _name; }
	uint16_t readGroupId() const { return _readGroupID; }
	size_t parsedLength() const { return _alignedPosition.size(); }
	uint32_t length() const { return _cigar.lengthRead(); }
	int32_t insertSize() const { return _insertSize_TLEN; }
	coretools::PhredInt mappingQuality() const { return _mappingQuality; }
	uint16_t flags() const { return _flags.asInt(); }

	const TCigar &cigar() const { return _cigar; }
	void removeMappedRight(size_t Length) {_cigar.removeMappedRight(Length);}
	void removedMappedLeft(size_t Length) { advanceOnRef(_cigar.removeMappedLeft(Length)); }

	TSequencedData &operator[](size_t internalPos) noexcept {
		assert(internalPos < _data.size());
		return _data[internalPos];
	}
	const TSequencedData &operator[](size_t internalPos) const noexcept {
		assert(internalPos < _data.size());
		return _data[internalPos];
	}

	const std::string& sequence() const;
	const std::string& qualities() const;
	bool isEmpty() const noexcept { return _empty; }
	bool isParsed() const noexcept { return _parsed; }
	bool isReverseStrand() const noexcept { return _flags.isReverseStrand(); }
	bool isPaired() const noexcept { return _flags.isPaired(); }
	bool isProperPair() const noexcept { return _flags.isProperPair(); }
	bool isSecondMate() const noexcept { return _flags.isSecondMate(); }
	Mate mate() const noexcept {return static_cast<Mate>(_flags.isSecondMate());}
	Strand strand() const noexcept {return static_cast<Strand>(_flags.isReverseStrand());}

	// looping
	auto begin() noexcept { return _data.begin(); }
	auto end() noexcept { return _data.end(); }
	auto rbegin() noexcept { return _data.rbegin(); }
	auto rend() noexcept { return _data.rend(); }
	auto begin() const noexcept { return _data.begin(); }
	auto end() const noexcept { return _data.end(); }
	auto rbegin() const noexcept { return _data.rbegin(); }
	auto rend() const noexcept { return _data.rend(); }
	
	size_t size() const noexcept { assert(_data.empty() || _data.size() == length()); return length(); }

	// filters and other functions to modify data
	template<typename Filter> void filter(const Filter &F) {
		// set quality = 0 and base = N if outside quality filter
		for (auto &b : _data) {
			if (!F.pass(b)) {
				b.base         = genometools::Base::N;
				b.recalQuality = coretools::PhredInt::highest();
			}
		}
	}
	void trimRead(uint64_t trimmingLength3Prime, uint64_t trimmingLength5Prime);
	void trimSoftClips();
	void trimSoftClips(size_t maxNumberOfSoftClippedBases);
	void binQualityScoresIllumina();
	void recalibrateWithPMD(const GenotypeLikelihoods::TErrorModels &GLCalculator);
	void downsampleAlignment(coretools::Probability fraction);

	friend constexpr bool operator<(const TAlignment& lhs, const TAlignment& rhs) {
		return lhs.from() < rhs.from();
	}
};

} // namespace BAM

#endif /* TALIGNMENT_H_ */
