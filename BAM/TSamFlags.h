/*
 * TSamFlags.h
 *
 *  Created on: Jun 9, 2020
 *      Author: wegmannd
 */

#ifndef BAM_TSAMFLAGS_H_
#define BAM_TSAMFLAGS_H_

#include "coretools/Containers/TBitSet.h"

namespace BAM {

//-----------------------------------------------------
// TSamFlags
// A class to store, access and manipulate SAM flags.
// SAM flags are bitwise flags and this class has the only purpose to provided named access
//-----------------------------------------------------

class TSamFlags {
private:
	enum class ReadProperty : size_t {
		min,
		Paired = min,
		ProperPair,
		Unmapped,
		MateUnmapped,
		ReverseStrand,
		MateReverseStrand,
		FirstMate,
		SecondMate,
		Secondary,
		QCFailed,
		Duplicate,
		Supplementary,
		max
	};
	using bitset = coretools::TStrongBitSet<ReadProperty>;
	bitset _flags{0};

public:
	TSamFlags() noexcept = default;
	TSamFlags(const TSamFlags&) noexcept = default;
	TSamFlags& operator=(const TSamFlags&) noexcept = default;
	TSamFlags(uint16_t Flags) noexcept : _flags(Flags) {}
	TSamFlags& operator=(uint16_t Flags) noexcept { set(Flags); return *this; }

	// getters
	uint16_t asInt() const noexcept { return _flags.to_ulong(); }
	void reset() noexcept { _flags.reset(); }
	bool isPaired() const noexcept { return _flags.get<ReadProperty::Paired>(); }
	bool isProperPair() const noexcept { return _flags.get<ReadProperty::ProperPair>(); }
	bool isUnmapped() const noexcept { return _flags.get<ReadProperty::Unmapped>(); }
	bool mateIsUnmapped() const noexcept { return _flags.get<ReadProperty::MateUnmapped>(); }

	bool isReverseStrand() const noexcept { return _flags.get<ReadProperty::ReverseStrand>(); };
	bool mateIsReverseStrand() const noexcept { return _flags.get<ReadProperty::MateReverseStrand>(); };
	bool isFirstMate() const noexcept { return _flags.get<ReadProperty::FirstMate>(); };  // first segment in template
	bool isSecondMate() const noexcept { return _flags.get<ReadProperty::SecondMate>(); }; // last segment in template

	bool isSecondary() const noexcept { return _flags.get<ReadProperty::Secondary>(); };
	bool isQCFailed() const noexcept { return _flags.get<ReadProperty::QCFailed>(); };
	bool isDuplicate() const noexcept { return _flags.get<ReadProperty::Duplicate>(); };
	bool isSupplementary() const noexcept { return _flags.get<ReadProperty::Supplementary>(); };

	// valid combination?
	bool isValid() const noexcept {
		if (!isPaired() &&
		    (isProperPair() || mateIsUnmapped() || mateIsReverseStrand() || isFirstMate() || isSecondMate())) {
			return false;
		}
		if (isUnmapped() &&
		    (isProperPair() || isSecondary() || isSupplementary())) {
			return false;
		}
		return true;
	}

	void makeSingle() noexcept {
		_flags.set<ReadProperty::Paired>(false);
		_flags.set<ReadProperty::ProperPair>(false);
		_flags.set<ReadProperty::MateUnmapped>(false);
		_flags.set<ReadProperty::MateReverseStrand>(false);
		_flags.set<ReadProperty::FirstMate>(false);
		_flags.set<ReadProperty::SecondMate>(false);
	}

	// setters
	void set(uint16_t Flags) noexcept { _flags = bitset(Flags); };
	void set(ReadProperty RP, bool ok) noexcept { _flags[RP] = ok;}

	void setIsPaired(bool ok) noexcept { _flags.set<ReadProperty::Paired>(ok); };
	void setIsProperPair(bool ok) noexcept { _flags.set<ReadProperty::ProperPair>(ok); };
	void setIsUnmapped(bool ok) noexcept { _flags.set<ReadProperty::Unmapped>(ok); };
	void setMateIsUnmapped(bool ok) noexcept { _flags.set<ReadProperty::MateUnmapped>(ok); };

	void setIsReverseStrand(bool ok) noexcept  { _flags.set<ReadProperty::ReverseStrand>(ok); };
	void setMateIsReverseStrand(bool ok) noexcept  { _flags.set<ReadProperty::MateReverseStrand>(ok); };
	void setIsSecondMate(bool ok) noexcept { _flags.set<ReadProperty::SecondMate>(ok); _flags.set<ReadProperty::FirstMate>(!ok);};

	void setIsSecondary(bool ok) noexcept  { _flags.set<ReadProperty::Secondary>(ok); };
	void setIsQCFailed(bool ok) noexcept  { _flags.set<ReadProperty::QCFailed>(ok); };
	void setIsDuplicate(bool ok) noexcept  { _flags.set<ReadProperty::Duplicate>(ok); };
	void setIsSupplementary(bool ok) noexcept  { _flags.set<ReadProperty::Supplementary>(ok); };
};
static_assert(sizeof(TSamFlags) == 2);

}; // namespace BAM

#endif /* BAM_TSAMFLAGS_H_ */
