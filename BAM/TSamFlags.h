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
	using bitset = coretools::TBitSet<12>;
	bitset _flags{0};
public:
	TSamFlags() noexcept = default;
	TSamFlags(const TSamFlags&) noexcept = default;
	TSamFlags& operator=(const TSamFlags&) noexcept = default;
	TSamFlags(uint16_t Flags) noexcept : _flags(Flags) {}
	TSamFlags& operator=(uint16_t Flags) noexcept { set(Flags); return *this; }

	// getters
	bool operator[](uint16_t index) noexcept { return _flags[index]; }
	uint16_t asInt() const noexcept { return _flags.to_ulong(); }
	void reset() noexcept { _flags.reset(); }
	bool isPaired() const noexcept { return _flags.get<0>(); }
	bool isProperPair() const noexcept { return _flags.get<1>(); }
	bool isUnmapped() const noexcept { return _flags.get<2>(); }
	bool mateIsUnmapped() const noexcept { return _flags.get<3>(); }

	bool isReverseStrand() const noexcept { return _flags.get<4>(); };
	bool mateIsReverseStrand() const noexcept { return _flags.get<5>(); };
	bool isFirstMate() const noexcept { return _flags.get<6>(); };  // first segment in template
	bool isSecondMate() const noexcept { return _flags.get<7>(); }; // last segment in template

	bool isSecondary() const noexcept { return _flags.get<8>(); };
	bool isQCFailed() const noexcept { return _flags.get<9>(); };
	bool isDuplicate() const noexcept { return _flags.get<10>(); };
	bool isSupplementary() const noexcept { return _flags.get<11>(); };

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

	// setters
	void set(uint16_t Flags) noexcept { _flags = bitset(Flags); };
	void setIsPaired(bool ok) noexcept { _flags.set<0>(ok); };
	void setIsProperPair(bool ok) noexcept { _flags.set<1>(ok); };
	void setIsUnmapped(bool ok) noexcept { _flags.set<2>(ok); };
	void setMateIsUnmapped(bool ok) noexcept { _flags.set<3>(ok); };

	void setIsReverseStrand(bool ok) noexcept  { _flags.set<4>(ok); };
	void setMateIsReverseStrand(bool ok) noexcept  { _flags.set<5>(ok); };
	void setIsSecondMate(bool ok) noexcept { _flags.set<7>(ok); _flags.set<6>(!ok);};

	void setIsSecondary(bool ok) noexcept  { _flags.set<8>(ok); };
	void setIsQCFailed(bool ok) noexcept  { _flags.set<9>(ok); };
	void setIsDuplicate(bool ok) noexcept  { _flags.set<10>(ok); };
	void setIsSupplementary(bool ok) noexcept  { _flags.set<11>(ok); };
};
static_assert(sizeof(TSamFlags) == 2);

}; // namespace BAM

#endif /* BAM_TSAMFLAGS_H_ */
