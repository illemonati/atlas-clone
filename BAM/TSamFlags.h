/*
 * TSamFlags.h
 *
 *  Created on: Jun 9, 2020
 *      Author: wegmannd
 */

#ifndef BAM_TSAMFLAGS_H_
#define BAM_TSAMFLAGS_H_

#include <bitset>

namespace BAM{

//-----------------------------------------------------
//TSamFlags
//A class to store, access and manipulate SAM flags.
//SAM flags are bitwise flags and this class has the only purpose to provided named access
//-----------------------------------------------------

class TSamFlags{
private:
	std::bitset<16> flags;

public:
	TSamFlags(const uint16_t Flags){
		flags = std::bitset<16>(Flags);
	};

	//getters
	uint16_t asInt() const{ return flags.to_ulong(); };
	bool operator[](const uint16_t index){
		return flags[index];
	};
	bool isPaired() const{ return flags[0]; };
	bool isProperPair() const{ return flags[1]; };
	bool isUnmapped() const{ return flags[2]; };
	bool mateIsUnmapped() const{ return flags[3]; };

	bool isReverseStrand() const{ return flags[4]; };
	bool mateIsReverseStrand() const{ return flags[5]; };
	bool isFirstMate() const{ return flags[6]; }; //first segment in template
	bool isSecondMate() const{ return flags[7]; }; //last segment in template

	bool isSecondary() const{ return flags[8]; };
	bool isQCFailed() const{ return flags[9]; };
	bool isDuplicate() const{ return flags[10]; };
	bool isSupplementary() const{ return flags[11]; };

	//setters
	void set(const uint16_t & Flags){ flags = std::bitset<16>(Flags); };
	void operator=(const uint16_t & Flags){ flags = std::bitset<16>(Flags); };
	void operator=(const TSamFlags & other){ flags = other.flags; };
	void setIsPaired(const bool ok){flags[0] = ok; };
	void setIsProperPair(const bool ok){flags[1] = ok; };
	void setIsUnmapped(const bool ok){flags[2] = ok; };
	void setMateIsUnmapped(const bool ok){flags[3] = ok; };

	void setIsReverseStrand(const bool ok){flags[4] = ok; };
	void setMateIsReverseStrand(const bool ok){flags[5] = ok; };
	void setIsRead1(const bool ok){flags[6] = ok; };
	void setIsRead2(const bool ok){flags[7] = ok; };

	void setIsSecondary(const bool ok){flags[8] = ok; };
	void setIsQCFailed(const bool ok){flags[9] = ok; };
	void setIsDuplicate(const bool ok){flags[10] = ok; };
	void setIsSupplementary(const bool ok){flags[11] = ok; };
};

}; /7end namespace

#endif /* BAM_TSAMFLAGS_H_ */
