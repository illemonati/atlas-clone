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
	std::bitset<16> _flags;

public:
	TSamFlags(const uint16_t Flags){
		_flags = std::bitset<16>(Flags);
	};

	//getters
	uint16_t asInt() const{ return _flags.to_ulong(); };
	bool operator[](const uint16_t index){
		return _flags[index];
	};
	void reset(){ _flags.reset(); };
	bool isPaired() const{ return _flags[0]; };
	bool isProperPair() const{ return _flags[1]; };
	bool isUnmapped() const{ return _flags[2]; };
	bool mateIsUnmapped() const{ return _flags[3]; };

	bool isReverseStrand() const{ return _flags[4]; };
	bool mateIsReverseStrand() const{ return _flags[5]; };
	bool isFirstMate() const{ return _flags[6]; }; //first segment in template
	bool isSecondMate() const{ return _flags[7]; }; //last segment in template

	bool isSecondary() const{ return _flags[8]; };
	bool isQCFailed() const{ return _flags[9]; };
	bool isDuplicate() const{ return _flags[10]; };
	bool isSupplementary() const{ return _flags[11]; };

	//setters
	void set(const uint16_t & Flags){ _flags = std::bitset<16>(Flags); };
	void operator=(const uint16_t & Flags){ _flags = std::bitset<16>(Flags); };
	void operator=(const TSamFlags & other){ _flags = other._flags; };
	void setIsPaired(const bool ok){_flags[0] = ok; };
	void setIsProperPair(const bool ok){_flags[1] = ok; };
	void setIsUnmapped(const bool ok){_flags[2] = ok; };
	void setMateIsUnmapped(const bool ok){_flags[3] = ok; };

	void setIsReverseStrand(const bool ok){_flags[4] = ok; };
	void setMateIsReverseStrand(const bool ok){_flags[5] = ok; };
	void setIsRead1(const bool ok){_flags[6] = ok; };
	void setIsRead2(const bool ok){_flags[7] = ok; };

	void setIsSecondary(const bool ok){_flags[8] = ok; };
	void setIsQCFailed(const bool ok){_flags[9] = ok; };
	void setIsDuplicate(const bool ok){_flags[10] = ok; };
	void setIsSupplementary(const bool ok){_flags[11] = ok; };
};

}; //end namespace

#endif /* BAM_TSAMFLAGS_H_ */
