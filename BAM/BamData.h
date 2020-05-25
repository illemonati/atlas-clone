/*
 * TSamFlags.h
 *
 *  Created on: May 20, 2020
 *      Author: phaentu
 */

#ifndef BAM_BAMDATA_H_
#define BAM_BAMDATA_H_

#include <bitset>
#include <vector>

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
	TSamFlags(){};
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

//-----------------------------------------------------
//TCigar
//A class to store, access and manipulate CIGAR operators
//-----------------------------------------------------
struct CigarOperator {
    char     type;
    uint32_t length;

    CigarOperator(){
    	type = '\0';
    	length = 0;
    };
    CigarOperator(const char & Type, const uint32_t & Length){
    	type = Type;
    	length = Length;
    };
};

class TCigar{
private:
	std::vector<CigarOperator> cigar;
	uint16_t _lengthAligned;
	uint16_t _lengthInserted;
	uint16_t _lengthDeleted;
	uint16_t _lengthSoftClipped;

public:
	TCigar(){
		clear();
	};

	void clear(){
		cigar.clear();
		_lengthAligned = 0;
		_lengthInserted = 0;
		_lengthDeleted = 0;
		_lengthSoftClipped = 0;
	};

	std::vector<CigarOperator>::const_iterator begin() const { return cigar.begin(); }
	std::vector<CigarOperator>::const_iterator end() const { return cigar.end(); }

    void add(const char & Type, const uint32_t & Length){
    	if(Type == 'M' || Type == '=' || Type == 'X'){
    		_lengthAligned += Length;
    	} else if(Type == 'I'){
    		_lengthInserted += Length;
    	} else if(Type =='D'){
    		_lengthDeleted += Length;
    	} else if(Type == 'S'){
    		_lengthSoftClipped += Length;
    	} else if(Type != 'N' && Type != 'H' && Type != 'P'){
    		throw "Unknown CIGAR operation '" + Type + "'!";
    	}

    	//add to vector
    	cigar.emplace_back(Type, Length);
    }

    uint16_t lengthAligned() const{ return _lengthAligned; };
    uint16_t lengthInserted() const{ return _lengthInserted; };
    uint16_t lengthDeleted() const{ return _lengthDeleted; };
    uint16_t lengthSoftClipped() const{ return lengthSoftClipped; };
    uint16_t lengthSequenced() const{ return _lengthAligned + _lengthInserted; };
};


}; //end namespace

#endif /* BAM_BAMDATA_H_ */
