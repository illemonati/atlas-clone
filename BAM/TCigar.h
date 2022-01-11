/*
 * TSamFlags.h
 *
 *  Created on: May 20, 2020
 *      Author: phaentu
 */

#ifndef BAM_TCIGAR_H_
#define BAM_TCIGAR_H_

#include <vector>
#include <cstdint>

namespace BAM{

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
    CigarOperator(char Type, uint32_t Length){
    	type = Type;
    	length = Length;
    };
};

class TCigar{
private:
	std::vector<CigarOperator> _cigar;
	uint32_t _lengthAligned;
	uint32_t _lengthInserted;
	uint32_t _lengthDeleted;
	uint32_t _lengthSkipped;
	uint32_t _lengthSoftClippedLeft;
	uint32_t _lengthSoftClippedRight;
	bool _addSoftClippedLeft;
	bool _softClippedAddedOnRight;

public:
	TCigar(){
		clear();
	};

	void clear();
	std::vector<CigarOperator>::const_iterator begin() const { return _cigar.begin(); }
	std::vector<CigarOperator>::const_iterator end() const { return _cigar.end(); }
    void add(char Type, uint32_t Length);
    void removeSoftClips();

    uint32_t lengthAligned() const { return _lengthAligned; };
    uint32_t lengthInserted() const { return _lengthInserted; };
    uint32_t lengthDeleted() const { return _lengthDeleted; };
    uint32_t lengthSoftClippedLeft() const { return _lengthSoftClippedLeft; };
    uint32_t lengthSoftClippedRight() const { return _lengthSoftClippedRight; };
    uint32_t lengthSoftClipped() const { return _lengthSoftClippedLeft + _lengthSoftClippedRight; };
    uint32_t lengthSequenced() const { return _lengthAligned + _lengthInserted; };
    uint32_t lengthRead() const { return _lengthAligned + _lengthInserted + _lengthSoftClippedLeft + _lengthSoftClippedRight; };
    uint32_t lengthMapped() const { return _lengthAligned + _lengthDeleted + _lengthSkipped; };
};

}; //end namespace

#endif /* BAM_TCIGAR_H_ */
