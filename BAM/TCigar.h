/*
 * TSamFlags.h
 *
 *  Created on: May 20, 2020
 *      Author: phaentu
 */

#ifndef BAM_TCIGAR_H_
#define BAM_TCIGAR_H_

#include <vector>
#include <stdint.h>

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
    CigarOperator(const char & Type, const uint32_t & Length){
    	type = Type;
    	length = Length;
    };
};

class TCigar{
private:
	std::vector<CigarOperator> _cigar;
	uint16_t _lengthAligned;
	uint16_t _lengthInserted;
	uint16_t _lengthDeleted;
	uint16_t _lengthSoftClippedLeft;
	uint16_t _lengthSoftClippedRight;
	bool _addSoftClippedLeft;

public:
	TCigar(){
		clear();
	};

	void clear();
	std::vector<CigarOperator>::const_iterator begin() const { return _cigar.begin(); }
	std::vector<CigarOperator>::const_iterator end() const { return _cigar.end(); }
    void add(const char & Type, const uint32_t & Length);
    void removeSoftClips();

    uint16_t lengthAligned() const{ return _lengthAligned; };
    uint16_t lengthInserted() const{ return _lengthInserted; };
    uint16_t lengthDeleted() const{ return _lengthDeleted; };
    uint16_t lengthSoftClippedLeft() const{ return _lengthSoftClippedLeft; };
    uint16_t lengthSoftClippedRight() const{ return _lengthSoftClippedRight; };
    uint16_t lengthSoftClipped() const{ return _lengthSoftClippedLeft + _lengthSoftClippedRight; };
    uint16_t lengthSequenced() const{ return _lengthAligned + _lengthInserted; };
    uint16_t lengthRead() const{ return _lengthAligned + _lengthInserted + _lengthSoftClippedLeft + _lengthSoftClippedRight; };
};

}; //end namespace

#endif /* BAM_TCIGAR_H_ */
