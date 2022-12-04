/*
 * TSamFlags.h
 *
 *  Created on: May 20, 2020
 *      Author: phaentu
 */

#ifndef BAM_TCIGAR_H_
#define BAM_TCIGAR_H_

#include <cstdint>
#include <vector>
#include <string>

namespace BAM {

//-----------------------------------------------------
// TCigar
// A class to store, access and manipulate CIGAR operators
//-----------------------------------------------------
struct CigarOperator {
	char type       = '\0';
	uint32_t length = 0;

	CigarOperator(char Type, uint32_t Length) : type(Type), length(Length) {}
};

class TCigar {
private:
	std::vector<CigarOperator> _cigar;
	uint32_t _lengthAligned          = 0;
	uint32_t _lengthInserted         = 0;
	uint32_t _lengthDeleted          = 0;
	uint32_t _lengthSkipped          = 0;
	uint32_t _lengthSoftClippedLeft  = 0;
	uint32_t _lengthSoftClippedRight = 0;
	void _flipCigar();

public:
	TCigar()=default;
	TCigar(TCigar cigar, uint16_t overlapLength, bool isFirst, size_t &mappedBasesClipped);
	void clear();
	std::vector<CigarOperator>::const_iterator begin() const noexcept { return _cigar.begin(); }
	std::vector<CigarOperator>::const_iterator end() const noexcept { return _cigar.end(); }
	std::vector<CigarOperator>::const_reverse_iterator rbegin() const noexcept { return _cigar.rbegin(); }
	std::vector<CigarOperator>::const_reverse_iterator rend() const noexcept { return _cigar.rend(); }
	void add(char Type, uint32_t Length);
	void removeSoftClips();

	constexpr uint32_t lengthAligned() const noexcept { return _lengthAligned; };
	constexpr uint32_t lengthInserted() const noexcept { return _lengthInserted; };
	constexpr uint32_t lengthDeleted() const noexcept { return _lengthDeleted; };
	constexpr uint32_t lengthSoftClippedLeft() const noexcept { return _lengthSoftClippedLeft; };
	constexpr uint32_t lengthSoftClippedRight() const noexcept { return _lengthSoftClippedRight; };
	constexpr uint32_t lengthSoftClipped() const noexcept { return _lengthSoftClippedLeft + _lengthSoftClippedRight; };
	constexpr uint32_t lengthSequenced() const noexcept { return _lengthAligned + _lengthInserted; };
	constexpr uint32_t lengthRead() const noexcept {
		return _lengthAligned + _lengthInserted + _lengthSoftClippedLeft + _lengthSoftClippedRight;
	};
	uint32_t lengthMapped() const noexcept { return _lengthAligned + _lengthDeleted + _lengthSkipped; };
	std::string compileString() const;
	explicit operator std::string() const;
};

}; // namespace BAM

#endif /* BAM_TCIGAR_H_ */
