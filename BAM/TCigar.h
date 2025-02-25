/*
 * TSamFlags.h
 *
 *  Created on: May 20, 2020
 *      Author: phaentu
 */

#ifndef BAM_TCIGAR_H_
#define BAM_TCIGAR_H_

#include <vector>
#include <string>

namespace BAM {

//-----------------------------------------------------
// TCigar
// A class to store, access and manipulate CIGAR operators
//-----------------------------------------------------
struct CigarOperator {
	char type       = '\0';
	size_t length = 0;

	CigarOperator(char Type, size_t Length) : type(Type), length(Length) {}
};

class TCigar {
private:
	std::vector<CigarOperator> _cigar;
	size_t _lengthAligned          = 0;
	size_t _lengthInserted         = 0;
	size_t _lengthDeleted          = 0;
	size_t _lengthSkipped          = 0;
	size_t _lengthSoftClippedLeft  = 0;
	size_t _lengthSoftClippedRight = 0;
	void _flipCigar();

public:
	TCigar()=default;
	TCigar(TCigar cigar, size_t overlapLength, bool isFirst, size_t &mappedBasesClipped);
	void clear();
	std::vector<CigarOperator>::const_iterator begin() const noexcept { return _cigar.begin(); }
	std::vector<CigarOperator>::const_iterator end() const noexcept { return _cigar.end(); }
	std::vector<CigarOperator>::const_reverse_iterator rbegin() const noexcept { return _cigar.rbegin(); }
	std::vector<CigarOperator>::const_reverse_iterator rend() const noexcept { return _cigar.rend(); }
	void add(char Type, size_t Length);
	void removeSoftClips();
	void removeSoftClips(size_t maxNumberOfSoftClippedBases);
	void addSoftClipsLeft(size_t Length);
	void addSoftClipsRight(size_t Length);
	void setAllSoftClipped();

	constexpr size_t lengthAligned() const noexcept { return _lengthAligned; };
	constexpr size_t lengthInserted() const noexcept { return _lengthInserted; };
	constexpr size_t lengthDeleted() const noexcept { return _lengthDeleted; };
	constexpr size_t lengthSoftClippedLeft() const noexcept { return _lengthSoftClippedLeft; };
	constexpr size_t lengthSoftClippedRight() const noexcept { return _lengthSoftClippedRight; };
	constexpr size_t lengthSoftClipped() const noexcept { return _lengthSoftClippedLeft + _lengthSoftClippedRight; };
	constexpr size_t lengthSequenced() const noexcept { return _lengthAligned + _lengthInserted; };
	constexpr size_t lengthRead() const noexcept {
		return _lengthAligned + _lengthInserted + _lengthSoftClippedLeft + _lengthSoftClippedRight;
	};
	size_t lengthMapped() const noexcept { return _lengthAligned + _lengthDeleted + _lengthSkipped; };
	std::string compileString() const;
	explicit operator std::string() const;
};

}; // namespace BAM

#endif /* BAM_TCIGAR_H_ */
