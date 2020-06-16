/*
 * TGenomePosition.h
 *
 *  Created on: Jun 9, 2020
 *      Author: phaentu
 */

#ifndef BAM_TGENOMEPOSITION_H_
#define BAM_TGENOMEPOSITION_H_

#include <cstdint>
#include <set>
#include "TFile.h"

namespace BAM{


//-----------------------------------------------------
// TGenomePosition
//-----------------------------------------------------
//forward declerations;

class TGenomeWindow;

class TGenomePosition{
private:
	uint32_t _refID;
	uint32_t _position;

public:
	TGenomePosition(){ clear(); };
	TGenomePosition(const uint32_t RefID, const uint32_t Position);
	TGenomePosition(const TGenomePosition & other) = default;

	void clear();
	uint32_t refID() const{ return _refID; };
	uint32_t position() const{ return _position; };

	void update(const uint32_t RefID, const uint32_t Position);
	void update(const TGenomePosition & other);

	void operator=(const TGenomePosition & other);
	TGenomePosition operator+(const uint32_t length) const;
	TGenomePosition operator-(const uint32_t length) const;
	int32_t operator-(const TGenomePosition other);
	void operator+=(const uint32_t length);
	void operator-=(const uint32_t length);

	bool sameChr(const TGenomePosition other) const;

	bool operator<(const TGenomePosition other) const;
	bool operator<(const TGenomeWindow other) const;
	bool operator>(const TGenomeWindow other) const;
	TOutputFile& operator<<(TOutputFile & out) const;
};

//-----------------------------------------------------
// TGenomeWindow
//-----------------------------------------------------
class TGenomeWindow{
private:
	uint32_t _refID;
	uint32_t _start, _end; //end not included: [start, end)

public:
	TGenomeWindow(){ clear(); };
	TGenomeWindow(const uint32_t RefID, const uint32_t Start, const uint32_t End);
	TGenomeWindow(const uint32_t RefID, const uint32_t Start);
	TGenomeWindow(const TGenomePosition & position);
	TGenomeWindow(const TGenomeWindow & other) = default;
	void clear();

	uint32_t refID() const{ return _refID; };
	uint32_t start() const{ return _start; };
	uint32_t end() const{ return _end; };
	uint32_t size() const{ return _end - _start; };
	TGenomePosition startPos(){ return TGenomePosition(_refID, _start); };

	void update(const uint32_t RefID, const uint32_t Start, const uint32_t End);
	void update(const TGenomeWindow & other);

	void operator=(const TGenomeWindow & other){ update(other); };
	TGenomeWindow operator+(const uint32_t length) const;
	TGenomeWindow operator-(const uint32_t length) const;

	bool sameChr(const TGenomeWindow other) const{ return _refID == other.refID(); };
	bool within(const TGenomePosition other) const;
	bool within(const TGenomeWindow other) const;
	bool overlaps(const TGenomeWindow other) const;
	TGenomeWindow merge(const TGenomeWindow other);
	bool mergeWith(const TGenomeWindow & other);

	void operator+=(const uint32_t length);
	void operator-=(const uint32_t length);
	bool operator<(const TGenomePosition & pos) const;
	bool operator>(const TGenomePosition & pos) const;
	bool operator<(const TGenomeWindow & other) const;
	bool operator>(const TGenomeWindow & other) const;

	TOutputFile& operator<<(TOutputFile & out) const;
};


}; //end namespace

#endif /* BAM_TGENOMEPOSITION_H_ */
