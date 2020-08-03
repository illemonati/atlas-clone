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
//forward declarations;

class TGenomeWindow;

class TGenomePosition{
protected:
	uint32_t _refID;
	uint32_t _position;

public:
	TGenomePosition(){ clear(); };
	TGenomePosition(const uint32_t& RefID, const uint32_t& Position);
	TGenomePosition(const TGenomePosition & other) = default;

	void clear(){  _refID = 0; _position = 0; };

	uint32_t refID() const{ return _refID; };
	uint32_t position() const{ return _position; };

	void move(const uint32_t& RefID, const uint32_t& Position);
	void move(const TGenomePosition & other);

	TGenomePosition& operator=(const TGenomePosition & other);
	TGenomePosition operator+(const uint32_t & length) const;
	TGenomePosition operator-(const uint32_t & length) const;
	int32_t operator-(const TGenomePosition & other) const;
	void operator+=(const uint32_t & length);
	void operator-=(const uint32_t & length);
	void operator++();
	void operator--();

	bool sameChr(const TGenomePosition & other) const;

	bool operator==(const TGenomePosition & other) const;
	bool operator<(const TGenomePosition & other) const;
	bool operator>(const TGenomePosition & other) const;
	bool operator>=(const TGenomePosition & other) const;
	bool operator<=(const TGenomePosition & other) const;
	bool operator<(const TGenomeWindow & other) const;
	bool operator>(const TGenomeWindow & other) const;
};

std::ostream& operator<<(std::ostream& os, const TGenomePosition & position);

//-----------------------------------------------------
// TGenomeWindow
//-----------------------------------------------------
class TGenomeWindow{
protected:
	TGenomePosition _from, _to; //NOTE: position _to is NOT included. Window is [from, to).

public:
	TGenomeWindow(){ clear(); };
	virtual ~TGenomeWindow(){ _from.move(0,0); _to.move(0, 1); };
	TGenomeWindow(const uint32_t& RefID, const uint32_t& From, const uint32_t& To);
	TGenomeWindow(const uint32_t& RefID, const uint32_t& From);
	TGenomeWindow(const TGenomePosition & position);
	TGenomeWindow(const TGenomePosition & From, const TGenomePosition & To);
	TGenomeWindow(const TGenomeWindow & other) = default;

	void clear();

	uint32_t refID() const{ return _from.refID(); };
	TGenomePosition from() const{ return _from; };
	TGenomePosition to() const{ return _to; };
	uint32_t fromOnChr() const{ return _from.position(); };
	uint32_t toOnChr() const{ return _to.position(); };
	uint32_t size() const{ return _to.position() - _from.position(); };

	virtual void move(const uint32_t& RefID, const uint32_t& Start, const uint32_t& End);
	virtual void move(const TGenomePosition & From, const uint32_t & Length);
	virtual void move(const TGenomePosition & From, const TGenomePosition & To);
	virtual void move(const TGenomeWindow & other);

	TGenomeWindow& operator=(const TGenomeWindow & other){ move(other); return *this; };
	TGenomeWindow operator+(const uint32_t& length) const;
	TGenomeWindow operator-(const uint32_t& length) const;

	bool sameChr(const TGenomeWindow & other) const{ return refID() == other.refID(); };
	bool within(const TGenomePosition & other) const;
	bool contains(const TGenomeWindow & other) const;
	bool overlaps(const TGenomeWindow & other) const;
	bool overlapsOrExtends(const TGenomeWindow & other) const;
	bool mergeWith(const TGenomeWindow & other);

	//move / expand
	virtual void operator+=(const uint32_t & length);
	virtual void operator-=(const uint32_t & length);
	virtual void resize(const uint32_t & newLength);

	bool operator<(const TGenomePosition & pos) const;
	bool operator>(const TGenomePosition & pos) const;
	bool operator<(const TGenomeWindow & other) const;
	bool operator>(const TGenomeWindow & other) const;
};

TGenomeWindow merge(const TGenomeWindow & first, const TGenomeWindow & second);
std::ostream& operator<<(std::ostream& os, const TGenomeWindow & window);

}; //to namespace

#endif /* BAM_TGENOMEPOSITION_H_ */
