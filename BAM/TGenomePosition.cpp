/*
 * TGenomePosition.cpp
 *
 *  Created on: Jun 10, 2020
 *      Author: phaentu
 */

#include "TGenomePosition.h"

namespace BAM{


//-----------------------------------------------------
// TGenomePosition
//-----------------------------------------------------
TGenomePosition::TGenomePosition(const uint32_t& RefID, const uint32_t& Position){
	_refID = RefID;
	_position = Position;
};

void TGenomePosition::move(const uint32_t& RefID, const uint32_t& Position){
	_refID = RefID;
	_position = Position;
};

void TGenomePosition::move(const TGenomePosition & other){
	_refID = other._refID;
	_position = other._position;
};

TGenomePosition& TGenomePosition::operator=(const TGenomePosition & other){
    // correct way to overload = is to return lhs by reference, see
    // https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md#c62-make-copy-assignment-safe-for-self-assignment
	move(other);
    return *this;
};

TGenomePosition TGenomePosition::operator+(const uint32_t & length) const{
	return TGenomePosition(_refID, _position + length);
};

TGenomePosition TGenomePosition::operator-(const uint32_t & length) const{
	if(length > _position){
		return TGenomePosition(_refID, 0);
	} else {
		return TGenomePosition(_refID, _position - length);
	}
};

int32_t TGenomePosition::operator-(const TGenomePosition & other) const{
	if(!sameChr(other)){
		throw std::runtime_error("TGenomePosition operator-(const TGenomePosition other): positions are not on same chromosome!");
	}
	if (other._position > _position)
	    return 0;
	return _position - other._position;
};

void TGenomePosition::operator+=(const uint32_t & length){
	_position += length;
};

void TGenomePosition::operator-=(const uint32_t & length){
	if(length > _position){
		_position = 0;
	} else {
		_position -= length;
	}
};

void TGenomePosition::operator++(){
	++_position;
};

void TGenomePosition::operator--(){
	if(_position > 0)
		--_position;
};

bool TGenomePosition::sameChr(const TGenomePosition & other) const{
	return _refID == other.refID();
};

bool TGenomePosition::operator==(const TGenomePosition & other) const{
	if(this->_refID == other._refID && this->_position == other._position){
		return true;
	} else {
		return false;
	}
};

bool TGenomePosition::operator<(const TGenomePosition & other) const{
	//are they on the same chromosome?
	if(this->_refID < other._refID){
		return true;
	} else if(this->_refID > other._refID){
		return false;
	} else {
		//on same chromosome: check position
		return this->_position < other._position;
	}
};

bool TGenomePosition::operator>(const TGenomePosition & other) const{
	//are they on the same chromosome?
	if(this->_refID > other._refID){
		return true;
	} else if(this->_refID < other._refID){
		return false;
	} else {
		//on same chromosome: check position
		return this->_position > other._position;
	}
};

bool TGenomePosition::operator>=(const TGenomePosition & other) const{
	if(*this == other){
		return true;
	}
	return *this > other;
};

bool TGenomePosition::operator<=(const TGenomePosition & other) const{
	if(*this == other){
		return true;
	}
	return *this < other;
};

bool TGenomePosition::operator<(const TGenomeWindow & other) const{
	return other > *this;
};

bool TGenomePosition::operator>(const TGenomeWindow & other) const{
	return other < *this;
};

std::ostream& operator<<(std::ostream& os, const TGenomePosition & position){
	os << position.refID() << ":" << position.position();
	return os;
};

//-----------------------------------------------------
// TGenomeWindow
//-----------------------------------------------------
TGenomeWindow::TGenomeWindow(const uint32_t& RefID, const uint32_t& From, const uint32_t& To){
	move(RefID, From, To);
};

//inserts a window of length one
TGenomeWindow::TGenomeWindow(const uint32_t& RefID, const uint32_t& From){
	move(TGenomePosition(RefID, From), 1);
};

TGenomeWindow::TGenomeWindow(const TGenomePosition & position){
	_from = position;
	_to = position + 1;
};

TGenomeWindow::TGenomeWindow(const TGenomePosition & From, const TGenomePosition & To){
	move(From, To);
};

void TGenomeWindow::clear(){
	_from.clear();
	_to = _from;
};

void TGenomeWindow::move(const uint32_t& RefID, const uint32_t& From, const uint32_t& To){
	if(To <= From){
		throw std::runtime_error("TGenomeWindow(const uint32_t RefID, const uint32_t Start, const uint32_t End): To <= From!");
	}
	_from.move(RefID, From);
	_to.move(RefID, To);
};

void TGenomeWindow::move(const TGenomePosition & From, const uint32_t & Length){
	_from = From;
	_to = _from + Length;
};

void TGenomeWindow::move(const TGenomePosition & From, const TGenomePosition & To){
	if(From.refID() != To.refID()){
		throw std::runtime_error("TGenomeWindow(const TGenomePosition & From, const TGenomePosition & To): Chromosomes do not match!");
	}
	if(To <= From){
		throw std::runtime_error("TGenomeWindow(const TGenomePosition & From, const TGenomePosition & To): To <= From!");
	}
	_from = From;
	_to = To;
};

void TGenomeWindow::move(const TGenomeWindow & other){
	_from = other.from();
	_to = other.to();
};

TGenomeWindow TGenomeWindow::operator+(const uint32_t& length) const{
	return TGenomeWindow(_from + length, _to + length);
};

TGenomeWindow TGenomeWindow::operator-(const uint32_t& length) const{
	if(length > _from.position()){
		if(length > _to.position()){
			return TGenomeWindow(_from.refID(), 0, 1);
		} else {
			return TGenomeWindow(_from.refID(), 0, _to.position() - length);
		}
	} else {
		return TGenomeWindow(_from.refID(), _from.position() - length, _to.position() - length);
	}
};

bool TGenomeWindow::within(const TGenomePosition & other) const{
	//checks if other position is entirely within
	return _from <= other && _to >= other;
};

bool TGenomeWindow::contains(const TGenomeWindow & other) const{
    //checks if other window is entirely within
	return _from.refID() == other.refID() && _from <= other.from() && _to >= other.to();
};

bool TGenomeWindow::overlaps(const TGenomeWindow & other) const{
	//check if other window overlaps
	if(_from.refID() != other.refID()){
		return false;
	}
	//on same chr: do they overlap?
	return (_from <= other.from() && _to > other.from()) || (other.from() <= _from && other.to() > _from);
};

bool TGenomeWindow::overlapsOrExtends(const TGenomeWindow & other) const{
	//std::cout << "this: [" << _from << ", " << _to << ") --- other: [" << other.from() << ", " << other.to() << ")" << std::endl;

	//check if other window overlaps, or are consecutive (no gap)
	if(_from.refID() != other.refID()){
		return false;
	}
	//on same chr: do they overlap?
	return (_from <= other.from() && _to >= other.from())
		|| (other.from() <= _from && other.to() >= _from);
};

bool TGenomeWindow::mergeWith(const TGenomeWindow & other){
	if(!overlapsOrExtends(other)){
		return false;
	}

	if(_from > other.from()){
		_from = other.from();
	}
	if(_to < other.to()){
		_to = other.to();
	}
	return true;
};

void TGenomeWindow::operator+=(const uint32_t & length){
	_from += length;
	_to += length;
};

void TGenomeWindow::operator-=(const uint32_t & length){
	if(length > _from.position()){
		_from.move(_from.refID(), 0);
		if(length > _to.position()){
			_to.move(_to.refID(), 1);
		} else {
			_to -= length;
		}
	} else {
		_from -= length;
		_to -= length;
	}
};

void TGenomeWindow::resize(const uint32_t & newLength){
	_to = _from + newLength;
};

bool TGenomeWindow::operator<(const TGenomePosition & pos) const{
	//checks if position is after window end
	return _to <= pos;
};

bool TGenomeWindow::operator>(const TGenomePosition & pos) const{
	//checks if position is before window start
	return pos < _from;
};

bool TGenomeWindow::operator<(const TGenomeWindow & other) const{
	//checks if from is < other.from
	return _from < other.from();
};

bool TGenomeWindow::operator>(const TGenomeWindow & other) const{
	//checks if from is > other.from
	return _from > other.from();
};

TGenomeWindow merge(const TGenomeWindow & first, const TGenomeWindow & second){
	if(!first.overlapsOrExtends(second)){
		throw std::runtime_error("TGenomeWindow merge(const TGenomeWindow other): windows do not overlap!");
	}
	return TGenomeWindow(std::min(first.from(), second.from()), std::max(first.to(), second.to()));
};

std::ostream& operator<<(std::ostream& os, const TGenomeWindow & window){
	os << window.from() << "-" << window.to();
	return os;
};

}; //end namespace


