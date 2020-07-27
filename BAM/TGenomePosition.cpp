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
TGenomePosition::TGenomePosition(const uint32_t RefID, const uint32_t Position){
	_refID = RefID;
	_position = Position;
}

TGenomePosition::TGenomePosition(const TGenomePosition & other) = default;

void TGenomePosition::clear(){
	_refID = 0;
	_position = 0;
};

void TGenomePosition::update(const uint32_t RefID, const uint32_t Position){
	_refID = RefID;
	_position = Position;
};

void TGenomePosition::update(const TGenomePosition & other){
	_refID = other._refID;
	_position = other._position;
};

void TGenomePosition::operator=(const TGenomePosition & other){
	update(other);
};

TGenomePosition TGenomePosition::operator+(const uint32_t length) const{
	return TGenomePosition(_refID, _position + length);
};

TGenomePosition TGenomePosition::operator-(const uint32_t length) const{
	if(length > _position){
		return TGenomePosition(_refID, 0);
	} else {
		return TGenomePosition(_refID, _position - length);
	}
};

int32_t TGenomePosition::operator-(const TGenomePosition other){
	if(!sameChr(other)){
		throw std::runtime_error("TGenomePosition operator-(const TGenomePosition other): positions are not on same chromosome!");
	}
	return _position - other._position;
};

void TGenomePosition::operator+=(const uint32_t length){
	_position += length;
};

void TGenomePosition::operator-=(const uint32_t length){
	if(length > _position){
		_position = 0;
	} else {
		_position -= length;
	}
};

bool TGenomePosition::sameChr(const TGenomePosition other) const{
	return _refID == other.refID();
};

bool TGenomePosition::operator<(const TGenomePosition other) const{
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

bool TGenomePosition::operator<(const TGenomeWindow other) const{
	return other > *this;
};

bool TGenomePosition::operator>(const TGenomeWindow other) const{
	return other < *this;
};

TOutputFile& TGenomePosition::operator<<(TOutputFile & out) const{
	out << _refID << _position;
	return out;
};

//-----------------------------------------------------
// TGenomeWindow
//-----------------------------------------------------
TGenomeWindow::TGenomeWindow(const uint32_t RefID, const uint32_t Start, const uint32_t End){
	update(RefID, Start, End);
};

//inserts a window of length one
TGenomeWindow::TGenomeWindow(const uint32_t RefID, const uint32_t Start){
	update(RefID, Start, Start+1);
};

TGenomeWindow::TGenomeWindow(const TGenomePosition & position){
	_refID = position.refID();
	_start = position.position();
	_end = position.position() + 1;
};

void TGenomeWindow::clear(){
	_refID = 0;
	_start=0;
	_end = 1;
};

void TGenomeWindow::update(const uint32_t RefID, const uint32_t Start, const uint32_t End){
	_refID = RefID;
	if(End <= Start){
		throw std::runtime_error("TGenomeWindow(const uint32_t RefID, const uint32_t Start, const uint32_t End): End <= Start!");
	}
	_start = Start;
	_end = End;
};

void TGenomeWindow::update(const TGenomeWindow & other){
	_refID = other._refID;
	_start = other._start;
	_end = other._end();
};

TGenomeWindow TGenomeWindow::operator+(const uint32_t length) const{
	return TGenomeWindow(_refID, _start + length, _end + length);
};

TGenomeWindow TGenomeWindow::operator-(const uint32_t length) const{
	if(length > _start){
		if(length > _end){
			return TGenomeWindow(_refID, 0, 1);
		} else {
			return TGenomeWindow(_refID, 0, _end - length);
		}
	} else {
		return TGenomeWindow(_refID, _start - length, _end - length);
	}
};

bool TGenomeWindow::within(const TGenomePosition other) const{
	//checks if other window is entirely within
	return _refID == other.refID() && _start <= other.position() && _end >= other.position();
};

bool TGenomeWindow::within(const TGenomeWindow other) const{
	//checks if other window is entirely within
	return _refID == other.refID() && _start <= other.start() && _end >= other.end();
};

bool TGenomeWindow::overlaps(const TGenomeWindow other) const{
	//check if other window overlaps
	if(refID != other.refID()){
		return false;
	}
	//on same chr: do they overlap?
	return (_start <= other._end && _start > other._end)
		|| (other._start <= _end && other._start > _start)
		|| within(other)
		|| other.within(*this);
};

TGenomeWindow TGenomeWindow::merge(const TGenomeWindow other){
	if(!overlaps(other)){
		throw std::runtime_error("TGenomeWindow merge(const TGenomeWindow other): windows do not overlap!");
	}
	return TGenomeWindow(_refID, std::min(_start, other._start), std::max(_end, other._end));
};

bool TGenomeWindow::mergeWith(const TGenomeWindow & other){
	if(!overlaps(other)){
		return false;
	}
	_start = std::min(_start, other._start);
	_end = std::max(_end, other._end);
	return true;
};

void TGenomeWindow::operator+=(const uint32_t length){
	_start += length;
	_end += length;
};

void TGenomeWindow::operator-=(const uint32_t length){
	if(length > _start){
		_start = 0;
		if(length > _end){
			_end = 1;
		} else {
			_end -= length;
		}
	} else {
		_start -= length;
		_end -= length;
	}
};

bool TGenomeWindow::operator<(const TGenomePosition & pos) const{
	//checks if position is after window end
	//are they on the same chromosome?
	if(this->_refID < pos.refID()){
		return true;
	} else if(this->_refID > pos.refID()){
		return false;
	} else {
		//on same chromosome: check position
		return this->_end <= pos._position;
	}
};

bool TGenomeWindow::operator>(const TGenomePosition & pos) const{
	//checks if position is before window start
	//are they on the same chromosome?
	if(this->_refID < pos._refID){
		return true;
	} else if(this->_refID > pos._refID){
		return false;
	} else {
		//on same chromosome: check position
		return this->_start > pos._position;
	}
};

bool TGenomeWindow::operator<(const TGenomeWindow & other) const{
	//checks if start is < other.start
	//are they on the same chromosome?
	if(this->_refID < other._refID){
		return true;
	} else if(this->_refID > other._refID){
		return false;
	} else {
		//on same chromosome: check position
		return this->_start < other._start;
	}
};

bool TGenomeWindow::operator>(const TGenomeWindow & other) const{
	//chesk if start is > other.start
	//are they on the same chromosome?
	if(this->_refID < other._refID){
		return true;
	} else if(this->_refID > other._refID){
		return false;
	} else {
		//on same chromosome: check position
		return this->_start > other._start;
	}
};

TOutputFile& TGenomeWindow::operator<<(TOutputFile & out) const{
	out << _refID << _start << _end;
	return out;
};

}; //end namespace


