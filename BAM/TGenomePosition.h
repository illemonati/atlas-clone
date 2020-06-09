/*
 * TGenomePosition.h
 *
 *  Created on: Jun 9, 2020
 *      Author: phaentu
 */

#ifndef BAM_TGENOMEPOSITION_H_
#define BAM_TGENOMEPOSITION_H_

#include <cstdint>
#include "TFile.h"
#include "TChromosomes.h"

namespace BAM{

class TGenomeWindow;

//-----------------------------------------------------
// TGenomePosition
//-----------------------------------------------------
class TGenomePosition{
private:
	uint32_t _refID;
	uint32_t _position;

public:
	TGenomePosition(){
		clear();
	};

	TGenomePosition(const uint32_t RefID, const uint32_t Position){
		_refID = RefID;
		_position = Position;
	}

	TGenomePosition(const TGenomePosition & other) = default;

	void clear(){
		_refID = 0;
		_position = 0;
	};

	uint32_t refId() const{ return _refID; };

	uint32_t position() const{ return _position; };

	void update(const uint32_t RefID, const uint32_t Position){
		_refID = RefID;
		_position = Position;
	};
	void update(const TGenomePosition & other){
		_refID = other._refID;
		_position = other._position;
	};

	void operator=(const TGenomePosition & other){
		update(other);
	};

	TGenomePosition operator+(const uint32_t length) const{
		return TGenomePosition(_refID, _position + length);
	};

	TGenomePosition operator-(const uint32_t length) const{
		if(length > _position){
			return TGenomePosition(_refID, 0);
		} else {
			return TGenomePosition(_refID, _position - length);
		}
	};

	int32_t operator-(const TGenomePosition other){
		if(!sameChr(other)){
			throw std::runtime_error("TGenomePosition operator-(const TGenomePosition other): positions are not on same chromosome!");
		}
		return _position - other._position;
	};

	bool sameChr(const TGenomePosition other) const{
		return _refID == other.refId();
	};

	void operator+=(const uint32_t length){
		_position += length;
	};

	void operator-=(const uint32_t length){
		if(length > _position){
			_position = 0;
		} else {
			_position -= length;
		}
	};

	friend bool operator<(const TGenomePosition other) const{
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

	friend bool operator<(const TGenomeWindow other) const{
		//are they on the same chromosome?
		if(this->_refID < other._refID){
			return true;
		} else if(this->_refID > other._refID){
			return false;
		} else {
			//on same chromosome: check position
			return this->_position < other._start;
		}
	};

	TOutputFile& operator<<(TOutputFile & out) const{
	    out << _refID << _position;
	    return out;
	};

	TOutputFile& operator<<(TOutputFile & out, const TChromosomes & Chromosomes) const{
		out << Chromosomes.name(_refID) << _position;
		return out;
	};
};

//-----------------------------------------------------
// TGenomeWindow
//-----------------------------------------------------
class TGenomeWindow{
private:
	uint32_t _refID;
	uint32_t _start, _end; //end not included: [start, end)

public:
	TGenomeWindow(){
		clear();
	};

	TGenomeWindow(const uint32_t RefID, const uint32_t Start, const uint32_t End){
		update(RefID, Start, End);
	};

	TGenomeWindow(const TGenomeWindow & other) = default;

	void clear(){
		_refID = 0;
		_start=0;
		_end = 1;
	};

	uint32_t refId() const{ return _refID; };
	uint32_t start() const{ return _start; };
	uint32_t end() const{ return _end; };

	void update(const uint32_t RefID, const uint32_t Start, const uint32_t End){
		_refID = RefID;
		if(End <= Start){
			throw std::runtime_error("TGenomeWindow(const uint32_t RefID, const uint32_t Start, const uint32_t End): End <= Start!");
		}
		_start = Start;
		_end = End;
	};

	void update(const TGenomeWindow & other){
		_refID = other._refID;
		_start = other._start;
		_end = other._end();
	};

	void operator=(const TGenomeWindow & other){
		update(other);
	};

	TGenomeWindow operator+(const uint32_t length) const{
		return TGenomeWindow(_refID, _start + length, _end + length);
	};

	TGenomeWindow operator-(const uint32_t length) const{
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

	bool sameChr(const TGenomeWindow other) const{
		return _refID == other.refId();
	};

	void operator+=(const uint32_t length){
		_start += length;
		_end += length;
	};

	void operator-=(const uint32_t length){
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

	friend bool operator<(const TGenomePosition & pos) const{
		//are they on the same chromosome?
		if(this->_refID < pos._refID){
			return true;
		} else if(this->_refID > pos._refID){
			return false;
		} else {
			//on same chromosome: check position
			return this->_end <= pos._position;
		}
	};

	friend bool operator<(const TGenomeWindow & other) const{
		//are they on the same chromosome?
		if(this->_refID < other._refID){
			return true;
		} else if(this->_refID > other._refID){
			return false;
		} else {
			//on same chromosome: check position
			return this->start <= other._start;
		}
	};

	TOutputFile& operator<<(TOutputFile & out) const{
	    out << _refID << _start << _end;
	    return out;
	};

	TOutputFile& operator<<(TOutputFile & out, const TChromosomes & Chromosomes) const{
		out << Chromosomes.name(_refID) << _start << _end;
		return out;
	};
};

}; //end namespace

#endif /* BAM_TGENOMEPOSITION_H_ */
