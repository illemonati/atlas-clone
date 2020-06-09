/*
 * TGenomePosition.h
 *
 *  Created on: Jun 9, 2020
 *      Author: phaentu
 */

#ifndef BAM_TGENOMEPOSITION_H_
#define BAM_TGENOMEPOSITION_H_

#include <cstdint>

namespace BAM{

//-----------------------------------------------------
// TGenomePosition
//-----------------------------------------------------
class TGenomePosition{
private:
	uint32_t _refID;
	uint32_t _position;
	bool _isSet;
public:
	TGenomePosition(){
		clear();
	};
	TGenomePosition(const uint32_t RefID, const uint32_t Position){
		_refID = RefID;
		_position = Position;
		_isSet = true;
	};
	void clear(){
		_refID = 0;
		_position = 0;
		_isSet = false;
	};
	bool isSet(){ return _isSet; };
	uint32_t refId() const{ return _refID; };
	uint32_t position() const{ return _position; };
	void update(const uint32_t RefID, const uint32_t Position){
		_refID = RefID;
		_position = Position;
		_isSet = true;
	};
	friend bool operator<(const TGenomePosition & other){
		if(!_isSet){
			throw "Can not compare TGenomicPosition: position not set!";
		}
		//are they on the same chromosome?
		if(this->_refID < other._refID){
			return true;
		} else if(this->_refID > other._refID){
			return false;
		} else {
			//on same chromosome: check posotion
			return this->_position < other._position;
		}
	};
};

}; //end namespace

#endif /* BAM_TGENOMEPOSITION_H_ */
