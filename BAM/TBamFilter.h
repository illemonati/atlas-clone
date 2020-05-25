/*
 * TBamFilter.h
 *
 *  Created on: May 23, 2020
 *      Author: phaentu
 */

#ifndef BAM_TBAMFILTER_H_
#define BAM_TBAMFILTER_H_

#include "TLog.h"
#include "TAlignmentBlacklist.h"

namespace BAM{

//-----------------------------------------------------
//TBamFileFilter
//-----------------------------------------------------
class TBamFileFilter_base{
protected:
	bool _keep;
	uint64_t _counter;
	std::string _errorString;
	bool _updateBlacklist;
	TAlignmentBlacklist* _blacklist;

	void _filterOut(const std::string & alignmentName, const bool & isReverseStrand){
		++_counter;
		if(_updateBlacklist){
			_blacklist->addToBlacklist(alignmentName, isReverseStrand, _errorString);
		}
	};


public:
	TBamFileFilter_base(){
		_counter = 0;
		_keep = true;
		_updateBlacklist = false;
		_blacklist = nullptr;
	};

	void keep(){
		_keep = true;
	};

	void setBlacklist(const TAlignmentBlacklist* Blacklist){
		_blacklist = Blacklist;
		_updateBlacklist = true;
	};

	void summary(TLog* logfile){
		if(!_keep){
			logfile->list(_errorString + ": " + toString(_counter));
		}
	};
};

class TBamFileFilterBool:public TBamFileFilter_base{
public:
	TBamFileFilterBool(){};

	void filter(const std::string ErrorString){
		_keep = false;
		_errorString = ErrorString;
	};

	bool pass(const bool state, const std::string & alignmentName, const bool & isReverseStrand){
		if(state && !_keep){
			_filterOut(alignmentName, isReverseStrand);
			return false;
		}
		return true;
	};
};

class TBamFileFilterRange:public TBamFileFilterBool{
private:
	uint32_t _min, _max;
public:
	TBamFileFilterRange(){
		_min = 0;
		_max = UINT32_MAX;
	};

	void filter(const uint32_t Min, const uint32_t Max, const std::string ErrorString){
		_keep = false;
		_min = Min;
		_max = Max;
		_errorString = ErrorString;
	};

	bool pass(const uint32_t value, const std::string & alignmentName, const bool & isReverseStrand){
		if(!_keep && (value < _min || value > _max)){
			_filterOut(alignmentName, isReverseStrand);
		}
		return true;
	};

};

}; //end namespace

#endif /* BAM_TBAMFILTER_H_ */
