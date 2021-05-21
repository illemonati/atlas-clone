/*
 * TBamFilter.h
 *
 *  Created on: May 23, 2020
 *      Author: phaentu
 */

#ifndef BAM_TBAMFILTER_H_
#define BAM_TBAMFILTER_H_

#include <set>
#include "TLog.h"
#include "TFile.h"
#include "TNumericRange.h"

namespace BAM{

//-----------------------------------------------------
//TBamFileLog
//-----------------------------------------------------
class TBamFileLog{
private:
	TOutputFile _log;

public:
	TBamFileLog(const std::string filename);
	void write(const std::string & alignmentName, const bool & isReverseStrand, const std::string & reason);
};

//-----------------------------------------------------
//TBamFileFilter
//-----------------------------------------------------
class TBamFileFilter{
protected:
	bool _keep;
	uint64_t _counter;
	std::string _reason; //used for reporting
	bool _updateLog;
	TBamFileLog* _log;

public:
	TBamFileFilter();
	void keep();
	bool filters() const{ return !_keep; };
	void setReason(const std::string reason);
	void setLog(TBamFileLog* Log);
	void filterOut(const std::string & alignmentName, const bool & isReverseStrand);
	void summary(TLog* logfile, const uint64_t & total);
	uint64_t numFiltered() const { return _counter; }
};

class TBamFileFilterBool:public TBamFileFilter{
public:
	TBamFileFilterBool(){};
	void filter(const std::string Reason);
	bool pass(const bool state, const std::string & alignmentName, const bool & isReverseStrand);
};

template <typename T>
class TBamFileFilterRange:public TBamFileFilterBool{
private:
	TNumericRange<T> _range;
public:
	TBamFileFilterRange(){};

	void filter(const TNumericRange<T> & Range, const std::string Reason){
		_keep = false;
		_range = Range;
		_reason = Reason;
	};

	/*
	void filter(const T & Min, const T & Max, const std::string Reason){
		_keep = false;
		_range.set(Min,  true,  Max,  true);
		_reason = Reason;
	};
	*/

	std::string rangeString() const{
		return _range.rangeString();
	};

	bool pass(const T & value, const std::string & alignmentName, const bool & isReverseStrand){
		if(!_keep && !_range.within(value)){
			filterOut(alignmentName, isReverseStrand);
			return false;
		}
		return true;
	};
};

//-----------------------------------------------------
//TAlignmentBlacklist
//-----------------------------------------------------
class TAlignmentList{
private:
	std::set<std::string> _list;

public:
	TAlignmentList(){};
	TAlignmentList(const std::string filename);
	void addFromFile(const std::string filename);
	void add(const std::string & alignment);
	void remove(const std::string & alignment);
	void clear();
	bool isInBlacklist(const std::string & alignment) const;
};

}; //end namespace

#endif /* BAM_TBAMFILTER_H_ */
