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
	void summary(TLog* logfile);
};

class TBamFileFilterBool:public TBamFileFilter{
public:
	TBamFileFilterBool(){};
	void filter(const std::string Reason);
	bool pass(const bool state, const std::string & alignmentName, const bool & isReverseStrand);
};

class TBamFileFilterRange:public TBamFileFilterBool{
private:
	uint32_t _min, _max;
public:
	TBamFileFilterRange();
	void filter(const uint32_t Min, const uint32_t Max, const std::string Reason);
	bool pass(const uint32_t value, const std::string & alignmentName, const bool & isReverseStrand);
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
