/*
 * TBamFilter.h
 *
 *  Created on: May 23, 2020
 *      Author: phaentu
 */

#ifndef BAM_TBAMFILTER_H_
#define BAM_TBAMFILTER_H_

#include <TMateFinder.h>
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
class TBamFileFilter_base{
protected:
	bool _keep;
	uint64_t _counter;
	std::string _reason; //used for reporting
	bool _updateLog;
	TBamFileLog* _log;

	void _filterOut(const std::string & alignmentName, const bool & isReverseStrand);

public:
	TBamFileFilter_base();
	void keep();
	void setLog(const TBamFileLog* Log);
	void summary(TLog* logfile);
};

class TBamFileFilterBool:public TBamFileFilter_base{
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
class TAlignmentBlacklist{
private:
	std::set<std::string> _blacklist;

public:
	TAlignmentBlacklist(){};
	TAlignmentBlacklist(const std::string filename);
	void addFromFile(const std::string filename);
	void add(const std::string & alignment);
	void remove(const std::string & alignment);
	bool isInBlacklist(const std::string & alignment);
};

}; //end namespace

#endif /* BAM_TBAMFILTER_H_ */
