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
#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TSequencedBase.h"
#include "TParameters.h"
#include "TLog.h"

namespace BAM{

using coretools::TLog;
using coretools::TNumericRange;

//-----------------------------------------------------
//TBamFileLog
//-----------------------------------------------------
class TBamFileLog{
private:
	coretools::TOutputFile _log;

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
	void summary(TLog* logfile, uint64_t total);
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

	const TNumericRange<T> range() const {
		return _range;
	};

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

//-------------------------------------
// TBaseFilter
//-------------------------------------
class TBaseFilter{
protected:
	bool _filter;

public:
	explicit constexpr TBaseFilter() : _filter(false) {};
	virtual ~TBaseFilter() = default;

	constexpr operator bool() const{
		return _filter;
	};

	virtual bool pass(const TSequencedBase & base) const = 0;
};

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
class TQualityFilter : public TBaseFilter{
private:
	coretools::TNumericRange<genometools::PhredIntProbability> _range;

	genometools::PhredIntProbability _minPhredInt, _maxPhredInt;

	void _default();

public:
	explicit constexpr TQualityFilter(){
		_default();
	};

	TQualityFilter(coretools::TParameters & params, coretools::TLog* logfile){
		set(params, logfile);
	};

	~TQualityFilter() = default;

	void set(coretools::TParameters & params, coretools::TLog* logfile);

	//TODO: check if we filter on TSequencedBase, and if yes, on which error rate (original or recal)
	bool pass(const TSequencedBase & base) const override{
		return _range.within(base.recalibratedQualityAsPhredInt);
	};

	bool pass(const genometools::BaseQuality & qual) const{
		return _range.within(genometools::PhredIntProbability(qual));
	};
};

//-------------------------------------
// TContextFilter
//-------------------------------------
class TContextFilter : public TBaseFilter{
private:
	std::array<bool, static_cast<uint8_t>(genometools::BaseContext::max) + 1> _keptContexts;

public:
	explicit TContextFilter(){
		_keptContexts.fill(true);
	};
	~TContextFilter() = default;

	void set(coretools::TParameters & params, coretools::TLog* logfile);

	bool pass(const TSequencedBase & base) const override;
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
