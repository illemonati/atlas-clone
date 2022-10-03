/*
 * TBamFilter.cpp
 *
 *  Created on: May 24, 2020
 *      Author: phaentu
 */

#include "TBamFilter.h"

#include <algorithm>
#include <exception>
#include <ostream>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "coretools/Files/TFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringFunctions.h"

namespace BAM{

using coretools::TLog;

//-----------------------------------------------------
//TBamFileLog
//-----------------------------------------------------
void TBamFileLog::write(const std::string & alignmentName, const bool & isReverseStrand, const std::string & reason){
	_log.writeln(alignmentName, isReverseStrand, reason);
};

//-----------------------------------------------------
//TBamFileFilter_base
//-----------------------------------------------------
TBamFileFilter::TBamFileFilter(){
	_keep = true;
	_updateLog = false;
	_log = nullptr;
};

void TBamFileFilter::filterOut(const std::string & alignmentName, const bool & isReverseStrand, const uint16_t readGroup){
	//counts filtered reads per read group and filter
	_counter.add(readGroup);
	if(_updateLog){
		_log->write(alignmentName, isReverseStrand, _reason);
	}
};

void TBamFileFilter::keep(){
	_keep = true;
};

void TBamFileFilter::setReason(const std::string reason){
	_reason = reason;
};

void TBamFileFilter::setLog(std::shared_ptr<TBamFileLog> & Log){
	_log = Log;
	_updateLog = true;
};

void TBamFileFilter::summary(TLog* logfile, uint64_t total, const uint16_t readGroup){
	if(!_keep && _counter[readGroup]  > 0){
		logfile->list(_reason + ": ", _counter[readGroup], " (" + coretools::str::toPercentString(_counter[readGroup], total, 3) + "%)");
	}
};

void TBamFileFilter::fillHeader(std::vector<std::string> &header){
	if (!getReason().empty()){
			header.push_back(getReason());
		}
};

void TBamFileFilter::printCounts(coretools::TOutputFile &out, uint16_t rg_ID) {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		coretools::TCountDistribution FilterCount = numFiltered();
		out << FilterCount[rg_ID];
	}
};

void TBamFileFilter::printCombinedCounts(coretools::TOutputFile &out) {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		coretools::TCountDistribution FilterCount = numFiltered();
		out << FilterCount.counts();
	}
};

//-----------------------------------------------------
//TBamFileFilterBool
//-----------------------------------------------------
void TBamFileFilterBool::filter(const std::string Reason){
	_keep = false;
	_reason = Reason;
};

bool TBamFileFilterBool::pass(const bool state, const std::string & alignmentName, const bool & isReverseStrand, const uint16_t readGroup){
	if(!state && !_keep){
		filterOut(alignmentName, isReverseStrand, readGroup);
		return false;
	}
	return true;
};


//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
void TQualityFilter::_default(){
	//default values according to SAM specifications
	_range.set(genometools::PhredIntProbability(1), true, genometools::PhredIntProbability(93), true);
};

void TQualityFilter::set(coretools::TParameters & params, coretools::TLog* logfile){
	if(params.parameterExists("filterBaseQual")){
		params.fillParameter("filterBaseQual", _range);
		if (_range.within(genometools::PhredIntProbability(0))){ throw "Base quality filter of 0 is not allowed (parameter 'filterBaseQual')"; }
		logfile->list("Will filter out bases with quality outside the range " + _range.rangeString() + " (parameter 'filterBaseQual')");
	} else {
		_default();
		logfile->list("Will filter out bases with quality outside the range " + _range.rangeString() + ". (use 'filterBaseQual' to change)");
	}
	_filter = true;
};

//-------------------------------------
// TContextFilter
//-------------------------------------
void TContextFilter::set(coretools::TParameters & params, coretools::TLog* logfile){
	using namespace genometools;
	_filter = false;
	if(params.parameterExists("ignoreContexts")){
		std::vector<std::string> contexts;
		params.fillParameterIntoContainer("ignoreContexts", contexts, ',');

		if(contexts.size() > 0){
			for(auto& c : contexts){
				if(c.size() != 2){
					throw "Context " + c + " does not consist of two bases! (parameter 'ignoreContexts')";
				}

				const Base first  = char2base(c[0]);
				const Base second = char2base(c[1]);

				if(base2char(first) != c[0] || base2char(second) != c[1]){
					throw "Unable to understand context '" + c + "'!  (parameter 'ignoreContexts')";
				}

				//save context
				_keptContexts[baseContext(first, second)] = false;
			}

			std::vector<std::string> rep;
			for(auto i = BaseContext::min; i <= BaseContext::max; ++i){ //including max
				if(!_keptContexts[i]){
					rep.push_back(toString(i));
				}
			}
			logfile->list("Will ignore the following contexts: " + coretools::str::concatenateString(rep, ", ")  + ". (parameter 'ignoreContexts')");
			_filter = true;
		}
	}

	if(!_filter){
		logfile->list("Will keep bases regardless of base context. (use 'ignoreContexts' to filter)");
	}
};

bool TContextFilter::pass(const TSequencedBase & base) const{
	return _keptContexts[base.context()];
};

//-----------------------------------------------------
//TAlignmentBlacklist
//-----------------------------------------------------
TAlignmentList::TAlignmentList(const std::string filename){
	addFromFile(filename);
};

void TAlignmentList::addFromFile(const std::string filename){
	coretools::TInputFile in(filename, 1);
	std::vector<std::string> vec;
	while(in.read(vec)){
		add(vec[0]);
	}
};

void TAlignmentList::add(const std::string & alignment){
	_list.insert(alignment);
};

void TAlignmentList::remove(const std::string & alignment){
	_list.erase(alignment);
};

void TAlignmentList::clear(){
	_list.clear();
};

bool TAlignmentList::isInBlacklist(const std::string & alignment) const{
	return _list.find(alignment) != _list.end();
};

}; //end namespace
