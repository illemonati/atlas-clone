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

using coretools::instances::parameters;
using coretools::instances::logfile;

//-----------------------------------------------------
//TBamFileLog
//-----------------------------------------------------
void TBamFileLog::write(const std::string & alignmentName, const bool & isSecondMate, const std::string & reason){
	_log.writeln(alignmentName, isSecondMate, reason);
};

//-----------------------------------------------------
//TBamFileFilter_base
//-----------------------------------------------------
TBamFileFilter::TBamFileFilter(){
	_keep = true;
	_updateLog = false;
	_log = nullptr;
};

void TBamFileFilter::filterOut(const std::string & alignmentName, const bool & isSecondMate, const uint16_t readGroup, const uint32_t chromosomeID){
	//counts filtered reads per read group and filter
	_counter.add(readGroup, chromosomeID);
	if(_updateLog){
		_log->write(alignmentName, isSecondMate, _reason);
	}
};

void TBamFileFilter::keep(){
	_keep = true;
};

void TBamFileFilter::resizeCounter(uint16_t numRG, uint32_t numChrom){
	_counter.resize(numRG);
	_counter.resizeDistributions(numChrom);
}

void TBamFileFilter::setReason(const std::string reason){
	_reason = reason;
};

void TBamFileFilter::setLog(std::shared_ptr<TBamFileLog> & Log){
	_log = Log;
	_updateLog = true;
};

void TBamFileFilter::summary(uint64_t total, const uint16_t readGroup){
	if(!_keep && _counter[readGroup].counts()  > 0){
		logfile().list(_reason + ": ", _counter[readGroup].counts(), " (" + coretools::str::toPercentString(_counter[readGroup].counts(), total, 3) + "%)");
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
		coretools::TCountDistributionVector FilterCount = numFiltered();
		out << FilterCount[rg_ID].counts();
	}
};

void TBamFileFilter::printCountsPerChromosome(coretools::TOutputFile &out, uint32_t ref_ID) {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		coretools::TCountDistributionVector FilterCount = numFiltered();
		out << FilterCount.horizontalCounts(ref_ID);
	}
};

void TBamFileFilter::printCombinedCounts(coretools::TOutputFile &out) {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		coretools::TCountDistributionVector FilterCount = numFiltered();
		out << FilterCount.counts();
	}
};

size_t TBamFileFilter::getCounts(uint16_t rg_ID) {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		coretools::TCountDistributionVector FilterCount = numFiltered();
		return(FilterCount[rg_ID].counts());
	} return 0;
};

size_t TBamFileFilter::getCountsPerChromosome(uint32_t ref_ID) {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		coretools::TCountDistributionVector FilterCount = numFiltered();
		return(FilterCount.horizontalCounts(ref_ID));
	} return 0;
};

size_t TBamFileFilter::getCountsAtReadGroupAndChromosome(uint16_t rg_ID, uint32_t ref_ID) {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		coretools::TCountDistributionVector FilterCount = numFiltered();
		return(FilterCount[rg_ID][ref_ID]);
	} return 0;
};

size_t TBamFileFilter::getCombinedCounts() {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		coretools::TCountDistributionVector FilterCount = numFiltered();
		return(FilterCount.counts());
	} return 0;
};

//-----------------------------------------------------
//TBamFileFilterBool
//-----------------------------------------------------
void TBamFileFilterBool::filter(const std::string Reason, uint16_t numRG, uint32_t numChrom){
	_keep = false;
	_reason = Reason;
	resizeCounter(numRG, numChrom);
};

bool TBamFileFilterBool::pass(const bool state, const std::string & alignmentName, const bool & isSecondMate, const uint16_t readGroup, const uint32_t chromosomeID){
	if(!state && !_keep){
		filterOut(alignmentName, isSecondMate, readGroup, chromosomeID);
		return false;
	}
	return true;
};


//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
TQualityFilter::TQualityFilter() {
	if(parameters().parameterExists("filterBaseQual")){
		parameters().fillParameter("filterBaseQual", _range);
		if (_range.within(genometools::PhredIntProbability(0))){ UERROR("Base quality filter of 0 is not allowed (parameter 'filterBaseQual')"); }
		logfile().list("Will filter out bases with quality outside the range " + _range.rangeString() + " (parameter 'filterBaseQual')");
	} else {
		_range.set(genometools::PhredIntProbability(1), true, genometools::PhredIntProbability(93), true);
		logfile().list("Will filter out bases with quality outside the range " + _range.rangeString() + ". (use 'filterBaseQual' to change)");
	}
	_filter = true;
};

//-------------------------------------
// TContextFilter
//-------------------------------------
TContextFilter::TContextFilter(){
	using namespace genometools;

	_keptContexts.fill(true);
	_filter = false;
	if(parameters().parameterExists("ignoreContexts")){
		std::vector<std::string> contexts;
		parameters().fillParameterIntoContainer("ignoreContexts", contexts, ',');

		if(contexts.size() > 0){
			for(auto& c : contexts){
				if(c.size() != 2){
					UERROR("Context ", c, " does not consist of two bases! (parameter 'ignoreContexts')");
				}

				const Base first  = char2base(c[0]);
				const Base second = char2base(c[1]);

				if(base2char(first) != c[0] || base2char(second) != c[1]){
					UERROR("Unable to understand context '", c, "'!  (parameter 'ignoreContexts')");
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
			logfile().list("Will ignore the following contexts: " + coretools::str::concatenateString(rep, ", ")  + ". (parameter 'ignoreContexts')");
			_filter = true;
		}
	}

	if(!_filter){
		logfile().list("Will keep bases regardless of base context. (use 'ignoreContexts' to filter)");
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
