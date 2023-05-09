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

#include "coretools/Files/TOutputFile.h"
#include "coretools/devtools.h"
#include "genometools/GenotypeTypes.h"
#include "coretools/Files/TFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringFunctions.h"

namespace BAM{

using coretools::instances::parameters;
using coretools::instances::logfile;

//-----------------------------------------------------
//TBamFileFilter_base
//-----------------------------------------------------

void TBamFileFilter::filterOut(std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID){
	//counts filtered reads per read group and filter
	_counter.add(readGroup, chromosomeID);
	if(_log){
		_log->writeln(alignmentName, isSecondMate, _reason);
	}
};

void TBamFileFilter::keep(){
	_keep = true;
};

void TBamFileFilter::resizeCounter(size_t numRG, size_t numChrom){
	_counter.resize(numRG);
	_counter.resizeDistributions(numChrom);
}

void TBamFileFilter::setReason(std::string_view reason){
	_reason = reason;
};

void TBamFileFilter::setLog(coretools::TOutputFile& Log){
	_log = &Log;
};

void TBamFileFilter::summary(size_t total, size_t readGroup) const {
	if(!_keep && readGroup < _counter.size() && _counter[readGroup].counts()  > 0){
		logfile().list(_reason + ": ", _counter[readGroup].counts(), " (" + coretools::str::toPercentString(_counter[readGroup].counts(), total, 3) + "%)");
	}
};

void TBamFileFilter::fillHeader(std::vector<std::string> &header) const {
	std::string tmp = getReason();
	if (!tmp.empty()){
		std::replace(tmp.begin(), tmp.end(), ' ', '_');
		header.push_back(tmp);
	}
};

void TBamFileFilter::printCounts(coretools::TOutputFile &out, size_t rg_ID) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		out << getCounts(rg_ID);
	}
};

void TBamFileFilter::printCountsPerChromosome(coretools::TOutputFile &out, size_t ref_ID) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		out << numFiltered().horizontalCounts(ref_ID);
	}
};

void TBamFileFilter::printCombinedCounts(coretools::TOutputFile &out) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		out << numFiltered().counts();
	}
};

size_t TBamFileFilter::getCounts(size_t rg_ID) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty() && rg_ID < numFiltered().size()){
		return numFiltered()[rg_ID].counts();
	}
	return 0;
};

size_t TBamFileFilter::getCountsPerChromosome(size_t ref_ID) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		return numFiltered().horizontalCounts(ref_ID);
	} return 0;
};

size_t TBamFileFilter::getCountsAtReadGroupAndChromosome(size_t rg_ID, size_t ref_ID) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty() && rg_ID < numFiltered().size()){
		return numFiltered()[rg_ID][ref_ID];
	} return 0;
};

size_t TBamFileFilter::getCombinedCounts() const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		return numFiltered().counts();
	} return 0;
};

//-----------------------------------------------------
//TBamFileFilterBool
//-----------------------------------------------------
void TBamFileFilterBool::filter(std::string_view Reason, size_t numRG, size_t numChrom){
	_keep = false;
	_reason = Reason;
	resizeCounter(numRG, numChrom);
};

bool TBamFileFilterBool::pass(bool state, std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID){
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
TAlignmentList::TAlignmentList(std::string_view filename){
	addFromFile(filename);
};

void TAlignmentList::addFromFile(std::string_view filename){
	coretools::TInputFile in(filename, 1);
	std::vector<std::string> vec;
	while(in.read(vec)){
		add(vec[0]);
	}
};

void TAlignmentList::add(std::string_view alignment){
	_list.emplace(alignment);
};

void TAlignmentList::remove(const std::string& alignment){
	_list.erase(alignment);
};

void TAlignmentList::clear(){
	_list.clear();
};

bool TAlignmentList::isInBlacklist(const std::string & alignment) const{
	return _list.find(alignment) != _list.end();
};

}; //end namespace
