/*
 * TBamFilter.cpp
 *
 *  Created on: May 24, 2020
 *      Author: phaentu
 */

#include "TBamFilter.h"

#include <algorithm>
#include <vector>

#include "coretools/Files/TFile.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringFunctions.h"
#include "genometools/GenotypeTypes.h"

namespace BAM{

using coretools::instances::logfile;

//-----------------------------------------------------
//TBamFileFilter_base
//-----------------------------------------------------

void TBamFilter::filterOut(std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID){
	//counts filtered reads per read group and filter
	_counter.add(readGroup, chromosomeID);
	if(_log){
		_log->writeln(alignmentName, isSecondMate, _reason);
	}
};

void TBamFilter::keep(){
	_keep = true;
};

void TBamFilter::resizeCounter(size_t numRG, size_t numChrom){
	_counter.resize(numRG);
	_counter.resizeDistributions(numChrom);
}

void TBamFilter::setReason(std::string_view reason){
	_reason = reason;
};

void TBamFilter::setLog(coretools::TOutputFile& Log){
	_log = &Log;
};

void TBamFilter::summary(size_t total, size_t readGroup) const {
	if(!_keep && readGroup < _counter.size() && _counter[readGroup].counts()  > 0){
		logfile().list(_reason + ": ", _counter[readGroup].counts(), " (" + coretools::str::toPercentString(_counter[readGroup].counts(), total, 3) + "%)");
	}
};

void TBamFilter::fillHeader(std::vector<std::string> &header) const {
	std::string tmp = getReason();
	if (!tmp.empty()){
		std::replace(tmp.begin(), tmp.end(), ' ', '_');
		header.push_back(tmp);
	}
};

void TBamFilter::printCounts(coretools::TOutputFile &out, size_t rg_ID) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		out << getCounts(rg_ID);
	}
};

void TBamFilter::printCombinedCounts(coretools::TOutputFile &out) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		out << numFiltered().counts();
	}
};

size_t TBamFilter::getCounts(size_t rg_ID) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty() && rg_ID < numFiltered().size()){
		return numFiltered()[rg_ID].counts();
	}
	return 0;
};

size_t TBamFilter::getCountsPerChromosome(size_t ref_ID) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		return numFiltered().horizontalCounts(ref_ID);
	} return 0;
};

size_t TBamFilter::getCountsAtReadGroupAndChromosome(size_t rg_ID, size_t ref_ID) const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty() && rg_ID < numFiltered().size()){
		return numFiltered()[rg_ID][ref_ID];
	} return 0;
};

size_t TBamFilter::getCombinedCounts() const {
	//Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed reads per read group are printed here
	if (!getReason().empty()){
		return numFiltered().counts();
	} return 0;
};

//-----------------------------------------------------
//TBamFileFilterBool
//-----------------------------------------------------
void TBamFilterBool::filter(std::string_view Reason, size_t numRG, size_t numChrom){
	_keep = false;
	_reason = Reason;
	resizeCounter(numRG, numChrom);
};

bool TBamFilterBool::pass(bool state, std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID){
	if(!state && !_keep){
		filterOut(alignmentName, isSecondMate, readGroup, chromosomeID);
		return false;
	}
	return true;
};


}; //end namespace
