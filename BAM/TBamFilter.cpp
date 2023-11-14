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

namespace BAM {

using coretools::instances::logfile;

void TBamFilter::filterOut(std::string_view alignmentName, bool isSecondMate, size_t readGroup, int64_t chromosomeID, coretools::TOutputFile* const log) {
	// counts filtered reads per read group and filter
	_counter.add(readGroup, chromosomeID);
	if (log) { log->writeln(alignmentName, isSecondMate, _reason); }
}

void TBamFilter::summary(size_t total, size_t readGroup) const {
	if (filters() && readGroup < _counter.size() && _counter[readGroup].counts() > 0) {
		logfile().list(_reason + ": ", _counter[readGroup].counts(),
					   " (" + coretools::str::toPercentString(_counter[readGroup].counts(), total, 3) + "%)");
	}
}

void TBamFilter::fillHeader(std::vector<std::string> &header) const {
	if (filters()) {
		std::string tmp = _reason;
		std::replace(tmp.begin(), tmp.end(), ' ', '_');
		header.push_back(tmp);
	}
}

void TBamFilter::printCounts(coretools::TOutputFile &out, size_t rg_ID) const {
	// Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed
	// reads per read group are printed here
	if (filters()) { out << getCounts(rg_ID); }
}

void TBamFilter::printCombinedCounts(coretools::TOutputFile &out) const {
	// Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed
	// reads per read group are printed here
	if (filters()) { out << _counter.counts(); }
}

size_t TBamFilter::getCounts(size_t rg_ID) const {
	// Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed
	// reads per read group are printed here
	if (filters() && rg_ID < _counter.size()) { return _counter[rg_ID].counts(); }
	return 0;
}

size_t TBamFilter::getCountsPerChromosome(size_t ref_ID) const {
	// Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed
	// reads per read group are printed here
	if (filters()) { return _counter.horizontalCounts(ref_ID); }
	return 0;
}

size_t TBamFilter::getCountsAtReadGroupAndChromosome(size_t rg_ID, size_t ref_ID) const {
	// Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed
	// reads per read group are printed here
	if (filters() && rg_ID < _counter.size()) { return _counter[rg_ID][ref_ID]; }
	return 0;
}

size_t TBamFilter::getCombinedCounts() const {
	// Reason is only set if filter is applied (see TBamFile::setFilters), in which case reason and number of removed
	// reads per read group are printed here
	if (filters()) { return _counter.counts(); }
	return 0;
}

}; // namespace BAM
