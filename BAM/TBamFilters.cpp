#include "TBamFilters.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TNumericRange.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/enum.h"

namespace BAM {

using coretools::TNumericRange;
using coretools::instances::logfile;
using coretools::instances::parameters;

TBamFilters::TBamFilters(bool Enable) {
	enable(FilterType::MappedLength, coretools::TNumericRange<size_t>(0, true, 500, true),
		   "MappedLength outside [0, 500]");

	// MappedLength doesn't count as enabled
	_enabled = false;

	if (!Enable) return;

	// alignment filters
	logfile().startIndent("Will use the following filters on reads:");

	// mapping length
	//--------------
	// is relevant for storage
	// print error if reads are longer and filter is default
	TNumericRange<size_t> mappingLengthRange;
	if (parameters().exists("filterMappingLength")) {
		parameters().fill("filterMappingLength", mappingLengthRange);
	} else {
		// set default
		mappingLengthRange.set(0, true, 500, true);
	}
	enable(FilterType::MappedLength, mappingLengthRange, "MappedLengthOutside" + mappingLengthRange.rangeString());

	logfile().list("Mapped length: restrict to range " + mappingLengthRange.rangeString() +
				   ". (parameter 'filterMappingLength')");
	if (mappingLengthRange.max() > 100000) {
		logfile().warning("The chosen mapping length filter allows for reads to span >100kb of the reference genome. "
						  "This may affect performance in case of paired-end reads.");
	}

	// keep all otherwise?
	//-------------------
	if (parameters().exists("keepAllReads")) {
		logfile().list("Will keep all reads. (parameter 'keepAllReads', overrules any other QC filter except "
					   "filterMappingLength)");
	} else {
		// duplicates
		if (parameters().exists("keepDuplicates")) {
			disable(FilterType::Duplicate);
			logfile().list("Duplicate reads: keep. (parameter 'keepDuplicates')");
		} else {
			enable(FilterType::Duplicate, "Duplicate");
			logfile().list("Duplicate reads: filter out. (use 'keepDuplicates' to keep)");
		}

		// soft clips
		if (parameters().exists("filterSoftClips")) {
			enable(FilterType::SoftClippedRation, "Soft clipped");
			if (parameters().get("filterSoftClips").empty()) {
				_softClipRatio = 0.;
				logfile().list("Soft clipped reads: filter out. (parameter 'filterSoftClips')");
			} else {
				_softClipRatio = parameters().get<double>("filterSoftClips");
				logfile().list("Soft clipped reads: filter out if softClipLength/readLength > ", _softClipRatio,
							   ". (parameter 'filterSoftClips')");
			}
		} else {
			disable(FilterType::SoftClippedRation);
			logfile().list("Soft clipped reads: keep. (use 'filterSoftClips' to filter out)");
		}

		// improper pairs
		if (parameters().exists("keepImproperPairs")) {
			disable(FilterType::ImproperPairs);
			logfile().list("Improper pairs: keep. (parameter 'keepImproperPairs')");
		} else {
			enable(FilterType::ImproperPairs, "ImproperPair");
			logfile().list("Improper pairs: filter out. (use 'keepImproperPairs' to keep)");
		}

		// unmapped reads
		if (parameters().exists("keepUnmappedReads")) {
			disable(FilterType::Unmapped);
			logfile().list("Unmapped reads: keep. (parameter 'keepUnmappedReads')");
		} else {
			enable(FilterType::Unmapped, "Unmapped");
			logfile().list("Unmapped reads: filter out. (use 'keepUnmappedReads' to keep)");
		}

		// failed QC
		if (parameters().exists("keepFailedQC")) {
			disable(FilterType::FailedQC);
			logfile().list("Failed QC: keep. (parameter 'keepFailedQC')");
		} else {
			enable(FilterType::FailedQC, "FailedQC");
			logfile().list("Failed QC: filter out. (use 'keepFailedQC' to keep)");
		}

		// secondary reads
		if (parameters().exists("keepSecondaryReads")) {
			disable(FilterType::Secondary);
			logfile().list("Secondary reads: keep. (parameter 'keepSecondaryReads')");
		} else {
			enable(FilterType::Secondary, "SecondaryAlignment");
			logfile().list("Secondary reads: filter out. (use 'keepSecondaryReads' to keep)");
		}

		// supplementary reads
		if (parameters().exists("keepSupplementaryReads")) {
			disable(FilterType::Supplementary);
			logfile().list("Supplementary reads: keep. (parameter 'keepSupplementaryReads')");
		} else {
			enable(FilterType::Supplementary, "SupplementaryAlignment");
			logfile().list("Supplementary reads: filter out. (use 'keepSupplementaryReads' to keep)");
		}

		// fragment length
		if (parameters().exists("filterReadsLongerThanFragment")) {
			enable(FilterType::LongerThanFragment, "Longer than fragment");
			logfile().list("Reads longer than fragment size: filter out. (parameter 'filterReadsLongerThanFragment')");
		} else {
			disable(FilterType::LongerThanFragment);
			logfile().list(
				"Reads longer than fragment size: keep. (use 'filterReadsLongerThanFragment' to filter out)");
		}

		// strand
		if (parameters().exists("keepOnlyFwd")) {
			disable(FilterType::FwdStrand);
			enable(FilterType::RevStrand, "Reverse strand");
			logfile().list("Strand: keep only forward. (parameter 'keepOnlyFwd')");
		} else if (parameters().exists("keepOnlyRev")) {
			enable(FilterType::FwdStrand, "Forward strand");
			disable(FilterType::RevStrand);
			logfile().list("Strand: keep only reverse. (parameter 'keepOnlyRev')");
		} else {
			disable(FilterType::FwdStrand);
			disable(FilterType::RevStrand);
			logfile().list("Strand: keep forward and reverse. (use 'keepOnlyFwd' or 'keepOnlyRev' to limit)");
		}

		// mate
		if (parameters().exists("keepOnlyFirst")) {
			enable(FilterType::FirstMate, "Second mate");
			disable(FilterType::SecondMate);
			logfile().list("Mate: keep only first. (parameter 'keepOnlyFirst')");
		} else if (parameters().exists("keepOnlySecond")) {
			disable(FilterType::FirstMate);
			enable(FilterType::SecondMate, "First mate");
			logfile().list("Mate: keep only second. (parameter 'keepOnlySecond')");
		} else {
			disable(FilterType::FirstMate);
			disable(FilterType::SecondMate);
			logfile().list("Mate: keep first and second. (use 'keepOnlyFirst' or 'keepOnlySecond' to limit)");
		}

		// blacklist
		if (parameters().exists("blacklist")) {
			std::string blacklistFilename = parameters().get("blacklist");
			logfile().list("Will filter out reads present in the file '" + blacklistFilename +
						   "'. (parameter 'blacklist')");
			_blacklist.addFromFile(blacklistFilename);
			enable(FilterType::Blacklist, "Was in provided blacklist");
		} else {
			disable(FilterType::Blacklist);
			logfile().list("Blacklist: keep all. (use 'blacklist' to provide a list and filter specific reads)");
		}

		// Mapping quality filter
		if (parameters().exists("filterMQ")) {
			TNumericRange<size_t> Range;
			parameters().fill("filterMQ", Range);

			enable(FilterType::MappingQuality, Range, "MappingQualityOutside" + Range.rangeString());
			logfile().list("Mapping quality: restrict to range " + Range.rangeString() + ". (parameter 'filterMQ')");
		} else {
			disable(FilterType::MappingQuality);
			logfile().list("Mapping quality: keep all. (use 'filterMQ' to limit)");
		}

		// Read length filter
		if (parameters().exists("filterReadLength")) {
			TNumericRange<size_t> Range;
			parameters().fill("filterReadLength", Range);

			enable(FilterType::ReadLength, Range, "Read length outside " + Range.rangeString());
			logfile().list("Read length: restrict to range " + Range.rangeString() +
						   ". (parameter 'filterReadLength')");
		} else {
			disable(FilterType::ReadLength);
			logfile().list("Read length: keep all. (use 'filterReadLength' to limit)");
		}

		// Fragment length filter
		if (parameters().exists("filterFragmentLength")) {
			TNumericRange<size_t> Range;
			parameters().fill("filterFragmentLength", Range);

			enable(FilterType::FragmentLength, Range, "Fragment length outside " + Range.rangeString());
			logfile().list("Fragment length: restrict to range " + Range.rangeString() +
						   ". (parameter 'filterFragmentLength')");
		} else {
			disable(FilterType::FragmentLength);
			logfile().list("Fragment length: keep all. (use 'filterFragmentLength' to limit)");
		}
	}
	logfile().endIndent();

	// log filtered reads?
}

void TBamFilters::resize(size_t numRG, size_t numChrom, std::string_view Filename) {
	_numRG    = numRG;
	_numChrom = numChrom;

	if (parameters().exists("bamLog")) {
		std::string logFilename = parameters().get<std::string>("bamLog");
		if (logFilename.empty()) {
			logFilename = coretools::str::readBeforeLast(Filename, ".");
			logFilename += ".bamlog.txt.gz";
		}
		logfile().list("Will write all filtered out reads to '" + logFilename + "'.");
		_log.open(logFilename, 3);
	}
}

void TBamFilters::clone(const TBamFilters& Filters) {
	for (auto t = FilterType::min; t < FilterType::max; ++t) {
		_filters[t].clone(Filters[t], _numRG, _numChrom);
	}
	_ranges        = Filters._ranges;
	_softClipRatio = Filters._softClipRatio;
	_blacklist     = Filters._blacklist;
	_enabled       = Filters._enabled;
}

} // namespace BAM
