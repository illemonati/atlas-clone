#include "TReadTraverser.h"

#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringConversions.h"

namespace BAM {

namespace impl {
constexpr size_t million = 1'000'000;
std::string millionReads(size_t N) { return coretools::str::toStringWithPrecision((double) N / million, 1); }
}

using coretools::instances::logfile;
using coretools::instances::parameters;

TReadTraverser::TReadTraverser(bool EnableFilters)
    : TReadTraverser(parameters().get<std::string>("bam"), EnableFilters, 0) {}
TReadTraverser::TReadTraverser(std::string_view Name, bool EnableFilters, size_t i)
    : _bamFile(Name, i, EnableFilters), _rgInfo(_bamFile.readGroups()), _errorModels(_rgInfo) {
	if (parameters().exists("out")) {
		_outputName = parameters().get("out");
		if (i > 0) {
			_outputName += coretools::str::toString("_", i);
		} else {
			logfile().list("Writing output files with prefix '" + _outputName + "'. (parameter 'out')");
		}
	} else {
		// guess from BAM filename.
		_outputName = coretools::str::readBeforeLast(_bamFile.filename(), ".");
		logfile().list("Writing output files with prefix '" + _outputName + "'. (specify with 'out')");
	}
}

TReadTraverser::~TReadTraverser() {
	if (_rgInfo.isParsed()) _rgInfo.write(_outputName + "_RGInfo.json");
}

void TReadTraverser::_log() {
	if (_nextPrint == _noLog) return;

	if (_nextPrint == 0) {
		_timer.start();
		logfile().startIndent("Parsing through BAM file ", _bamFile.filename(), ":");
		logfile().list("Parsing Chromomsome ", chromosomes().front().name(), ".");
		_nextPrint = impl::million;
	}

	if (bamFile().numAlignments() > _nextPrint) {
		logfile().list("Parsed ", impl::millionReads(bamFile().numAlignments()), " million reads (est. ",
		               coretools::str::toStringWithPrecision(bamFile().filePercentage(), 2) + "%) in ",
		               _timer.formattedTime());
		_nextPrint += impl::million;
	}

	if (_eor) {
		logfile().list("Reached end of BAM file ", bamFile().filename(), " in " + _timer.formattedTime() + ':');
		logfile().conclude("Parsed a total of ", impl::millionReads(bamFile().numAlignments()),  " million reads in " + _timer.formattedTime() + '.');
		logfile().endIndent();
		bamFile().printSummary(_outputName);
	} else if (bamFile().chrChanged()) {
		logfile().list("Parsing Chromomsome ", curChr().name(), ".");
	}
}

void TReadTraverser::nextRead() {
	_eor = !bamFile().readNextAlignmentThatPassesFilters();
	_log();
}

bool TReadTraverser::endOfReads() {
	// only at first call, so filters and other settings can be set after constructor
	if (bamFile().atStart()) { nextRead(); }
	return _eor;
}
} // namespace BAM
