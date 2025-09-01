#include "TReadTraverser.h"

#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringConversions.h"

namespace BAM {

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

void TReadTraverser::nextRead() {
	constexpr size_t million = 1'000'000;
	static size_t _nextPrint = million;

	if (bamFile().numAlignments() > _nextPrint) {
		logfile().list("Parsed ", bamFile().numAlignments() / million, " million reads (est. ",
					   coretools::str::toStringWithPrecision(bamFile().filePercentage(), 2) + "%) in ",
					   _timer.formattedTime());
		_nextPrint += million;
	}
	const auto refID = curChr().refID();
	_eor = !bamFile().readNextAlignmentThatPassesFilters();
	if (_eor) {
		logfile().list("Reached end of BAM file ", bamFile().filename(), " in " + _timer.formattedTime() + ':');
		logfile().conclude("Parsed a total of ", bamFile().numAlignments()/million, " million reads in " + _timer.formattedTime() + '.');
		logfile().endIndent();
		bamFile().printSummary(_outputName);
	} else if (refID != curChr().refID()) {
		logfile().list("Parsing Chromomsome ", curChr().name(), ".");
	}
}

bool TReadTraverser::endOfReads() {
	if (bamFile().atStart()) {
		_timer.start();
		logfile().startIndent("Parsing through BAM file ", _bamFile.filename() , ":");
		logfile().list("Parsing Chromomsome ", chromosomes().front().name(), ".");
		nextRead();
	}
	return _eor;
}
} // namespace BAM
