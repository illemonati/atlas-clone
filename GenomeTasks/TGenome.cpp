/*
 * TGenome.cpp
 */

#include "TGenome.h"

#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Strings/toString.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;

TGenome::TGenome(const BAM::TBamFilters &Filters) : TGenome(parameters().get<std::string>("bam"), Filters, 0) {}

TGenome::TGenome(std::string_view Name, const BAM::TBamFilters &Filters, size_t i)
    : _bamFile(Name, i), _rgInfo(_bamFile.readGroups()), _errorModels(_rgInfo) {
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
	_bamFile.setFilters(Filters);
}

TGenome::~TGenome() {
	if (_rgInfo.isParsed()) _rgInfo.write(_outputName + "_RGInfo.json");
}

}; // namespace GenomeTasks
