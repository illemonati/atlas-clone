#include "TReadTraverser.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {

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
} // namespace GenomeTasks
