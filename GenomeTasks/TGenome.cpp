/*
 * TGenome.cpp
 */

#include "TGenome.h"

#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;

TGenome::TGenome()
	: _bamFile(parameters().get<std::string>("bam")), _rgInfo(_bamFile.readGroups()) {
	// outputname
	if (parameters().exists("out")) {
		_outputName = parameters().get<std::string>("out");
		logfile().list("Writing output files with prefix '" + _outputName + "'. (parameter 'out')");
	} else {
		// guess from BAM filename.
		_outputName = coretools::str::readBeforeLast(_bamFile.filename(), ".");
		logfile().list("Writing output files with prefix '" + _outputName + "'. (specify with 'out')");
	}
}

TGenome::~TGenome() {
	if (_rgInfo.isParsed()) _rgInfo.write(_outputName + "_RGInfo.json");
}

}; // namespace GenomeTasks
