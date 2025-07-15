#include "TPileup2Bed.h"
#include "coretools/Files/TInputRcpp.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/stringManipulations.h"
#include "genometools/GenomePositions/TGenomePosition.h"

namespace GenomeTasks {
using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::str::readBeforeLast;

TPileup2Bed::TPileup2Bed()
	: _pileup(parameters().get("pileup"), coretools::FileType::Header),
	  _bedName(parameters().get("out", readBeforeLast(_pileup.name(), ".txt")) + ".bed"),
	  _depthRange(parameters().get("depth")) {
	logfile().list("Opening pileup file '", _pileup.name(), "'.");
	logfile().list("Keeping depths in  ", _depthRange.rangeString(), ".");
}

void TPileup2Bed::run() {
	genometools::TBed bed;
	const auto ids = _pileup.indices({"chromosome", "position", "depth"});

	for (; !_pileup.empty(); _pileup.popFront()) {
		const auto depth = _pileup.get<size_t>(ids.back());
		if (_depthRange.within(depth)) {
			bed.add(_pileup.get(ids[0]), _pileup.get<size_t>(ids[1]) -1);
		}
	}
	logfile().list("Writing bed file ", _bedName, ".");
	bed.write(_bedName);
}
}
