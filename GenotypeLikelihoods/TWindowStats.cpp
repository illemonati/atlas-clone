#include "TWindowStats.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {

using coretools::instances::parameters;

void TWindowStats::_traverseWindows() {
	for (; !_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		for (; !_siteTraverser.endOfCurChr(); _siteTraverser.nextWindow()) {
			const auto &window = _siteTraverser.window();
			_out.writeln(_siteTraverser.curChr().name(), window.from().position(), window.to().position(), 1,
						 window.depth(), window.size(), window.numSites(), window.numSitesWithData());
		}
	}
}

void TWindowStats::_traverseGenome() {
	for (; !_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		for (; !_siteTraverser.endOfCurChr(); _siteTraverser.nextSite()) {
		}
	}
	_out.writeln(1, _siteTraverser.depth(), _siteTraverser.chromosomes().totLength(), _siteTraverser.numSites(), _siteTraverser.numSitesWithData());
}

TWindowStats::TWindowStats() {
	_genomeWide = parameters().exists("genomeWide");
}

void TWindowStats::run() {
	const auto fName = _siteTraverser.outputName() + "_winStats.txt.gz";
	coretools::instances::logfile().list("Writing windows summary statistics to file '", fName, "'.");

	if (_genomeWide) {
		_out.open(_siteTraverser.outputName() + "_winStats.txt.gz",
				  {"Prob", "depth", "NSitesTot", "NSitesUsed", "NSitesWithData"});
		_traverseGenome();
	} else {
		_out.open(_siteTraverser.outputName() + "_winStats.txt.gz",
				  {"Chr", "Start", "End", "Prob", "depth", "NSitesTot", "NSitesUsed", "NSitesWithData"});
		_traverseWindows();
	}
}
}
