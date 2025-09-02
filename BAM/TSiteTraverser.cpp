#include "TSiteTraverser.h"

#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/stringConversions.h"
#include "coretools/Strings/stringProperties.h"
#include "genometools/GenomePositions/TGenomePosition.h"

#include <filesystem>

namespace BAM{

namespace impl {
constexpr size_t millions = 1'000'000;
std::string millionsReads(size_t N) { return coretools::str::toStringWithPrecision((double) N / millions, 1); }
}

using coretools::instances::parameters;
using coretools::instances::logfile;

void TSiteTraverser::_skipShinkFill() {
	if (endOfCurChr()) return;

	// skip empty windows
	while ((regions() && !regions().overlaps(_window)) || (alleles() && !alleles().overlaps(_window))) {
		_advanceWindow();
		if (endOfCurChr()) return;
	}

	// shrink window if needed
	_window.move(_window.from(), std::min(_window.to(), curChr().to()));
	_fillWindow();
}

void TSiteTraverser::_fillWindow() {
	_bamWindow.move(_window, _alnTraverser.reference(), _filterCpG);

	for (; !_alnTraverser.endOfAlignments(); _alnTraverser.nextAlignment()) {
		const auto& aln = _alnTraverser.alignment();
		if (aln.from() > window().to()) break;
		if (aln.to() < window().from()) continue;
		
		_bamWindow.add(aln);
	}

	_filterFindI();

	if (_i >= _bamWindow.size()) {
		nextWindow();
	}
}

void TSiteTraverser::setDepthFilter(size_t Min, size_t Max) noexcept {
	_depthFilter.set(Min, true, Max, true);
	_filterFindI();
}

void TSiteTraverser::_filterFindI() {
	_i = -1;
	for (size_t i = 0; i < _bamWindow.size(); ++i) {
		if (_depthFilter.outside(_bamWindow[i].depth())) {
			_bamWindow.mask(i);
		} else if (_i == size_t(-1)) {
			_i = i;
		}
	}
}

void TSiteTraverser::_initChr(size_t RefID) {
	if (RefID >= chromosomes().size()) {
		// end of chromosomes
		_window.move(RefID, 0, _wSize);
		return;
	}

	logfile().startIndent("Parsing chromosome '", curChr().name(), "':");

	if (!chromosomes()[RefID].inUse() || (regions() && regions().empty(RefID)) ||
	           (alleles() && alleles().empty(RefID)) || (_windowList && _windowList.empty(RefID))) {
		// skip chromosome
		_window.move(chromosomes()[RefID].to(), 1);
		logfile().list("Chromsome is not used.");
		return;
	}

	if (_windowList) {
		// take from window list
		_window   = _windowList[RefID].front();
		_iWindows = 0;
	} else {
		// use window size
		_window.move(RefID, 0, _wSize);
	}
	_skipShinkFill();
}

    TSiteTraverser::TSiteTraverser() : _bamWindow(chromosomes()), _filterCpG(parameters().exists("filterCpG")) {
	const auto sWindows = parameters().get("window", "");
	if (!sWindows.empty()) {
		if (std::filesystem::exists(sWindows)) {
			logfile().list("Using user-defined window-list in file '", sWindows, "'. (parameter 'window')");
			_windowList.parse(sWindows, chromosomes(), false);
		} else if (coretools::str::stringIsProbablyANumber(sWindows)) {
			logfile().list("Setting user-defined window-size: ", sWindows, ". (parameter 'window')");
			coretools::str::fromString<true>(sWindows, _wSize);
			coretools::user_assert(_wSize > 0, "Window size must be > 0!");
		} else {
			logfile().list("Using user-defined window-list defined as follows '", sWindows, "'. (parameter 'window')");
			_windowList.parse(sWindows, chromosomes(), false);
		}
	} else {
		logfile().list("Using default window-size: ", sWindows, ". (set with 'window')");
	}

	if (_filterCpG) {
		logfile().list("Will filter out CpG sites. (parameter 'filterCpG')");
		_alnTraverser.openReference(true);
	} else {
		logfile().list("Will keep CpG sites. (use 'filterCpG' to remove)");
		_alnTraverser.openReference(false);
	}

	if (parameters().exists("filterDepth")) {
		parameters().fill("filterDepth", _depthFilter);
		logfile().list("Will filter out sites with sequencing depth outside ", _depthFilter,
					   ". (parameters 'filterDepth')");
	} else {
		logfile().list("Will keep all sites with data. (use 'filterDepth' to filter)");
	}
	_alnTraverser.setSilent();
}

void TSiteTraverser::requireReference() const {
	coretools::user_assert(_alnTraverser.reference().isOpen(),
						   "No reference provided! (Use parameter fasta to provide a reference)");
}

bool TSiteTraverser::endOfChrs() {
	if (_window.size() == 0) {
		DEBUG_ASSERT(_window.from() == genometools::TGenomePosition{});

		logfile().startIndent("Traversing through BAM file ", _alnTraverser.bamFile().filename(), " in windows:");
		_timer.start();
		_initChr(0);
	}
	return refID() >= chromosomes().size();
}

void TSiteTraverser::_advanceWindow() {
	if (_windowList) {
		++_iWindows;
		if (_windowList.windowsOnChr(refID()).size() > _iWindows) {
			_window = _windowList.windowsOnChr(refID())[_iWindows];
		} else {
			_window.move(curChr().to(), 1); // move outside of chr
		}
	} else {
		_window += _wSize;
	}
}

void TSiteTraverser::nextWindow() {
	_advanceWindow();
	_skipShinkFill();
}

void TSiteTraverser::nextSite() {
	do { ++_i; } while (_i < _bamWindow.size() && _depthFilter.outside(_bamWindow[_i].depth()));

	if (_i >= _bamWindow.size()) {
		nextWindow();
	}
	++_numSites;
	_log();
}

void TSiteTraverser::_log() {
	if (_numSites > _nextPrint) {
		logfile().list("Parsed ", impl::millionsReads(_numSites), " million sites (est. ",
					   coretools::str::toStringWithPrecision((double)_numSites / _window.size(), 2) + "%) in ",
					   _timer.formattedTime());
		_nextPrint += impl::millions;
	}
}

void TSiteTraverser::nextChr() {
	logfile().endIndent(); // from previous _initChr
	_initChr(refID() + 1);

	if (endOfChrs()) {
		logfile().list("Reached end of BAM file ", _alnTraverser.bamFile().filename(), " in " + _timer.formattedTime() + ':');
		logfile().conclude("Parsed a total of ", impl::millionsReads(_numSites),
						   " million sites in " + _timer.formattedTime() + '.');
		logfile().endIndent(); // from start in (endOfChrs)
		_alnTraverser.bamFile().printSummary(outputName());
	}
}

}
