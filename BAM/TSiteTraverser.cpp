#include "TSiteTraverser.h"

#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/stringProperties.h"
#include "genometools/GenomePositions/TGenomePosition.h"

#include <filesystem>

namespace BAM{
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
	if (_i >= _bamWindow.size()) nextWindow();
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
	if (_windowList) {
		if (RefID >= chromosomes().size()) {
			_window.move(RefID, 0, _wSize); // end of current chromosomes
		} else if (!_windowList.windowsOnChr(RefID).empty()) {
			_window   = _windowList.windowsOnChr(RefID).front();
			_iWindows = 0;
		} else {
			_window.move(chromosomes()[RefID].to(), 1); // move outside of chr
		}
	} else {
		_window.move(RefID, 0, _wSize);
	}
	if (!endOfChrs()) {
		_skipShinkFill();
	}
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
}

void TSiteTraverser::requireReference() const {
	coretools::user_assert(_alnTraverser.reference().isOpen(),
						   "No reference provided! (Use parameter fasta to provide a reference)");
}

bool TSiteTraverser::endOfChrs() {
	constexpr auto zerozero = genometools::TGenomePosition{}; 
	if (_window.size() == 0 && _window.from() == zerozero) {
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
}

void TSiteTraverser::nextChr() {
	_initChr(refID() + 1);
}

}
