#include "TSiteTraverser.h"
#include "coretools/Main/TParameters.h"

namespace BAM{
using coretools::instances::parameters;

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
	_bamWindow.move(_window);

	for (; !_alnTraverser.endOfAlignments(); _alnTraverser.nextAlignment()) {
		const auto& aln = _alnTraverser.alignment();
		if (aln.from() > window().to()) break;
		if (aln.to() < window().from()) continue;
		
		_bamWindow.add(aln);
	}

	_findFirstI();
	if (_i >= _bamWindow.size()) nextWindow();
}

void TSiteTraverser::_findFirstI() {
	for (_i = 0; _i < _bamWindow.size(); ++_i) {
		if (_bamWindow[_i].depth() >= 1) return;
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

TSiteTraverser::TSiteTraverser()
	: _bamWindow(chromosomes()), _windowList(parameters().get("windows", ""), chromosomes(), false) {}

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
	do { ++_i; } while (_i < _bamWindow.size() && _bamWindow[_i].depth() < 1);

	if (_i >= _bamWindow.size()) {
		nextWindow();
	}
}

void TSiteTraverser::nextChr() {
	_initChr(refID() + 1);
}

}
