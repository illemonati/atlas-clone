#include "TSiteTraverser.h"

#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/stringConversions.h"
#include "coretools/Strings/stringProperties.h"
#include "coretools/devtools.h"
#include "genometools/GenomePositions/TGenomePosition.h"

#include <filesystem>

namespace BAM{

namespace impl {
constexpr size_t million = 1'000'000;
std::string smillions(size_t N) { return coretools::str::toStringWithPrecision((double) N / million, 1); }

std::string bamNames(const std::vector<TAlignmentTraverser>& traversers) {
	DEV_ASSERT(!traversers.empty());

	std::string bamNames = "BAM file";
	if (traversers.size() == 1) {
		bamNames += " " + traversers.front().bamFile().filename();
	} else {
		bamNames += "s ";
		for (const auto &traverser : traversers) {
			bamNames += traverser.bamFile().filename() + ",";
		}
		// get rid of last ","
		bamNames.pop_back();
	}
	return bamNames;
}

std::vector<TAlignmentTraverser> initAln() {
	const auto bams   = coretools::instances::parameters().get<std::vector<std::string>>("bam");
	std::vector<TAlignmentTraverser> vec;
	vec.reserve(bams.size());
	vec.emplace_back(bams.front(), 0);
	for (size_t i = 1; i < bams.size(); ++i) {
		vec.emplace_back(bams[i], i, vec.front().bamFile().filters());
	}
	return vec;
}
} // namespace impl

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
	_bamWindow.move(_window, _traverser().reference());

	for (auto &traverser : _alnTraversers) {
		for (; !traverser.endOfAlignments(); traverser.nextAlignment()) {
			const auto &aln = traverser.alignment();
			if (aln.from() > _window.to()) break;
			if (aln.to() < _window.from()) continue;

			_bamWindow.add(aln);
		}
	}
	_bamWindow.filter();

	if (_bamWindow.endOfSites()) {
		nextWindow();
	}
	_winChanged      = true;
	_numSitesChr    += _bamWindow.numSites();
	_numWithDataChr += _bamWindow.numSitesWithData();
	_numBasesChr    += _bamWindow.numBases();
	_numReadsChr    += _bamWindow.numReads();
}

void TSiteTraverser::_initChr(size_t RefID) {
	_numSitesTot +=_numSitesChr;
	_numSitesChr = 0;

	_numWithDataTot += _numWithDataChr;
	_numWithDataChr = 0;

	_numBasesTot +=_numBasesChr;
	_numBasesChr = 0;

	_numReadsTot += _numReadsChr;
	_numReadsChr = 0;

	if (RefID >= chromosomes().size()) {
		// end of chromosomes
		_window.move(RefID, 0, _wSize);
		return;
	}

	logfile().startIndent("Parsing chromosome '", chromosomes()[RefID].name(), "'");

	if (!chromosomes()[RefID].inUse() || (regions() && regions().empty(RefID)) ||
	           (alleles() && alleles().empty(RefID)) || (_windowList && _windowList.empty(RefID))) {
		// skip chromosome
		_iWindows = -1;
		_window.move(chromosomes()[RefID].to(), 1);
		logfile().list("Chromsome is not used.");
		return;
	}

	_iWindows = 0;
	if (_windowList) {
		// take from window list
		_window   = _windowList[RefID].front();
	} else {
		// use window size
		_window.move(RefID, 0, _wSize);
	}
	_skipShinkFill();
}

    TSiteTraverser::TSiteTraverser(genometools::Morphic Mo) : _alnTraversers(impl::initAln()), _bamWindow(chromosomes(), Mo) {
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
		logfile().list("Using default window-size: ", _wSize, ". (set with 'window')");
	}
	_traverser().openReference(_bamWindow.filtersCpG());
	_traverser().setSilent();
}

void TSiteTraverser::filterRefN(bool Yes) noexcept {
	DEV_ASSERT(_atStart());
	DEV_ASSERT(_traverser().reference().isOpen());
	logfile().list("Will filter sites where Reference = N.");
	_bamWindow.filterRefN(Yes);
}

void TSiteTraverser::requireReference() const {
	DEV_ASSERT(_atStart());
	coretools::user_assert(_traverser().reference().isOpen(),
						   "No reference provided! (Use parameter fasta to provide a reference)");
}

void TSiteTraverser::requireSingleBAM() const {
	DEV_ASSERT(_atStart());
	coretools::user_assert(numSamples() == 1, "Cannot handle more than one BAM-file not allowed!");
}

bool TSiteTraverser::endOfChrs() {
	if (_atStart()) {
		DEBUG_ASSERT(_window.from() == genometools::TGenomePosition{});

		logfile().startIndent("Traversing sites of ", impl::bamNames(_alnTraversers), ":");
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

TBamWindow &TSiteTraverser::window() noexcept {
	static genometools::TGenomePosition lastPos(-1, -1);
	if (lastPos != _window.from()) {
		lastPos = _window.from();
		_logWindows();
	}

	return _bamWindow;
}

void TSiteTraverser::nextWindow() {
	_advanceWindow();
	_skipShinkFill();
}

void TSiteTraverser::nextSite() {
	_winChanged = false;
	_bamWindow.nextSite();

	if (_bamWindow.endOfSites()) nextWindow();

	_logSites();
}

void TSiteTraverser::_logWindows() const {
	if (_windowList) {
		logfile().startIndent("Parsing Window ", _iWindows + 1, " [", _window.fromOnChr(), ", ", _window.toOnChr(),
							  "] of ", _windowList[refID()].size(), " on ", curChr().name(), ":");
	} else {
		const auto nWindows = (curChr().length() - 1)/_wSize + 1;
		const auto i        = _window.fromOnChr()/_wSize + 1;
		logfile().startIndent("Parsing Window ", i, " [", _window.fromOnChr(), ", ", _window.toOnChr(), "] of ",
							  nWindows, " on ", curChr().name());
	}
	logfile().list("Total number of sites: ", _bamWindow.size());
	logfile().list("Fraction of filtered sites: ", 1. - (double)_bamWindow.numSites()/_bamWindow.size());
	logfile().conclude("Number of used sites: ", _bamWindow.numSites());
	logfile().list("Number of used reads: ", _bamWindow.numReads());
	logfile().list("Number of used bases: ", _bamWindow.numBases());
	logfile().conclude("Average depth: ", _bamWindow.depth());
	logfile().list("Number of used sites without Data: ", _bamWindow.numSites() - _bamWindow.numSitesWithData());
	logfile().conclude("Fraction of used sites without Data: ", 1. - (double)_bamWindow.numSitesWithData()/_bamWindow.numSites());
	logfile().endIndent();
}

void TSiteTraverser::_logSites() const {
	static const size_t NtotSites = [this]() {
		size_t tot = 0;
		if (regions()) {
			for (const auto& win: regions()) {
				tot += win.size();
			}
		} else if (alleles()) {
			tot = alleles().size();
		} else {
			for (const auto &chr : chromosomes()) {
				if (chr.inUse()) tot += chr.length();
			}
		}
		return tot;
	}();
	static size_t counter = 0;
	static size_t nextPrint = impl::million;

	if (counter > nextPrint) {
		logfile().list("Parsed ", impl::smillions(counter), " million sites (est. ",
					   coretools::str::toStringWithPrecision(counter*100. / NtotSites, 2) + "%) in ",
					   _timer.formattedTime());
		nextPrint += impl::million;
	}
	++counter;
}

void TSiteTraverser::nextChr() {
	if (_iWindows != size_t(-1)) {
		logfile().startIndent("Properties of chromosome '", curChr().name(), "':");
		logfile().list("Total number of sites: ", curChr().length());
		logfile().list("Fraction of filtered sites: ", 1. - (double)_numSitesChr/curChr().length());
		logfile().conclude("Number of used sites: ", _numSitesChr);
		logfile().list("Number of used reads: ", _numReadsChr);
		logfile().list("Number of used bases: ", _numBasesChr);
		logfile().conclude("Average depth: ", double(_numBasesChr)/_numSitesChr);
		logfile().list("Number of used sites without Data: ", _numSitesChr - _numWithDataChr);
		logfile().conclude("Fraction of used sites without Data: ", 1. - (double)_numWithDataChr/_numSitesChr);
		logfile().endIndent();
	}
	logfile().endIndent(); // _initChr
	_initChr(refID() + 1);

	if (endOfChrs()) {
		const auto totLength = chromosomes().totLength();
		logfile().startIndent("Reached end of ", impl::bamNames(_alnTraversers), " in " + _timer.formattedTime() + ':');
		logfile().list("Total number of sites: ", totLength);
		logfile().list("Fraction of filtered sites: ", 1. - (double)_numSitesTot/totLength);
		logfile().conclude("Number of used sites: ", _numSitesTot);
		logfile().list("Number of used reads: ", _numReadsTot);
		logfile().list("Number of used bases: ", _numBasesTot);
		logfile().conclude("Average depth: ", double(_numBasesTot) / _numSitesTot);
		logfile().list("Number of used sites without Data: ", _numWithDataTot);
		logfile().conclude("Fraction of used sites without Data: ", double(_numSitesTot - _numWithDataTot) / _numSitesTot);
		logfile().endIndent();
		for (const auto& traverser: _alnTraversers) {
			traverser.bamFile().printSummary(traverser.outputName());
		}
		logfile().endIndent(); // from start in (endOfChrs)
	}
}

}
