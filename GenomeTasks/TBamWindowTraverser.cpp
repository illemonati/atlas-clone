/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */

#include "TBamWindowTraverser.h"

#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TSubsamplePicker.h"
#include "coretools/Strings/stringFunctions.h"
#include <filesystem>

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;

	TBamWindowTraverser::TBamWindowTraverser() : _parser(_genome) {

	const BAM::TBamFilters filters{true};
	_genome.bamFile().setFilters(filters);
	// reading parameters regarding windows
	logfile().startIndent("Parsing window settings:");
	_setWindowParameters();
	_setParsingLimits();
	_setWindowFilters();
	_setMasks();
	_setSiteFilters();
	logfile().endIndent();
	_numWindowsOnChr = 0;
}

void TBamWindowTraverser::_setWindowParameters() {
	const auto sWindow = parameters().get<std::string>("window", "1000000");

	// check if it is a number
	if (std::filesystem::exists(sWindow)) {
		logfile().listFlush("Limiting analysis to windows defined in BED file '" + sWindow +
							"' (parameter window) ...");
		_predefinedWindows.add(sWindow, _genome.bamFile().chromosomes());
		logfile().done();
		logfile().conclude("Read ", _predefinedWindows.size(), " of cumulative length ", _predefinedWindows.length(),
						   " bp on ", _predefinedWindows.numChromosomesWithWindows(), " chromosomes.");
	} else {
		coretools::str::fromString(sWindow, _windowSize);
		logfile().list("Setting window size to ", _windowSize, ". (parameter 'window')");
	}
}

void TBamWindowTraverser::_setParsingLimits() {
	// limit windows
	_skipWindows = parameters().get<int>("skipWindows", 0);
	if (_skipWindows > 0)
		logfile().list("Will skip the first ", _skipWindows, " windows per chromosome. (parameter 'skipWindows')");
	_limitWindows = parameters().get<long>("limitWindows", 1000000000);
	if (parameters().exists("limitWindows"))
		logfile().list("Will limit analysis to the first ", _limitWindows,
					   " windows per chromosome. (parameter 'limitWindows')");
	if (_limitWindows <= _skipWindows) UERROR("limitWindows has to be larger than skipWindows!");
}

void TBamWindowTraverser::_setWindowFilters() {
	// filter for missing reference
	_maxMissing = parameters().get<double>("maxMissing", 1.0);
	if (_maxMissing < 0.0 || _maxMissing > 1.0) UERROR("maxMissing must be within [0, 1]!");
	if (_maxMissing < 1.0) {
		logfile().list("Will filter out windows with a missing data fraction > ", _maxMissing,
					   ". (parameter 'maxMissing')");
	} else {
		logfile().list("Will keep windows regardless of missingness. (use 'maxMissing' to filter)");
	}

	_maxRefN = parameters().get<double>("maxRefN", 1.0);
	if (_maxRefN < 0.0 || _maxRefN > 1.0) UERROR("maxRefN must be within interval [0,1]!");
	_parser.openReference();
	if (_maxRefN < 1.0 && !_parser.reference().isOpen())
		UERROR("Can only calculate percentage of reference bases that are 'N' in window if reference file is provided! "
			   "(use 'fasta' to provide a reference)");
	logfile().list("Will filter out windows with a fraction of 'N' in reference > ", _maxRefN,
				   ". (parameter 'maxRefN')");
}

void TBamWindowTraverser::_setSiteFilters() {
	// depth filter
	_readUpToDepth = parameters().get<size_t>("readUpToDepth", 1000);
	logfile().list("Will read data up to depth ", _readUpToDepth,
				   " and ignore additional bases. (parameter 'readUpToDepth')");

	// depth filter
	if (parameters().exists("filterDepth")) {
		parameters().fill("filterDepth", _depthFilter);
		_applyDepthFilter = true;
		logfile().list("Will filter out sites with sequencing depth outside ", _depthFilter,
					   ". (parameters 'filterDepth')");
	} else {
		_applyDepthFilter = false;
		logfile().list("Will keep sites regardless of depth. (use 'filterDepth' to filter)");
	}

	// downsample?
	_downsampleDepth = parameters().get<int>("downsample", 0);
	if (_downsampleDepth > 0) {
		logfile().list("Will downsample sites to a depth <= ", _downsampleDepth, ". (parameter 'downsample')");
		if (_depthFilter.larger(_downsampleDepth)) {
			logfile().warning("Downsample depth is >= max of depth filter: no downsampling will occur.");
		}
		_subsamplePicker = std::make_unique<coretools::TSubsamplePicker>(30);
	}

	// CpG filter
	if (parameters().exists("filterCpG")) {
		_filterCpG = true;
		logfile().list("Will filter out CpG sites. (parameter 'filterCpG')");
		_parser.openReference(true);
	} else {
		_filterCpG = false;
		logfile().list("Will keep CpG sites. (use 'filterCpG' to remove)");
	}
}

void TBamWindowTraverser::_setMasks() {
	// normal mask
	if (parameters().exists("mask") || parameters().exists("regions")) {
		std::string filename;
		if (parameters().exists("mask")) {
			// mask
			if (parameters().exists("regions")) UERROR("Cannot use mask and regions at the same time.");
			filename = parameters().get<std::string>("mask");
			logfile().startIndent("Will mask all sites listed in BED file '" + filename + "':");
			_doMasking       = true;
			_considerRegions = false;
		} else {
			// regions
			filename = parameters().get<std::string>("regions");
			logfile().startIndent("Will limit analysis to sites listed in BED file '" + filename +
								  "' (parameter 'regions'):");
			_doMasking       = false;
			_considerRegions = true;
		}

		// read file
		logfile().listFlush("Reading file ...");
		_mask.add(filename, _genome.bamFile().chromosomes());
		logfile().done();
		logfile().conclude("Read ", _mask.size(), " sites on ", _mask.numChromosomesWithWindows(), " chromosomes.");
		logfile().endIndent();
	} else {
		_doMasking       = false;
		_considerRegions = false;
	}
}

void TBamWindowTraverser::_openSiteSubset(const std::string &paramName, bool polymorphic) {
	//report
	if(polymorphic){
		logfile().startIndent("Limiting analysis to sites with known alleles (parameter '", paramName, "'):");
	} else {
		logfile().startIndent("Limiting analysis to sites with known allele (parameter '", paramName, "'):");
	}
	
	// only allow for one subset to be active
	if (_subsetPolymoprhic || _subsetMonomorphic) { DEVERROR("Site subset already initialized!"); }

	if (_considerRegions)
		UERROR("Site subsets (parameter '", paramName,
			   "') and regions (parameter 'regions') can not be used at the same time!");
	if (_doMasking)
		UERROR("Site subsets (parameter '", paramName,
			   "') and masks (parameter 'mask') can not be used at the same time!");

	const auto filename = parameters().get(paramName);

	if(polymorphic){
		_subsetPolymoprhic = std::make_unique<GenotypeLikelihoods::TSiteSubsetPolymorphic>(filename, _genome.bamFile().chromosomes());
	} else {
		_subsetMonomorphic = std::make_unique<GenotypeLikelihoods::TSiteSubsetMonomorphic>(filename, _genome.bamFile().chromosomes());
	}	
	logfile().endIndent();
}

void TBamWindowTraverser::_setCountersBeginningOfChromosome() {
	_chrChangedWindow = true;
	_windowNumber     = 1;
}

bool TBamWindowTraverser::_incrementWindow(GenotypeLikelihoods::TWindow &window) {
	// move to next
	window += _windowSize;
	++_windowNumber;

	// Move to next chromosome if 1) we are at begininning of BAM (_curchromosome at end), 2) we are beyond
	// _curChromosome or 3) reached window limit
	if (_curChromosome == _genome.bamFile().chromosomes().cend() || window.from() >= _curChromosome->to() ||
		_windowNumber > _limitWindows) {
		// move to next chromosome
		if (_curChromosome == _genome.bamFile().chromosomes().cend()) { // beginning of chromosome
			_curChromosome = _genome.bamFile().chromosomes().cbegin();
		} else {
			++_curChromosome;
		}

		// advance if this chromosome is not used
		while (_curChromosome != _genome.bamFile().chromosomes().cend() &&
			   (!_curChromosome->inUse() || _skipWindows * _windowSize > _curChromosome->length())) {
			++_curChromosome;
		}

		if (_curChromosome == _genome.bamFile().chromosomes().cend()) { return false; }

		_setCountersBeginningOfChromosome();
		_numWindowsOnChr = ceil(_curChromosome->length() / (double)_windowSize);

		// move window to beginning of chromosome
		const auto newFrom = _curChromosome->from() + _skipWindows * _windowSize;
		window.move(newFrom, _windowSize, _curChromosome->name());
	}

	return true;
}

bool TBamWindowTraverser::_moveToNextWindow(GenotypeLikelihoods::TWindow &window) {
	// move to next
	if (!_incrementWindow(window)) { return false; }

	// if sites defined, check if there are sites
	if (_subsetPolymoprhic) {
		while (!_subsetPolymoprhic->hasPositionsInWindow(window)) {
			if (!_incrementWindow(window)) { return false; }
		}
	} else if(_subsetMonomorphic){
		while (!_subsetMonomorphic->hasPositionsInWindow(window)) {
			if (!_incrementWindow(window)) { return false; }
		}
	}

	if (_considerRegions) {
		if (_mask.windowIsBeyond(window)) { return false; }
		while (!_mask.hasOverlapWith(window)) {
			if (!_incrementWindow(window)) { return false; }
		}
	}

	// make sure window does not go beyond chromosome end
	if (window.to() > _curChromosome->to()) { window.resize(_curChromosome->to() - window.from()); }

	return true;
}

bool TBamWindowTraverser::_incrementPredefinedWindow() {
	const auto oldRefID = _curPredefinedWindow->refID();

	++_curPredefinedWindow;
	if (_curPredefinedWindow == _predefinedWindows.end()) { return false; }

	// check if chromosome changed
	if (_curPredefinedWindow->refID() != oldRefID) {
		_setCountersBeginningOfChromosome();
	} else {
		++_windowNumber;
	}

	return true;
}

bool TBamWindowTraverser::_moveToNextPredefinedWindow(GenotypeLikelihoods::TWindow &window) {
	// if at beginning of BAM file: restart
	if (_curChromosome == _genome.bamFile().chromosomes().cend()) {
		_curPredefinedWindow = _predefinedWindows.begin();
		_curChromosome       = _genome.bamFile().chromosomes().cbegin(_curPredefinedWindow->refID());
		_setCountersBeginningOfChromosome();
	} else {
		_incrementPredefinedWindow();
	}

	// go to next accepted window
	while (_curPredefinedWindow != _predefinedWindows.end() &&
		   (!_genome.bamFile().chromosomes().inUse(_curPredefinedWindow->refID()) || _windowNumber < _skipWindows ||
			_windowNumber >= _limitWindows)) {
		_incrementPredefinedWindow();
	}

	// return false if at end
	if (_curPredefinedWindow == _predefinedWindows.end()) { return false; }

	// make sure we are on the right chromosome
	if (_curChromosome->refID() != _curPredefinedWindow->refID()) {
		_curChromosome = _genome.bamFile().chromosomes().cbegin(_curPredefinedWindow->refID());

		// update num windows per chromosome
		_numWindowsOnChr = _predefinedWindows.numWindowsOnChr(_curChromosome->refID());
	}

	// move coordinates of window
	window.move(*_curPredefinedWindow, _curChromosome->name());

	// move in BAM: should we jump or are we already close enough to next window
	if (_genome.bamFile().refID() == window.refID()) {
		// same chromosome: jump only if we are far away
		if (_genome.bamFile().curPosition() > window.from() ||
			_genome.bamFile().curPosition() < window.from() - _genome.bamFile().maxReadLength()) {
			_genome.bamFile().jump(window.from() - _genome.bamFile().maxReadLength());
		}
	} else {
		// different chromosome: jump
		_genome.bamFile().jump(window.from() - _genome.bamFile().maxReadLength());
	}

	// return true as we continue reading
	return true;
}

bool TBamWindowTraverser::_moveWindow(GenotypeLikelihoods::TWindow &window) {
	// returns false when end of genome is reached
	if (_predefinedWindows.empty()) {
		// no predefined windows: regular traversing
		if (!_moveToNextWindow(window)) { return false; }
	} else {
		// there are predefined windows
		if (!_moveToNextPredefinedWindow(window)) { return false; }
	}

	// report chromosome
	if (_hasWindowIndent) logfile().removeIndent();
	if (_chrChangedWindow) {
		if (_curChromosome->refID() > 0) { logfile().endIndent(); }
		logfile().startNumbering("Parsing chromosome '" + _curChromosome->name() + "':");
	}

	// report window
	logfile().number("Window [", window.from().position() + 1, ", ", window.to().position(), "] of ", _numWindowsOnChr,
					 " on '", _curChromosome->name(), "':");
	logfile().addIndent();
	_hasWindowIndent = true;

	return true;
}

//---------------------
// read data in windows
//---------------------

bool TBamWindowTraverser::_readAndParseAlignment(BAM::TAlignment &_curAlignment) {
	if (!_genome.bamFile().readNextAlignmentThatPassesFilters()) return false;

	_genome.bamFile().fill(_curAlignment);
	_parser.apply(_curAlignment);
	return true;
}

void TBamWindowTraverser::_readAlignmentsIntoWindow(GenotypeLikelihoods::TWindow &window) {
	// measure runtime
	logfile().listFlushTime("Reading data ...");

	// use last read from last window
	BAM::TAlignment _curAlignment;
	_genome.bamFile().fill(_curAlignment);
	_parser.apply(_curAlignment);

	do {
		if (_curAlignment >= window.to()) break;
		if (_curAlignment.lastAlignedPositionWithRespectToRef() >= window.from()) {
			window.addAlignment(_curAlignment);
		}
	} while (_readAndParseAlignment(_curAlignment));
	// _curAlignment now holds first alignment of next window, don't discard!

	// fill sites
	if (_subsetPolymoprhic) {
		window.fillSitesSubset(*_subsetPolymoprhic, _readUpToDepth);
		window.addReferenceBaseToSites(*_subsetPolymoprhic);
	} else if (_subsetMonomorphic) {
		window.fillSitesSubset(*_subsetMonomorphic, _readUpToDepth);
		window.addReferenceBaseToSites(*_subsetMonomorphic);
	} else {
		window.fillSites(_readUpToDepth);
		window.addReferenceBaseToSites(_parser.reference());
	}

	// report
	logfile().doneTime();

	// apply filters
	_applyWindowFilters(window);
}

void TBamWindowTraverser::_applyWindowFilters(GenotypeLikelihoods::TWindow &window) {
	// apply site-specific filters
	if (window.numReadsInWindow() > 0) {
		// apply masks and filters
		if (_doMasking) {
			logfile().listFlush("Masking sites ...");
			window.applyMask(_mask, _considerRegions);
			logfile().done();
		} else if (_considerRegions) {
			logfile().listFlush("Masking sites outside regions ...");
			window.applyMask(_mask, _considerRegions);
			logfile().done();
		}

		// filter sites
		if (_applyDepthFilter) { window.applyDepthFilter(_depthFilter); }
		if (_filterCpG) { window.maskCpG(_parser.reference()); }
		if (_downsampleDepth > 0) { window.downsample(_downsampleDepth, *_subsamplePicker); };
	}

	// apply filters on window
	window.filter(_maxMissing, _maxRefN);

	// report
	if (window.numReadsInWindow() > 0) {
		window.dataSummary();
	} else {
		logfile().conclude("No data in this window.");
	}
}

void TBamWindowTraverser::_traverseBAMWindows() {
	logfile().startIndent("Traversing BAM file in windows:");

	// initializing
	_hasWindowIndent = false;
	_curChromosome   = _genome.bamFile().chromosomes().cend(); // set chromosome to end to trigger restart.

	// iterate through windows
	GenotypeLikelihoods::TWindow window;
	_genome.bamFile().readNextAlignmentThatPassesFilters();
	while (_moveWindow(window)) {
		coretools::TTimer _windowTimer;
		_readAlignmentsIntoWindow(window);
		if (window.passedFilters()) {
			// do stuff in derived classes
			_handleWindow(window);

			// report end of window calculations
			logfile().list("Total computation time for this window was ", _windowTimer.formattedTime(), ".");
		}
		_chrChangedWindow = false;
	}

	if (_hasWindowIndent) {
		logfile().removeIndent();
		_hasWindowIndent = false;
	}

	logfile().endIndent();

	// report
	_genome.bamFile().printEndWithSummary(_genome.outputName());
}

} // namespace GenomeTasks
