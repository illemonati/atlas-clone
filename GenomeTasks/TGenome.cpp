/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */

#include "TGenome.h"

#include <cstdint>
#include <exception>
#include <iostream>
#include <map>
#include <math.h>
#include <memory>

#include "genometools/BED/TBed.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TSiteSubset.h"
#include "coretools/Math/TSubsamplePicker.h"
#include "coretools/Strings/stringFunctions.h"

namespace coretools {
class TRandomGenerator;
}

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;
using namespace coretools::str;

//---------------------------------------------------------------
// TGenome_basic
// A base class without filters and genotype likelihoods
//---------------------------------------------------------------

TGenome_basic::TGenome_basic() {
	// open bam file
	// TODO: deal with index in better way: let tasks decide if they need an index or not
	_bamFile.open(parameters().getParameter<std::string>("bam"), parameters().parameterExists("indexNotRequired"),
				  &logfile());
	_bamFile.setLimits(parameters(), &logfile());

	// outputname
	if (parameters().parameterExists("out")) {
		_outputName = parameters().getParameter<std::string>("out");
		logfile().list("Writing output files with prefix '" + _outputName + "'. (parameter 'out')");
	} else {
		// guess from BAM filename.
		_outputName = _bamFile.filename();
		_outputName = extractBeforeLast(_outputName, ".");
		logfile().list("Writing output files with prefix '" + _outputName + "'. (specify with 'out')");
	}
};

void TGenome_basic::_openBamForWriting(const std::string &filename, BAM::TOutputBamFile &outBam) {
	outBam.open(filename, _bamFile);
};

//---------------------------------------------------------------
// TGenome_filtered
// A base class without genotype likelihoods but BAM filters enabled
//---------------------------------------------------------------

TGenome_filtered::TGenome_filtered() { _bamFile.setFilters(parameters(), &logfile()); };

void TGenome_filtered::_traverseBAMPassedQC() {
	// parse through bam file
	_bamFile.startProgressReporting();
	while (_bamFile.readNextAlignmentThatPassesFilters()) {
		// handle alignment by derived classes
		_handleAlignment();

		// report
		_bamFile.printProgress();
	}

	// report
	_bamFile.printEndWithSummary(_outputName);
};

//---------------------------------------------------------------
// TGenome_parsed
// A base class with BAM filters and recalibration
//---------------------------------------------------------------
TGenome_parsed::TGenome_parsed() : _genotypeLikelihoodCalculator(&_bamFile.readGroupsMutable()) {
	// set parsing filters
	_setReadTrimming();
	_qualityFilter.set(parameters(), &logfile());
	_contextFilter.set(parameters(), &logfile());
	_bamFile.setFilters(parameters(), &logfile());
};

void TGenome_parsed::_openReference(bool required) {
	if (!_reference.hasReference()) {
		if (parameters().parameterExists("fasta")) {
			std::string fastaFile = parameters().getParameter<std::string>("fasta");
			int bufferSize = parameters().getParameterWithDefault<int>("fastaBuffer", 100'000);
			logfile().list("Reading reference sequence from '" + fastaFile + "'. (parameter fasta)");
			_reference.initialize(fastaFile, bufferSize);
		} else {
			if (required) { throw "No reference provided! (Use parameter fasta to provide a reference)"; }
		}
	}
	_bamFile.setFilters(parameters(), &logfile());
};

void TGenome_parsed::_setReadTrimming() {
	// trimming ends
	if (parameters().parameterExists("trim3") || parameters().parameterExists("trim5")) {
		_trimmingLength3Prime = parameters().getParameterWithDefault<int>("trim3", 0);
		if (_trimmingLength3Prime < 0) throw "trimming distance trim3 must be >= 0!";
		_trimmingLength5Prime = parameters().getParameterWithDefault<int>("trim5", 0);
		if (_trimmingLength5Prime < 0) throw "trimming distance trim5 must be >= 0!";
		if (_trimmingLength3Prime > 0 || _trimmingLength5Prime > 0) {
			logfile().list("Will trim first " + toString(_trimmingLength3Prime) + " and " +
						   toString(_trimmingLength5Prime) +
						   " bases from the 3' and 5' end, respectively. (parameters 'trim3', 'trim5')");
		}
		_trimReads = true;
	} else {
		_trimmingLength3Prime = 0;
		_trimmingLength5Prime = 0;
		_trimReads            = false;
	}
};

void TGenome_parsed::_parseAlignment(BAM::TAlignment &alignment) {
	// parse
	alignment.parse(_genotypeLikelihoodCalculator.getSequencingErrorModels());

	// apply filters
	if (_trimReads) { alignment.trimRead(_trimmingLength3Prime, _trimmingLength5Prime); }

	alignment.filter(_qualityFilter);
	alignment.filter(_contextFilter);

	if (_reference.hasReference()) { alignment.addReference(_reference); }
};

void TGenome_parsed::_traverseBAMPassedQC() {
	// parse through bam file
	_bamFile.startProgressReporting();
	while (_bamFile.readNextAlignmentThatPassesFilters(_alignment)) {
		// parse
		_parseAlignment(_alignment);

		// handle alignment by derived classes
		_handleAlignment();

		// report
		_bamFile.printProgress();
	}

	// report
	_bamFile.printEndWithSummary(_outputName);
};

//---------------------------------------------------------------
// TGenome_windows
// A base class to traverse a BAM file in windows
//---------------------------------------------------------------
TGenome_windows::TGenome_windows() : _chromosomes(_bamFile.chromosomes()) {
	// reading parameters regarding windows
	logfile().startIndent("Parsing window settings:");
	_setWindowParameters();
	_setParsingLimits();
	_setWindowFilters();
	_setMasks();
	_setSiteFilters();
	logfile().endIndent();
};

void TGenome_windows::_setWindowParameters() {
	if (!parameters().parameterExists("window") && parameters().parameterExists("windows")) {
		logfile().warning("Argument 'windows' specified, but unknown. Did you mean 'window'?");
	}
	std::string tmp = parameters().getParameterWithDefault<std::string>("window", "1000000");

	// check if it is a number
	if (stringIsProbablyANumber(tmp)) {
		_windowSize = convertString<int>(tmp);
		logfile().list("Setting window size to " + toString(_windowSize) + ". (parameter 'window')");
		if (_windowSize < _bamFile.maxReadLength()) {
			throw "Window size " + tmp + " out of range! Windows must be at least as large as the max read length (" +
				toString(_bamFile.maxReadLength()) + " bp). (use parameter 'filterReadLength' to change)!";
		}
	} else {
		logfile().listFlush("Limiting analysis to windows defined in BED file '" + tmp + "' (parameter window) ...");
		_predefinedWindows.add(tmp, _chromosomes);
		logfile().done();
		logfile().conclude("Read " + toString(_predefinedWindows.size()) + " of cumulative length " +
						   toString(_predefinedWindows.length()) + " bp on " +
						   toString(_predefinedWindows.numChromosomesWithWindows()) + " chromosomes.");
	}
	_numWindowsOnChr = 0;
};

void TGenome_windows::_setParsingLimits() {
	// limit windows
	_skipWindows = parameters().getParameterWithDefault<int>("skipWindows", 0);
	if (_skipWindows > 0)
		logfile().list("Will skip the first " + toString(_skipWindows) +
					   " windows per chromosome. (parameter 'skipWindows')");
	_limitWindows = parameters().getParameterWithDefault<long>("limitWindows", 1000000000);
	if (parameters().parameterExists("limitWindows"))
		logfile().list("Will limit analysis to the first " + toString(_limitWindows) +
					   " windows per chromosome. (parameter 'limitWindows')");
	if (_limitWindows <= _skipWindows) throw "limitWindows has to be larger than skipWindows!";
};

void TGenome_windows::_setWindowFilters() {
	// filter for missing reference
	_maxMissing = parameters().getParameterWithDefault<double>("maxMissing", 1.0);
	if (_maxMissing < 0.0 || _maxMissing > 1.0) throw "maxMissing must be within [0, 1]!";
	if (_maxMissing < 1.0) {
		logfile().list("Will filter out windows with a missing data fraction > " + toString(_maxMissing) +
					   ". (parameter 'maxMissing')");
	} else {
		logfile().list("Will keep windows regardless of missingness. (use 'maxMissing' to filter)");
	}

	_maxRefN = parameters().getParameterWithDefault<double>("maxRefN", 1.0);
	if (_maxRefN < 0.0 || _maxRefN > 1.0) throw "maxRefN must be within interval [0,1]!";
	_openReference();
	if (_maxRefN < 1.0 && !_reference.hasReference())
		throw "Can only calculate percentage of reference bases that are 'N' in window if reference file is provided! "
			  "(use 'fasta' to provide a reference)";
	logfile().list("Will filter out windows with a fraction of 'N' in reference > " + toString(_maxRefN) +
				   ". (parameter 'maxRefN')");
};

void TGenome_windows::_setSiteFilters() {
	// depth filter
	_readUpToDepth = parameters().getParameterWithDefault<uint32_t>("readUpToDepth", 1000);
	logfile().list("Will read data up to depth " + toString(_readUpToDepth) +
				   " and ignore additional bases. (parameter 'readUpToDepth')");

	// depth filter
	if (parameters().parameterExists("filterDepth")) {
		parameters().fillParameter("filterDepth", _depthFilter);
		_applyDepthFilter = true;
		logfile().list("Will filter out sites with sequencing depth outside ", _depthFilter,
					   ". (parameters 'filterDepth')");
	} else {
		_applyDepthFilter = false;
		logfile().list("Will keep sites regardless of depth. (use 'filterDepth' to filter)");
	}

	// downsample?
	_downsampleDepth = parameters().getParameterWithDefault<int>("downsample", 0);
	if (_downsampleDepth > 0) {
		logfile().list("Will downsample sites to a depth <= ", _downsampleDepth, ". (parameter 'downsample')");
		if (_depthFilter.larger(_downsampleDepth)) {
			logfile().warning("Downsample depth is >= max of depth filter: no downsampling will occur.");
		}
		subsamplePicker = std::make_unique<coretools::TSubsamplePicker>(30);
	}

	// CpG filter
	if (parameters().parameterExists("filterCpG")) {
		_filterCpG = true;
		logfile().list("Will filter out CpG sites. (parameter 'filterCpG')");
		_openReference(true);
	} else {
		_filterCpG = false;
		logfile().list("Will keep CpG sites. (use 'filterCpG' to remove)");
	}
};

void TGenome_windows::_setMasks() {
	// normal mask
	if (parameters().parameterExists("mask") || parameters().parameterExists("regions")) {
		std::string filename;
		if (parameters().parameterExists("mask")) {
			// mask
			if (parameters().parameterExists("regions")) throw "Cannot use mask and regions at the same time.";
			filename = parameters().getParameter<std::string>("mask");
			logfile().startIndent("Will mask all sites listed in BED file '" + filename + "':");
			_doMasking       = true;
			_considerRegions = false;
		} else {
			// regions
			filename = parameters().getParameter<std::string>("regions");
			logfile().startIndent("Will limit analysis to sites listed in BED file '" + filename +
								  "' (parameter 'regions'):");
			_doMasking       = false;
			_considerRegions = true;
		}

		// read file
		logfile().listFlush("Reading file ...");
		_mask.add(filename, _chromosomes);
		logfile().done();
		logfile().conclude("Read " + toString(_mask.size()) + " sites on " +
						   toString(_mask.numChromosomesWithWindows()) + " chromosomes.");
		logfile().endIndent();
	} else {
		_doMasking       = false;
		_considerRegions = false;
	}
};

void TGenome_windows::_openSiteSubset(const std::string &paramName) {
	if (_subset) { throw "Site subset already initialized!"; }

	std::string filename = parameters().getParameter<std::string>(paramName, true);

	if (_considerRegions)
		throw "Site subsets (parameter '" + paramName +
			"') and regions (parameter 'regions') can not be used at the same time!";
	if (_doMasking)
		throw "Site subsets (parameter '" + paramName +
			"') and masks (parameter 'mask') can not be used at the same time!";

	_subset = std::make_unique<GenotypeLikelihoods::TSiteSubset>(filename, _bamFile.chromosomes(), false, _reference);
};

void TGenome_windows::_setCountersBeginningOfChromosome() {
	_chrChangedWindow = true;
	_windowNumber     = 1;
};

bool TGenome_windows::_incrementWindow(GenotypeLikelihoods::TWindow_base &window) {
	// move to next
	window += _windowSize;
	++_windowNumber;

	// Move to next chromosome if 1) we are at begininning of BAM (_curchromosome at end), 2) we are beyond
	// _curChromosome or 3) reached window limit
	if (_curChromosome == _chromosomes.cend() || window.from() >= _curChromosome->chrEnd ||
		_windowNumber > _limitWindows) {
		// move to next chromosome
		if (_curChromosome == _chromosomes.cend()) { // beginning of chromosome
			_curChromosome = _chromosomes.cbegin();
		} else {
			++_curChromosome;
		}

		// advance if this chromosome is not used
		while (_curChromosome != _chromosomes.cend() &&
			   (!_curChromosome->inUse || _skipWindows * _windowSize > _curChromosome->length)) {
			++_curChromosome;
		}

		if (_curChromosome == _chromosomes.cend()) { return false; }

		_setCountersBeginningOfChromosome();
		_numWindowsOnChr = ceil(_curChromosome->length / (double)_windowSize);

		// move window to beginning of chromosome
		genometools::TGenomePosition newFrom = _curChromosome->chrStart + _skipWindows * _windowSize;
		window.move(newFrom, _windowSize, _curChromosome->name);
	}

	return true;
};

bool TGenome_windows::_moveToNextWindow(GenotypeLikelihoods::TWindow_base &window) {
	// move to next
	if (!_incrementWindow(window)) { return false; }

	// if sites defined, check if there are sites
	if (_subset) {
		while (!_subset->hasPositionsInWindow(window)) {
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
	if (window.to() > _curChromosome->chrEnd) { window.resize(_curChromosome->chrEnd - window.from()); }

	return true;
};

bool TGenome_windows::_incrementPredefinedWindow() {
	uint32_t oldRefID = _curPredefinedWindow->refID();

	++_curPredefinedWindow;
	if (_curPredefinedWindow == _predefinedWindows.end()) { return false; }

	// check if chromosome changed
	if (_curPredefinedWindow->refID() != oldRefID) {
		_setCountersBeginningOfChromosome();
	} else {
		++_windowNumber;
	}

	return true;
};

bool TGenome_windows::_moveToNextPredefinedWindow(GenotypeLikelihoods::TWindow_base &window) {
	// if at beginning of BAM file: restart
	if (_curChromosome == _chromosomes.cend()) {
		_curPredefinedWindow = _predefinedWindows.begin();
		_curChromosome       = _chromosomes.cbegin(_curPredefinedWindow->refID());
		_setCountersBeginningOfChromosome();
	} else {
		_incrementPredefinedWindow();
	}

	// go to next accepted window
	while (_curPredefinedWindow != _predefinedWindows.end() &&
		   (!_chromosomes.inUse(_curPredefinedWindow->refID()) || _windowNumber < _skipWindows ||
			_windowNumber >= _limitWindows)) {
		_incrementPredefinedWindow();
	}

	// return false if at end
	if (_curPredefinedWindow == _predefinedWindows.end()) { return false; }

	// make sure we are on the right chromosome
	if (_curChromosome->refID() != _curPredefinedWindow->refID()) {
		_curChromosome = _chromosomes.cbegin(_curPredefinedWindow->refID());

		// update num windows per chromosome
		_numWindowsOnChr = _predefinedWindows.numWindowsOnChr(_curChromosome->refID());
	}

	// move coordinates of window
	window.move(*_curPredefinedWindow, _curChromosome->name);

	// move in BAM: should we jump or are we already close enough to next window
	if (_bamFile.refID() == window.refID()) {
		// same chromosome: jump only if we are far away
		if (_bamFile.curPosition() > window.from() ||
			_bamFile.curPosition() < window.from() - _bamFile.maxReadLength()) {
			_bamFile.jump(window.from() - _bamFile.maxReadLength());
		}
	} else {
		// different chromosome: jump
		_bamFile.jump(window.from() - _bamFile.maxReadLength());
	}

	// return true as we continue reading
	return true;
};

bool TGenome_windows::_moveWindow(GenotypeLikelihoods::TWindow_base &window) {
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
		logfile().startNumbering("Parsing chromosome '" + _curChromosome->name + "':");
	}

	// report window
	logfile().number("Window [", window.from().position() + 1, ", ", window.to().position(), "] of ", _numWindowsOnChr,
					 " on '", _curChromosome->name, "':");
	logfile().addIndent();
	_hasWindowIndent = true;

	return true;
};

//---------------------
// read data in windows
//---------------------
bool TGenome_windows::_readDataInNextWindow(GenotypeLikelihoods::TWindow &window) {
	_windowTimer.start();

	// move window
	if (!_moveWindow(window)) {
		// reached end
		if (_hasWindowIndent) {
			logfile().removeIndent();
			_hasWindowIndent = false;
		}
		return false;
	}

	// read data
	_readAlignmentsIntoWindow(window);
	return true;
};

bool TGenome_windows::_readAndParseAlignment() {
	if (!_bamFile.readNextAlignmentThatPassesFilters(_curAlignment)) { return false; }

	_parseAlignment(_curAlignment);
	return true;
};

void TGenome_windows::_readAlignmentsIntoWindow(GenotypeLikelihoods::TWindow &window) {
	// measure runtime
	logfile().listFlushTime("Reading data ...");

	// use last read from last window
	if (_curAlignment.isEmpty() && !_readAndParseAlignment()) return;

	do {
		if (_curAlignment >= window.to()) break; 
		if (_curAlignment.lastAlignedPositionWithRespectToRef() >= window.from()) {
			window.addAlignment(_curAlignment);
		}
	} while (_readAndParseAlignment());
	// _curAlignment now holds first alignment of next window, don't discard!

	// fill sites
	if (_subset) {
		window.fillSitesSubset(*_subset, _readUpToDepth);
		window.addReferenceBaseToSites(*_subset);
	} else {
		window.fillSites(_readUpToDepth);
		window.addReferenceBaseToSites(_reference);
	}

	// report
	logfile().doneTime();

	// apply filters
	_applyWindowFilters(window);
};

void TGenome_windows::_applyWindowFilters(GenotypeLikelihoods::TWindow_base &window) {
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
		if (_filterCpG) { window.maskCpG(_reference); }
		if (_downsampleDepth > 0) { window.downsample(_downsampleDepth, *subsamplePicker); };
	}

	// apply filters on window
	window.filter(_maxMissing, _maxRefN, &logfile());

	// report
	if (window.numReadsInWindow() > 0) {
		window.dataSummary();
	} else {
		logfile().conclude("No data in this window.");
	}
};

void TGenome_windows::_traverseBAMWindows() {
	logfile().startIndent("Traversing BAM file in windows:");

	// initializing
	_hasWindowIndent = false;
	_curChromosome   = _chromosomes.cend(); // set chromosome to end to trigger restart.
	_curAlignment.clear();

	// iterate through windows
	while (_readDataInNextWindow(_window)) {
		if (_window.passedFilters()) {
			// do stuff in derived classes
			_handleWindow();

			// report end of window calculations
			logfile().list("Total computation time for this window was ", _windowTimer.formattedTime(), ".");
		}
		_chrChangedWindow = false;
	}

	logfile().endIndent();

	// report
	_bamFile.printEndWithSummary(_outputName);
};

}; // namespace GenomeTasks
