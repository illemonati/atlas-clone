/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */

#include "TBamWindowTraverser.h"

#include "TWindow.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TSubsamplePicker.h"
#include "coretools/Strings/stringFunctions.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenomePositions/TGenomeWindow.h"
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

void TBamWindowTraverser::_fillAlignments(GenotypeLikelihoods::TWindow &window) {
	// measure runtime
	logfile().listFlushTime("Reading data ...");

	// first, use last read from last window, before reading next
	do {
		const auto curPos = _genome.bamFile().curPosition();
		const auto maxLen = _genome.bamFile().curCIGAR().lengthSequenced();
		
		if (curPos >= window.to()) break; // too far
		if (curPos.position() + maxLen < window.fromOnChr()) continue; // too short

		BAM::TAlignment _curAlignment;
		_genome.bamFile().fill(_curAlignment);
		_parser.apply(_curAlignment);

		if (_curAlignment.lastAlignedPositionWithRespectToRef() >= window.from()) {
			window.addAlignment(_curAlignment);
		}
	} while (_genome.bamFile().readNextAlignmentThatPassesFilters());

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
	using GenotypeLikelihoods::TWindow;

	_genome.bamFile().startProgressReporting();
	logfile().startIndent("Traversing BAM file in windows:");
	_genome.bamFile().readNextAlignmentThatPassesFilters();
	for (const auto &chr : _genome.bamFile().chromosomes()) {
		if (!chr.inUse()) continue;
		logfile().startNumbering("Parsing chromosome '" + chr.name() + "':");
		_onChrChange(chr);

		const size_t nWindows = std::ceil((double)chr.to().position() / _windowSize);
		const genometools::TGenomePosition to(chr.refID(), std::min(chr.to().position(), _limitWindows*_windowSize));

		for (TWindow window(chr.name(), chr.refID(), _windowSize * _skipWindows, _windowSize); window.from() < to;
			 window += _windowSize) {
			if (window.to() > chr.to()) window.resize(chr.to() - window.from());
			if ((_subsetPolymoprhic && !_subsetPolymoprhic->hasPositionsInWindow(window)) ||
			    (_subsetMonomorphic && !_subsetMonomorphic->hasPositionsInWindow(window)) ||
			    (_considerRegions && !_mask.hasOverlapWith(window))) {
				continue;
			}

			logfile().number("Window [", window.from().position() + 1, ", ", window.to().position(), "] of ",
							 nWindows, " on '", chr.name(), "':");

			coretools::TTimer _windowTimer;
			_fillAlignments(window);
			if (window.passedFilters()) {
				// do stuff in derived classes
				_handleWindow(window);

				// report end of window calculations
				logfile().list("Total computation time for this window was ", _windowTimer.formattedTime(), ".");
			}
		}
		logfile().endIndent();
	}
	logfile().endIndent();

	// report
	_genome.bamFile().printEndWithSummary(_genome.outputName());
}

} // namespace GenomeTasks
