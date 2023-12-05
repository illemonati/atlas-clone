/*
 * TGenome.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: wegmannd
 */

#include "TBamWindowTraverser.h"

#include "TWindow.h"
#include "coretools/Files/TFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TSubsamplePicker.h"
#include "coretools/Strings/stringFunctions.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/GenomePositions/TGenomeWindow.h"
#include <algorithm>
#include <filesystem>

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;

	TBamWindowTraverser::TBamWindowTraverser()  {
	_genome.bamFile().setFilters(BAM::TBamFilters{true});

	logfile().startIndent("Parsing window settings:");
	_setParsingLimits();
	_setWindowFilters();
	_setMasks();
	_setSiteFilters();
	_setWindowParameters();
	logfile().endIndent();
}

void TBamWindowTraverser::_setWindowParameters() {
	const auto sWindow = parameters().get<std::string>("window", "1000000");
	size_t lTot        = 0;
	size_t nTot        = 0;
	size_t nUsed       = 0;

	_windows.resize(_genome.bamFile().chromosomes().size());

	if (std::filesystem::exists(sWindow)) {
		logfile().listFlush("Reading windows defined in BED file '", sWindow, "' (parameter window) ...");
		coretools::TInputFile iFile(sWindow, 3);

		std::vector<genometools::TGenomeWindow> windows;
		std::vector<std::string> vec;

		while (iFile.read(vec)) {
			size_t refId = _genome.bamFile().chromosomes().refID(vec[0]);
			size_t start = coretools::str::fromString<size_t, true>(vec[1]);
			size_t end   = coretools::str::fromString<size_t, true>(vec[2]);
			windows.emplace_back(refId, start, end-start);
		}
		std::sort(windows.begin(), windows.end());

		for (auto &window: windows) {
			const auto& chr = _genome.bamFile().chromosomes()[window.refID()];
			if (!chr.inUse() || (_subsetPolymoprhic && !_subsetPolymoprhic->hasPositionsInWindow(window)) ||
				(_subsetMonomorphic && !_subsetMonomorphic->hasPositionsInWindow(window)) ||
				(_considerRegions && !_mask.hasOverlapWith(window))) {
				logfile().list("Ignoring window [", window.from().position(), ", ", window.to().position(), "] on chr ", chr.name(), "!");
				continue;
			}
			if (window.to() > chr.to()) {
				logfile().list("Resizing window [", window.from().position(), ", ", window.to().position(), "] on chr ", chr.name(), "!");
				window.resize(chr.to() - window.from());
			}

			_windows[chr.refID()].push_back(window);
			lTot += window.size();
			++nTot;
		}
		nUsed = std::count_if(_windows.begin(), _windows.end(), [](const auto& chr){return !chr.empty();});

		logfile().done();
	} else {
		coretools::str::fromString(sWindow, _windowSize);
		logfile().list("Setting window size to ", _windowSize, ". (parameter 'window')");

		for (const auto& chr: _genome.bamFile().chromosomes()) {
			if (!chr.inUse()) continue;

			++nUsed;
			const genometools::TGenomePosition from(chr.refID(), _windowSize * _skipWindows);
			const genometools::TGenomePosition to(chr.refID(),
												  std::min(chr.to().position(), _limitWindows * _windowSize));

			for (genometools::TGenomeWindow window(from, _windowSize); window.from() < to; window += _windowSize) {
				if ((_subsetPolymoprhic && !_subsetPolymoprhic->hasPositionsInWindow(window)) ||
					(_subsetMonomorphic && !_subsetMonomorphic->hasPositionsInWindow(window)) ||
					(_considerRegions && !_mask.hasOverlapWith(window))) {
					continue;
				}
				if (window.to() > chr.to()) window.resize(chr.to() - window.from());

				_windows[chr.refID()].push_back(window);
				lTot += window.size();
				++nTot;
			}
		}
	}
	logfile().conclude("Will traverse ", nTot, " windows with cumulative length of ", lTot,
						   " bp on ", nUsed, " chromosomes.");
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
	BAM::TAlignment alignment;

	// first, use last read from last window, before reading next
	do {
		const auto curPos = _genome.bamFile().curPosition();
		const auto maxLen = _genome.bamFile().curCIGAR().lengthSequenced();
		
		if (curPos >= window.to()) break; // too far
		if (curPos.position() + maxLen < window.fromOnChr()) continue; // too short

		_parser.fill(_genome, alignment);

		if (alignment.lastAlignedPositionWithRespectToRef() >= window.from()) {
			window.addAlignment(alignment);
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
		_handleChromosome(chr);
		TWindow window(chr.name());

		for (const auto& gWindow: _windows[chr.refID()]) {
			window.move(gWindow);
			logfile().number("Window [", window.from().position() + 1, ", ", window.to().position(), "] of ",
							 _windows[chr.refID()].size(), " on '", chr.name(), "':");

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
