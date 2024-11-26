#include "TBamWindows.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;

TBamWindows::TBamWindows(const genometools::TChromosomes& Chromosomes)  {
	_parser.openReference(); // non-mandatory

	logfile().startIndent("Parsing window settings:");
	_setWindowFilters();
	_setMasks(Chromosomes);
	_setSiteFilters();
	_setWindowParameters(Chromosomes);
	logfile().endIndent();
}

void TBamWindows::requireReference() const {
	if (!_parser.reference().isOpen())  UERROR("No reference provided! (Use parameter fasta to provide a reference)");
}

void TBamWindows::_setWindowParameters(const genometools::TChromosomes& Chromosomes) {
	const auto sWindow = parameters().get<std::string>("window", "1000000");
	size_t lTot        = 0;
	size_t nTot        = 0;
	size_t nUsed       = 0;

	_windows.resize(Chromosomes.size());

	if (std::filesystem::exists(sWindow)) {
		logfile().listFlush("Reading windows defined in BED file '", sWindow, "' (parameter window) ...");

		std::vector<genometools::TGenomeWindow> windows;

		for (coretools::TInputFile iFile(sWindow, coretools::FileType::NoHeader); !iFile.empty(); iFile.popFront()) {
			const auto refId = Chromosomes.refID(iFile.get(0));
			const auto start = iFile.get<size_t>(1);
			const auto end   = iFile.get<size_t>(2);
			windows.emplace_back(refId, start, end - start);
		}
		std::sort(windows.begin(), windows.end());

		for (auto &window: windows) {
			const auto& chr = Chromosomes[window.refID()];
			if (!chr.inUse() || (_subsetPolymoprhic && !_subsetPolymoprhic->hasPositionsInWindow(window)) ||
				(_subsetMonomorphic && !_subsetMonomorphic->hasPositionsInWindow(window)) ||
				(_considerRegions && !_mask.overlaps(window))) {
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

		// limit windows
		const auto skip = parameters().get<size_t>("skipWindows", 0);
		if (skip > 0)
			logfile().list("Will skip the first ", skip, " windows per chromosome. (parameter 'skipWindows')");
		const auto limit = parameters().get<size_t>("limitWindows", 1000000000);
		if (parameters().exists("limitWindows"))
			logfile().list("Will limit analysis to the first ", limit,
						" windows per chromosome. (parameter 'limitWindows')");
		if (limit <= skip) UERROR("limitWindows has to be larger than skipWindows!");

		for (const auto& chr: Chromosomes) {
			if (!chr.inUse()) continue;

			++nUsed;
			const genometools::TGenomePosition from(chr.refID(), _windowSize * skip);
			const genometools::TGenomePosition to(chr.refID(),
												  std::min(chr.to().position(), limit * _windowSize));

			for (genometools::TGenomeWindow window(from, _windowSize); window.from() < to; window += _windowSize) {
				if ((_subsetPolymoprhic && !_subsetPolymoprhic->hasPositionsInWindow(window)) ||
					(_subsetMonomorphic && !_subsetMonomorphic->hasPositionsInWindow(window)) ||
					(_considerRegions && !_mask.overlaps(window))) {
					continue;
				}
				if (window.to() > chr.to()) window.resize(chr.to() - window.from());

				_windows[chr.refID()].push_back(window);
				lTot += window.size();
				++nTot;
			}
		}
	}
	if (parameters().exists("shuffleSites")) {
		logfile().list("Will shuffle bases at sites. (parameter 'shuffleSites')");
		_shuffleSites = true;
	} else {
		logfile().list("Will not shuffle bases at sites. (use 'shuffleSites')");
		_shuffleSites = false;
	}
	logfile().conclude("Will traverse ", nTot, " windows with cumulative length of ", lTot,
						   " bp on ", nUsed, " chromosomes.");
}

void TBamWindows::_setWindowFilters() {
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
	if (_maxRefN < 1.0 && !_parser.reference().isOpen())
		UERROR("Can only calculate percentage of reference bases that are 'N' in window if reference file is provided! "
			   "(use 'fasta' to provide a reference)");
	logfile().list("Will filter out windows with a fraction of 'N' in reference > ", _maxRefN,
				   ". (parameter 'maxRefN')");
}

void TBamWindows::_setSiteFilters() {
	// depth filter
	_upToDepth = parameters().get<size_t>("readUpToDepth", 1000);
	logfile().list("Will read data up to depth ", _upToDepth,
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
		requireReference();
	} else {
		_filterCpG = false;
		logfile().list("Will keep CpG sites. (use 'filterCpG' to remove)");
	}
}

void TBamWindows::_setMasks(const genometools::TChromosomes& Chromosomes) {
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
		_mask.parse(filename, Chromosomes);
		logfile().done();
		logfile().conclude("Read ", _mask.size(), " sites on ", _mask.NChrWindows(), " chromosomes.");
		logfile().endIndent();
	} else {
		_doMasking       = false;
		_considerRegions = false;
	}
}

void TBamWindows::openSiteSubset(const std::string &paramName, const genometools::TChromosomes& Chromosomes, bool polymorphic) {
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
		_subsetPolymoprhic = std::make_unique<GenotypeLikelihoods::TSiteSubsetPolymorphic>(filename, Chromosomes);
	} else {
		_subsetMonomorphic = std::make_unique<GenotypeLikelihoods::TSiteSubsetMonomorphic>(filename, Chromosomes);
	}	
	logfile().endIndent();
}

void TBamWindows::filter(GenotypeLikelihoods::TWindow &Window) {
	// apply site-specific filters
	if (Window.numReadsInWindow() > 0) {
		// apply masks and filters
		if (_doMasking) {
			const auto N = Window.applyMask(_mask, _considerRegions);
			logfile().list("Masking ", N, " sites.");
		} else if (_considerRegions) {
			const auto N = Window.applyMask(_mask, _considerRegions);
			logfile().list("Masking ", N, " sites outside regions.");
		}

		// filter sites
		if (_applyDepthFilter) { Window.applyDepthFilter(_depthFilter); }
		if (_filterCpG) { Window.maskCpG(_parser.reference()); }
		if (_downsampleDepth > 0) { Window.downsample(_downsampleDepth, *_subsamplePicker); };
	}

	// apply filters on window
	Window.filter(_maxMissing, _maxRefN);

	// report
	if (Window.numReadsInWindow() > 0) {
		Window.dataSummary();
	} else {
		logfile().conclude("No data in this window.");
	}
}

void TBamWindows::fillSites(GenotypeLikelihoods::TWindow &Window) {
	// fill sites
	if (_subsetPolymoprhic) {
		Window.fillSitesSubset(*_subsetPolymoprhic, _upToDepth);
		Window.addReferenceBaseToSites(*_subsetPolymoprhic);
	} else if (_subsetMonomorphic) {
		Window.fillSitesSubset(*_subsetMonomorphic, _upToDepth);
		Window.addReferenceBaseToSites(*_subsetMonomorphic);
	} else {
		Window.fillSites(_upToDepth);
		Window.addReferenceBaseToSites(_parser.reference());
	}
	if (_shuffleSites) {
		for (auto& site: Window) {
			site.shuffle();
		}
	}
}



}
