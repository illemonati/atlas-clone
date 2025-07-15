#include "TBamWindows.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/toString.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "genometools/TAlleles.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::user_assert;

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
	user_assert(_parser.reference().isOpen(), "No reference provided! (Use parameter fasta to provide a reference)");
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

		for (coretools::TInputFile bedFile(sWindow, coretools::FileType::NoHeader); !bedFile.empty(); bedFile.popFront()) {
			const auto refId = Chromosomes.refID(bedFile.get(0));
			const auto start = bedFile.get<size_t>(1);
			const auto end   = bedFile.get<size_t>(2);
			windows.emplace_back(refId, start, end - start);
		}
		std::sort(windows.begin(), windows.end());

		for (auto &window: windows) {
			const auto& chr = Chromosomes[window.refID()];
			if (!chr.inUse() || (_alleles && !_alleles.overlaps(window)) ||
				(considerRegions() && !_regions.overlaps(window))) {
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
		user_assert(limit > skip, "limitWindows has to be larger than skipWindows!");

		for (const auto& chr: Chromosomes) {
			if (!chr.inUse()) continue;

			++nUsed;
			const genometools::TGenomePosition from(chr.refID(), _windowSize * skip);
			const genometools::TGenomePosition to(chr.refID(),
												  std::min(chr.to().position(), limit * _windowSize));

			for (genometools::TGenomeWindow window(from, _windowSize); window.from() < to; window += _windowSize) {
				if ((_alleles && !_alleles.overlaps(window)) ||
					(considerRegions() && !_regions.overlaps(window))) {
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
	user_assert(_maxMissing >= 0.0 && _maxMissing <= 1.0, "maxMissing must be within [0, 1]!");
	if (_maxMissing < 1.0) {
		logfile().list("Will filter out windows with a missing data fraction > ", _maxMissing,
					   ". (parameter 'maxMissing')");
	} else {
		logfile().list("Will keep windows regardless of missingness. (use 'maxMissing' to filter)");
	}

	_maxRefN = parameters().get<double>("maxRefN", 1.0);
	user_assert(_maxRefN >= 0.0 && _maxRefN <= 1.0, "maxRefN must be within interval [0,1]!");
	user_assert(_maxRefN == 1.0 || _parser.reference().isOpen(), "Can only calculate percentage of reference bases that are 'N' in window if reference file is provided! "
			   "(use 'fasta' to provide a reference)");
	logfile().list("Will filter out windows with a fraction of 'N' in reference > ", _maxRefN,
				   ". (parameter 'maxRefN')");
}

void TBamWindows::_setSiteFilters() {
	// depth filter
	_upToDepth = parameters().get<size_t>("readUpToDepth", 1000);
	logfile().list("Will read data up to depth ", _upToDepth,
				   " and ignore additional bases. (parameter 'readUpToDepth')");

	constexpr std::string_view downsample = "downsampleSites";
	_downProb = parameters().get(downsample, coretools::P(0.));
	if (_downProb > 0.) {
		logfile().list("Will downsample sites with probability ", _downProb, ".(parameter '", downsample, "')");
	} else {
		logfile().list("Will not downsample sites.(use '", downsample, "')");
	}

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

void TBamWindows::_setMasks(const genometools::TChromosomes &Chromosomes) {
	// normal mask
	if (parameters().exists("regions")) { // regions
		const auto regionsFile = parameters().get("regions");

		logfile().list("Will limit analysis to sites listed in BED file '", regionsFile, "' (parameter 'regions'):");
		_regions.parse(regionsFile, Chromosomes);
	}
	if (parameters().exists("mask")) {
		const auto maskFile = parameters().get("mask");
		logfile().list("Will mask all sites listed in BED file '", maskFile, "':");
		if (!parameters().exists("regions")) {
			_regions.parse(maskFile, Chromosomes);

			// flip mask to get regions
			_regions.flip(Chromosomes);
		} else {
			genometools::TBed mask;
			mask.parse(maskFile, Chromosomes);
			_regions.addMask(mask);
			auto fName = coretools::str::toString(coretools::str::readBeforeLast(parameters().get("regions"), '.'),
												  "_masked.bed");
			logfile().list("Will write masked regions-file to '", fName, "'.");
			_regions.write(fName);
		}
	}
}

void TBamWindows::openSiteSubset(const std::string &paramName, const genometools::TChromosomes& Chromosomes, genometools::Morphic Morph) {
	DEV_ASSERT(_alleles.empty());
	//report
	if(Morph == genometools::Morphic::Poly || Morph == genometools::Morphic::Both){
		logfile().startIndent("Limiting analysis to sites with known alleles (parameter '", paramName, "'):");
	} else {
		logfile().startIndent("Limiting analysis to sites with known allele (parameter '", paramName, "'):");
	}
	
	// only allow for one subset to be active

	user_assert(!considerRegions(), "Site subsets (parameter '", paramName,
				"') and regions (parameter 'regions' or 'mask') can not be used at the same time!");

	const auto filename = parameters().get(paramName);
	_alleles.parse(filename, Chromosomes, Morph);

	logfile().endIndent();
}

void TBamWindows::filter(GenotypeLikelihoods::TWindow &Window) {
	// apply site-specific filters
	if (Window.numReadsInWindow() > 0) {
		// filter sites
		if (_applyDepthFilter) { Window.applyDepthFilter(_depthFilter); }
		if (_filterCpG) { Window.maskCpG(_parser.reference()); }
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
	if (_alleles) Window.limitSites(_alleles);
	Window.addReferenceBaseToSites(_parser.reference());
	if (_downProb != 0.) Window.downsampleSites(_downProb);
	Window.limitDepth(_upToDepth, _shuffleSites);
}

}
