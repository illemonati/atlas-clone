#include "TBamWindow.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TParameters.h"
#include <cstddef>

namespace BAM {

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::str::readBeforeLast;
using coretools::str::toString;

void TBamWindow::_makeBedOrAlleles(const genometools::TChromosomes &Chromosomes) {
	if (parameters().exists("regions")) { _regions.parse(parameters().get("regions"), Chromosomes); }
	if (parameters().exists("mask")) {
		if (_regions.empty()) {
			_regions.parse(parameters().get("mask"), Chromosomes);
			_regions.flip(Chromosomes);
		} else {
			genometools::TBed mask(parameters().get("mask"), Chromosomes);
			_regions.addMask(mask);
			auto fName = toString(readBeforeLast(parameters().get("regions"), '.'), "_masked.bed");
			logfile().list("Will write masked regions-file to '", fName, "'.");
			_regions.write(fName);
		}
	}
	if (parameters().exists("alleles")) {
		coretools::user_assert(_regions.empty(), "Cannot define regions and alleles at the same time!");
		_alleles.parse(parameters().get("alleles"), Chromosomes);
	}
}

TBamWindow::TBamWindow(const genometools::TChromosomes &Chromosomes)
	: _upToDepth(parameters().get<size_t>("readUpToDepth", 1000)), _filterCpG(parameters().exists("filterCpG")) {
	_makeBedOrAlleles(Chromosomes);

	logfile().list("Will read data up to depth ", _upToDepth,
	               " and ignore additional bases. (parameter 'readUpToDepth')");

	if (_filterCpG) {
		logfile().list("Will filter out CpG sites. (parameter 'filterCpG')");
	} else {
		logfile().list("Will keep CpG sites. (use 'filterCpG' to remove)");
	}

	if (parameters().exists("filterDepth")) {
		parameters().fill("filterDepth", _depthFilter);
		logfile().list("Will filter out sites with sequencing depth outside ", _depthFilter,
					   ". (parameters 'filterDepth')");
	} else {
		logfile().list("Will keep all sites with data. (use 'filterDepth' to filter)");
	}

	constexpr std::string_view downsample = "downsampleSites";
	_downProb = parameters().get(downsample, coretools::P(1.));
	if (_downProb < 1.) {
		logfile().list("Will downsample sites with probability ", _downProb, ".(parameter '", downsample, "')");
	} else {
		logfile().list("Will not downsample sites.(use '", downsample, "')");
	}
}

void TBamWindow::move(genometools::TGenomeWindow Window, const genometools::TFastaReader& Reference) {
	using genometools::Base;
	_iSite = 0;
	_from  = Window.from();

	_entries.assign(Window.size(), {});
	if (_tagReads) _readIDs.assign(Window.size(), {});
	// clear everything
	_numBases   = 0;
	_numReads   = 0;
	_sitesData  = 0;
	_sites2Plus = 0;
	_numMasked  = 0;

	if (_regions) {
		_masked.assign(Window.size(), true);
		_numMasked = _masked.size();

		for (auto it = _regions.begin(Window); it != _regions.end() && Window.overlaps(*it); ++it) {
			const auto to = std::min(Window.to(), it->to());
			for (auto pos = it->from(); pos < to; ++pos) {
				_masked[pos - Window.from()] = false;
				_numMasked--;
			}
		}
	} else if (_alleles) {
		if (Reference) logfile().warning("Allele-Files will overwrite ref-values given by reference!");
		_masked.assign(Window.size(), true);

		for (auto it = _alleles.begin(Window); it != _alleles.end() && Window.within(it->position); ++it) {
			const auto i = it->position - Window.from();
			_masked[i] = false;
			_entries[i].refBase = it->ref;
			_entries[i].altBase = it->alt;
			_numMasked--;
		}
	} else {
		_masked.assign(Window.size(), false);
	}

	if (Reference) {
		const auto view = Reference.view(window());
		for (size_t i = 0; i < size(); ++i) {
			_entries[i].refBase = view[i];
			if (_filterRefN && _entries[i].refBase == Base::N) {
				_masked[i] = true;
				++_numMasked;
			}
			if (_filterCpG && i > 0) {
				if (view[i - 1] == Base::C && view[i] == Base::G) {
					_masked[i - 1] = true; // C
					_masked[i]     = true; // G
					_numMasked += 2;
				}
			}
		}
	}

	if (_filterCpG) {
		// borders
		if (from().position() > 0 && Reference[from() - 1] == Base::C && Reference[from()] == Base::G) {
			_masked.front() = true;
			++_numMasked;
		}
		if (Reference[to() - 1] == Base::C && Reference[to()] == Base::G) {
			_masked.back() = true;
			++_numMasked;
		}
	}

	const auto copy = _overlap; // there could be alignments longer than window
	_overlap.clear();
	for (const auto& aln: copy) {
		add(aln);
	}
}

void TBamWindow::mask(size_t I) {
	_masked[I] = true;

	const auto d = _entries[I].depth();
	_entries[I].clear();

	// Housekeeping
	if (d > 1) _sites2Plus -= d;
	if (d > 0) _sitesData  -= d;
	_numBases -= d;
}

void TBamWindow::add(const TAlignment &Alignment) {
	for (size_t p = 0; p < Alignment.parsedLength(); ++p) {
		const auto& site = Alignment[p];
		if (!site.get<Flags::Aligned>() || site.base == genometools::Base::N || Alignment.positionInRef(p) < from()) continue;	

		if (Alignment.positionInRef(p) >= to()) {
			_overlap.push_back(Alignment);
			break;
		}
		const auto iWindow = Alignment.positionInRef(p) - from();

		if (_masked[iWindow]) continue;
		if (_entries[iWindow].depth() == _upToDepth) continue;

		++_numBases;
		if (_entries[iWindow].empty()) ++_sitesData;
		if (_entries[iWindow].depth() == 1) ++_sites2Plus;

		_entries[iWindow].add(site);
		if (_tagReads) _readIDs[iWindow].push_back(_numReads%maxReadID);
	}
	++_numReads;
}
void TBamWindow::filter() {
	using coretools::P;
	_iSite = -1;
	for (size_t i = 0; i < _entries.size(); ++i) {
		if (_masked[i]) continue;

		if (_downProb < P(1)) {
			_entries[i].downsample(_downProb);
		}
		if (_depthFilter.outside(_entries[i].depth())) {
			mask(i);
			continue;
		}
		if (_skipEmpty && _entries[i].empty()) {
			// do not mask, they will be counted as missing data!
			continue;
		}

		if (_iSite == size_t(-1)) {
			_iSite = i;
		}
	}
}

void TBamWindow::nextSite() {
	++_iSite;
	for (; _iSite < _entries.size(); ++_iSite) {
		if (_masked[_iSite]) continue;
		if (_skipEmpty && _entries[_iSite].empty()) continue;

		break;
	}
}

const GenotypeLikelihoods::TSite &TBamWindow::operator[](size_t I) const noexcept(coretools::noDebug) {
	DEBUG_ASSERT(I < size());
	return _entries[I];
}

GenotypeLikelihoods::TSite &TBamWindow::operator[](size_t I) noexcept(coretools::noDebug) {
	DEBUG_ASSERT(I < size());
	return _entries[I];
}

} // namespace BAM
