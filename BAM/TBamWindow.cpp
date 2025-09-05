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
	: _upToDepth(parameters().get<size_t>("readUpToDepth", 1000)) {
	_makeBedOrAlleles(Chromosomes);

	logfile().list("Will read data up to depth ", _upToDepth,
	               " and ignore additional bases. (parameter 'readUpToDepth')");
}

void TBamWindow::move(genometools::TGenomeWindow Window, const genometools::TFastaReader& Reference, bool FilterCpG) {
	using genometools::Base;
	_from = Window.from();

	_entries.assign(Window.size(), {});
	// clear everything
	_numBases   = 0;
	_sitesData  = 0;
	_sites2Plus = 0;
	_numMasked  = 0;

	if (Reference) {
		const auto view = Reference.view(window());
		for (size_t i = 0; i < size(); ++i) {
			_entries[i].refBase = view[i];
		}
	}

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

	if (FilterCpG) {
		for(size_t i = 1; i < size(); ++i) {
			if (Reference[from() + i - 1] == Base::C && Reference[from() + i] == Base::G) {
				_masked[i - 1] = true; // C
				_masked[i]     = true; // G
			}
		}

		// borders
		if (from().position() > 0 && Reference[from() - 1] == Base::C && Reference[from()] == Base::G)
			_masked.front() = true;
		if (Reference[to() - 1] == Base::C && Reference[to()] == Base::G) _masked.back() = true;
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
		DEBUG_ASSERT(_entries.size() > iWindow);
		DEBUG_ASSERT(_masked.size() > iWindow);
		if (_masked[iWindow]) continue;
		if (_entries[iWindow].depth() == _upToDepth) continue;

		++_numBases;
		if (_entries[iWindow].empty()) ++_sitesData;
		if (_entries[iWindow].depth() == 1) ++_sites2Plus;

		_entries[iWindow].add(site);
	}
}
}
