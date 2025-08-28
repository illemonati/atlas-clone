#include "TBamWindow.h"
#include "coretools/Main/TParameters.h"

namespace BAM {

using coretools::instances::parameters;

void TBamWindow::_makeBedOrAlleles(const genometools::TChromosomes &Chromosomes) {
	if (parameters().exists("regions")) { _regions.parse(parameters().get("regions"), Chromosomes); }
	if (parameters().exists("mask")) {
		if (_regions.empty()) {
			_regions.parse(parameters().get("mask"), Chromosomes);
			_regions.flip(Chromosomes);
		} else {
			genometools::TBed mask(parameters().get("mask"), Chromosomes);
			_regions.addMask(mask);
		}
	}
	if (parameters().exists("alleles")) {
		coretools::user_assert(_regions.empty(), "Cannot define regions and alleles at the same time!");
		_alleles.parse(parameters().get("alleles"), Chromosomes);
	}
}

TBamWindow::TBamWindow(const genometools::TChromosomes &Chromosomes) {
	_makeBedOrAlleles(Chromosomes);
}

void TBamWindow::move(genometools::TGenomeWindow Window) {
	_from = Window.from();
	const auto length   = Window.size();

	_entries.assign(length, {});
	// clear everything
	_masked.clear();
	_depthTot   = 0;
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
		_masked.assign(Window.size(), true);

		for (auto it = _alleles.begin(Window); it != _alleles.end() && Window.within(it->position); ++it) {
			_masked[it->position - Window.from()] = false;
			_numMasked--;
		}
	}

	const auto copy = _overlap;
	_overlap.clear();
	for (const auto& aln: copy) {
		add(aln);
	}
	//_overlap.clear();
}

void TBamWindow::add(const TAlignment& Alignment) {
	for (size_t p = 0; p < Alignment.parsedLength(); ++p) {
		const auto& site = Alignment[p];
		if (!site.get<Flags::Aligned>() || site.base == genometools::Base::N || Alignment.positionInRef(p) < from()) continue;	

		if (Alignment.positionInRef(p) >= to()) {
			_overlap.push_back(Alignment);
			break;
		}

		const auto iWindow = Alignment.positionInRef(p) - from();
		if (!_masked.empty() && _masked[p]) continue;

		++_depthTot;
		if (_entries[iWindow].empty()) ++_sitesData;
		if (_entries[iWindow].depth() == 1) ++_sites2Plus;
		_entries[iWindow].add(Alignment[p]);
	}
}

}
