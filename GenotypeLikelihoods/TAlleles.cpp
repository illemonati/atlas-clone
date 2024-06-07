#include "TAlleles.h"

#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/stringManipulations.h"
#include "genometools/GenomePositions/TChromosomes.h"

namespace SiteSubset {
using genometools::TChromosomes;
using genometools::Base;
using namespace coretools::str;
using coretools::instances::logfile;

void TAlleles::parse(std::string_view Filename, const TChromosomes &Chromosomes, Morphic M) {
	_chromosomes.assign(Chromosomes.size(), "");
	for (const auto &chr : Chromosomes) { _chromosomes[chr.refID()] = chr.name(); }
	parse(Filename, M);
}

void TAlleles::parse(std::string_view Filename, coretools::TConstView<std::string> Chromosomes, Morphic M) {
	_chromosomes.clear();
	_chromosomes.reserve(Chromosomes.size());
	for (const auto &chr : Chromosomes) { _chromosomes.push_back(chr); }
	parse(Filename, M);
}

void TAlleles::parse(std::string_view Filename, Morphic M) {
	_allSites.clear();
	_sites.clear();
	if (Filename.empty()) { return; }

	coretools::TInputFile in(Filename, coretools::FileType::Header);

    std::array<size_t, 4> idx;
	idx[0] = in.indexOfFirstMatch({"CHR", "Chr", "chr", "CHROMOSOME", "Chromosome", "chromosome", "#CHROM"});
	idx[1] = in.indexOfFirstMatch({"POS", "Pos", "pos", "POSITION", "Position", "position"});
	idx[2] = in.indexOfFirstMatch({"REF", "Ref", "ref", "REFERENCE", "Reference", "reference", "Allele", "Allele1"});

	if (M == Morphic::Poly) {
		idx[3] = in.indexOfFirstMatch({"ALT", "Alt", "alt", "ALTERNATIVE", "Alternative", "alternative", "Allele2",
									   "DERIVED", "Derived", "derived"});
	}

    std::vector<std::string> unknownChrom;
    for (; !in.empty(); in.popFront()) {
		const auto chr = in.get(idx[0]);
		const auto rID = refID(chr);
		if (rID == _chromosomes.size()) {
			unknownChrom.emplace_back(chr);
			continue;
		}
		
		const auto pos = in.get<size_t>(idx[1]) + 1;
		const auto ref = genometools::char2base(in.get(idx[2]).front());
		if (M == Morphic::Mono) {
			_allSites.emplace_back(rID, pos, ref, Base::N);
		} else {
			const auto alt = genometools::char2base(in.get(idx[3]).front());
			if (alt == ref) {
				UERROR("reference and alternative base are the same in sites file '", Filename, "': ", in.frontRaw(), "!");
			}
			_allSites.emplace_back(rID, pos, ref, alt);
		}
	}

    std::sort(_allSites.begin(), _allSites.end());

	_sites.assign(_chromosomes.size(), {end(), 0});
	_nChrWindows = 0;
	auto refID   = _allSites.front().position.refID();
	auto data    = _allSites.data();
	size_t n     = 1;
	for (size_t i = 1; i < _allSites.size(); ++i) {
		if (_allSites[i].position == _allSites[i - 1].position) {
			UERROR("Duplicates in sites file '", Filename, "': position ", name(_allSites[i].position.refID()), ":",
				   _allSites[i].position.position(), " is present multiple times!");
		}
		if (_allSites[i].position.refID() != refID) {
			_sites[refID] = {data, n};
			data          = _allSites.data() + i;
			n             = 1;
			refID         = _allSites[i].position.refID();
			++_nChrWindows;
		} else {
			n++;
		}
	}
	_sites[refID] = {data, n};
	++_nChrWindows;


	if (unknownChrom.size() > 0) {
		logfile().warning("Sites of the following chromosomes were ignored because these chromosomes were not in the "
						  "BAM header / GLF index: ",
						  unknownChrom, ".");
	}
}
}
