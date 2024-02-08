/*
 * TSiteSubset.cpp
 *
 *  Created on: Nov 25, 2019
 *      Author: vivian
 */

#include "TSiteSubset.h"

namespace GenotypeLikelihoods::SiteSubset {

//-----------------------------------------------
// TSitePolymorphic / TSiteMonomorphic
//-----------------------------------------------
TSitePolymorphic::TSitePolymorphic(uint32_t refID, uint32_t position, const std::vector<std::string_view> &Line,
                                   const genometools::TChromosomes &Chromosomes)
    : TGenomePosition(refID, position) {
	// parse first allele (reference)
	_ref = genometools::char2base(Line[2][0]);
	if (_ref == genometools::Base::N) {
		UERROR("Unknown allele1 '", Line[2], "' on ", asFormattedString(Chromosomes), "!");

		// parse second allele
		_alt = genometools::char2base(Line[3][0]);
		if (_alt == genometools::Base::N) {
			UERROR("Unknown allele2 '", Line[3], "' on ", asFormattedString(Chromosomes), "!");
		}
		if (_ref == _alt) { UERROR("Site on ", asFormattedString(Chromosomes), " is invariant!"); }
	}
}

void TSitePolymorphic::write(coretools::TOutputFile &out, const genometools::TChromosomes &Chromosomes) const {
	out.writeln(asFormattedString(Chromosomes, "\t"), _ref, _alt);
}

TSiteMonomorphic::TSiteMonomorphic(uint32_t refID, uint32_t position, const std::vector<std::string_view> &Line,
                                   const genometools::TChromosomes &Chromosomes)
    : TGenomePosition(refID, position) {
	// parse first allele (reference)
	_ref = genometools::char2base(Line[2][0]);
	if (_ref == genometools::Base::N) {
		UERROR("Unknown allele1 '", Line[2], "' on ", asFormattedString(Chromosomes), "!");
	}
}

void TSiteMonomorphic::write(coretools::TOutputFile &out, const genometools::TChromosomes &Chromosomes) const {
	out.writeln(asFormattedString(Chromosomes, "\t"), _ref);
}

} // namespace SiteSubset
