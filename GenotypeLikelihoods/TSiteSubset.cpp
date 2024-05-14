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
TSitePolymorphic::TSitePolymorphic(uint32_t refID, uint32_t position, char ref, char alt,
                                   const genometools::TChromosomes &Chromosomes)
    : TGenomePosition(refID, position) {
	// parse first allele (reference)
	_ref = genometools::char2base(ref);
	if (_ref == genometools::Base::N) {
		UERROR("Unknown allele1 '", ref, "' on ", asFormattedString(Chromosomes), "!");

		// parse second allele
		_alt = genometools::char2base(alt);
		if (_alt == genometools::Base::N) {
			UERROR("Unknown allele2 '", alt, "' on ", asFormattedString(Chromosomes), "!");
		}
		if (_ref == _alt) { UERROR("Site on ", asFormattedString(Chromosomes), " is invariant!"); }
	}
}

TSiteMonomorphic::TSiteMonomorphic(uint32_t refID, uint32_t position, char ref,
                                   const genometools::TChromosomes &Chromosomes)
    : TGenomePosition(refID, position) {
	// parse first allele (reference)
	_ref = genometools::char2base(ref);
	if (_ref == genometools::Base::N) {
		UERROR("Unknown allele1 '", ref, "' on ", asFormattedString(Chromosomes), "!");
	}
}

} // namespace SiteSubset
