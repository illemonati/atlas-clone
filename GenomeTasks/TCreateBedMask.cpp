/*
 * TCreateBedMask.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#include "TCreateBedMask.h"
#include "coretools/algorithms.h"

namespace GenomeTasks{

using coretools::str::toString;
using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::user_assert;

//--------------------------------------
// TCreateBedMask
//--------------------------------------
TCreateBedMask::TCreateBedMask() :  _bed(_siteTraverser.bamFile().chromosomes()){
	_minDepth = parameters().get<uint32_t>("minDepth", 2);
	_siteTraverser.requireSingleBAM();
	if (_minDepth > 0) _siteTraverser.filterEmpty();
};

void TCreateBedMask::_createMask(std::string_view fileTag){
	for (;!_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		for (;!_siteTraverser.endOfCurChr(); _siteTraverser.nextSite()) {
			_handleSite(_siteTraverser.site(), _siteTraverser.position());
			
		}
	}

	//write mask
	const auto filename = toString(_siteTraverser.outputName(), "_minDepth", _minDepth, "_", fileTag, ".bed");
	logfile().listFlush("Writing mask to BED file '" + filename + "' ...");

	_bed.write(filename);
	logfile().done();
};

//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
TCreateDepthBedMask::TCreateDepthBedMask():TCreateBedMask(){
	_maxDepth = parameters().get<uint32_t>("maxDepth", 1000000);
	logfile().list("Will create a mask for all sites with depth outside the range [" + toString(_minDepth) + ", " + toString(_maxDepth) + "].");

	user_assert(_maxDepth > _minDepth, "maxDepthForMask must be > minDepthForMask!");
	_siteTraverser.filterEmpty(false);
}

void TCreateDepthBedMask::_handleSite(const GenotypeLikelihoods::TSite &Site, genometools::TGenomePosition Position) {
	if (Site.depth() < _minDepth || Site.depth() > _maxDepth) { _bed.add(Position); }
}

void TCreateDepthBedMask::createDepthMask(){
	_createMask("maxDepth" + toString(_maxDepth, "Mask"));
}

//--------------------------------------
// TCreateInvariantBedMask
//--------------------------------------
TCreateInvariantBedMask::TCreateInvariantBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepth) + " (parameter 'minDepthForMask') for which a single allele was observed (invariant).");

	user_assert(_minDepth >= 2, "minDepthForMask must be >= 2 to assess variant / invariant status!");
}

void TCreateInvariantBedMask::_handleSite(const GenotypeLikelihoods::TSite &Site,
										  genometools::TGenomePosition Position) {
	if (Site.depth() >= _minDepth) {
		const auto bCounts = Site.countAlleles();
		if (coretools::numNonZero(bCounts) == 1) { _bed.add(Position); }
	}
}

void TCreateInvariantBedMask::createInvariantMask(){
	_createMask("invariantMask");
}

//--------------------------------------
// TCreateVariantBedMask
//--------------------------------------
TCreateVariantBedMask::TCreateVariantBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepth) + " (parameter 'minDepthForMask') for which multiple alleles were observed (variant).");

	user_assert (_minDepth >= 2, "minDepthForMask must be >= 2 to assess variant / invariant status!");
}

void TCreateVariantBedMask::_handleSite(const GenotypeLikelihoods::TSite &Site, genometools::TGenomePosition Position) {
	if (Site.depth() >= _minDepth) {
		const auto bCounts = Site.countAlleles();
		if (coretools::numNonZero(bCounts) > 1) { _bed.add(Position); }
	}
}

void TCreateVariantBedMask::createVariantMask(){
	_createMask("variantMask");
};

//--------------------------------------
// TCreateNonRefBedMask
//--------------------------------------
TCreateNonRefBedMask::TCreateNonRefBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepth) + " (parameter 'minDepthForMask') for which at least one non-ref allele was observed.");

	user_assert(_minDepth > 1, "maxDepthForMask must be > 1 to check for ref / non-ref status!");
	_siteTraverser.requireReference();
}

void TCreateNonRefBedMask::_handleSite(const GenotypeLikelihoods::TSite &Site, genometools::TGenomePosition Position) {
	if (Site.depth() >= _minDepth) {
		if (Site.refDepth() < Site.depth()) { _bed.add(Position); }
	}
}

void TCreateNonRefBedMask::createVariantMask(){
	_createMask("nonRefMask");
}

} // namespace GenomeTasks
