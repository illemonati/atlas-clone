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
	TCreateBedMask::TCreateBedMask() :  _bed(_genome.bamFile().chromosomes()){
	_minDepth = parameters().get<uint32_t>("minDepth", 2);
};

void TCreateBedMask::_createMask(const std::string fileTag){
	_traverseBAMWindows();

	//write mask
	std::string filename = _genome.outputName() + "_minDepth"+ toString(_minDepth) + "_" + fileTag + ".bed";
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

	user_assert(!parameters().exists("maxDepth") && !parameters().exists("minDepth"), "Cannot mask sites for sequencing depth (parameters 'minDepth' and 'maxDepth') while creating the mask!");
}

void TCreateDepthBedMask::_handleWindow(GenotypeLikelihoods::TWindow& window){
	uint32_t p = 0;
	for(auto& s : window){
		if(s.depth() < _minDepth || s.depth() > _maxDepth){
			_bed.add(window.from() + p);
		}
		++p;
	}
}

void TCreateDepthBedMask::createDepthMask(){
	_createMask("maxDepth" + toString(_maxDepth) + "_depthMask.bed");
}

//--------------------------------------
// TCreateInvariantBedMask
//--------------------------------------
TCreateInvariantBedMask::TCreateInvariantBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepth) + " (parameter 'minDepthForMask') for which a single allele was observed (invariant).");

	user_assert(_minDepth >= 2, "minDepthForMask must be >= 2 to assess variant / invariant status!");
}

void TCreateInvariantBedMask::_handleWindow(GenotypeLikelihoods::TWindow& window){
	uint32_t p = 0;
	for(auto& s : window){
		if(s.depth() >= _minDepth){
			const auto bCounts = s.countAlleles();
			if(coretools::numNonZero(bCounts) == 1){
				_bed.add(window.from() + p);
			}
		}
		++p;
	}
}

void TCreateInvariantBedMask::createInvariantMask(){
	_createMask("invariantMask.bed");
}

//--------------------------------------
// TCreateVariantBedMask
//--------------------------------------
TCreateVariantBedMask::TCreateVariantBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepth) + " (parameter 'minDepthForMask') for which multiple alleles were observed (variant).");

	user_assert (_minDepth >= 2, "minDepthForMask must be >= 2 to assess variant / invariant status!");
}

void TCreateVariantBedMask::_handleWindow(GenotypeLikelihoods::TWindow& window){
	uint32_t p = 0;
	for(auto& s : window){
		if(s.depth() >= _minDepth){
			const auto bCounts = s.countAlleles();
			if(coretools::numNonZero(bCounts) > 1){
				_bed.add(window.from() + p);
			}
		}
		++p;
	}
};

void TCreateVariantBedMask::createVariantMask(){
	_createMask("variantMask.bed");
};

//--------------------------------------
// TCreateNonRefBedMask
//--------------------------------------
TCreateNonRefBedMask::TCreateNonRefBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepth) + " (parameter 'minDepthForMask') for which at least one non-ref allele was observed.");

	user_assert(_minDepth > 1, "maxDepthForMask must be > 1 to check for ref / non-ref status!");
	_windows.requireReference();
}

void TCreateNonRefBedMask::_handleWindow(GenotypeLikelihoods::TWindow& window){
	uint32_t p = 0;
	for(auto& s : window){
		if(s.depth() >= _minDepth){
			if(s.refDepth() < s.depth()){
				_bed.add(window.from() + p);
			}
		}
		++p;
	}
};

void TCreateNonRefBedMask::createVariantMask(){
	_createMask("nonRefMask.bed");
};

}; // end namespace
