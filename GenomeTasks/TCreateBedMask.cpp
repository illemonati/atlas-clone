/*
 * TCreateBedMask.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#include "TCreateBedMask.h"
#include <vector>
#include "genometools/GenomePositions/TGenomePosition.h"
#include "TGenotypeData.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TSite.h"
#include "TWindow.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenomeTasks{

using coretools::str::toString;
using coretools::instances::parameters;
using coretools::instances::logfile;

//--------------------------------------
// TCreateBedMask
//--------------------------------------
TCreateBedMask::TCreateBedMask():TGenome_windows(){
	_minDepth = parameters().get<uint32_t>("minDepth", 2);
	_bed.addChromosomes(_chromosomes);
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

	if(_maxDepth < _minDepth){
		UERROR("maxDepthForMask must be > minDepthForMask!");
	}

	if(parameters().exists("maxDepth") || parameters().exists("minDepth"))
		UERROR("Cannot mask sites for sequencing depth (parameters 'minDepth' and 'maxDepth') while creating the mask!");
};

void TCreateDepthBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto& s : _window){
		if(s.depth() < _minDepth || s.depth() > _maxDepth){
			_bed.add(_window.from() + p);
		}
		++p;
	}
};

void TCreateDepthBedMask::createDepthMask(){
	_createMask("maxDepth" + toString(_maxDepth) + "_depthMask.bed");
};

//--------------------------------------
// TCreateInvariantBedMask
//--------------------------------------
TCreateInvariantBedMask::TCreateInvariantBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepth) + " (parameter 'minDepthForMask') for which a single allele was observed (invariant).");

	if(_minDepth < 2){
		UERROR("minDepthForMask must be >= 2 to assess variant / invariant status!");
	}
};

void TCreateInvariantBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto& s : _window){
		if(s.depth() >= _minDepth){
			const auto bCounts = s.countAlleles();
			if(coretools::numNonZero(bCounts) == 1){
				_bed.add(_window.from() + p);
			}
		}
		++p;
	}
};

void TCreateInvariantBedMask::createInvariantMask(){
	_createMask("invariantMask.bed");
};

//--------------------------------------
// TCreateVariantBedMask
//--------------------------------------
TCreateVariantBedMask::TCreateVariantBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepth) + " (parameter 'minDepthForMask') for which multiple alleles were observed (variant).");

	if(_minDepth < 2){
		UERROR("minDepthForMask must be >= 2 to assess variant / invariant status!");
	}
};

void TCreateVariantBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto& s : _window){
		if(s.depth() >= _minDepth){
			const auto bCounts = s.countAlleles();
			if(coretools::numNonZero(bCounts) > 1){
				_bed.add(_window.from() + p);
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

	if(_minDepth < 1){
		UERROR("maxDepthForMask must be > 1 to check for ref / non-ref status!");
	}

	_openReference(true);
};

void TCreateNonRefBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto& s : _window){
		if(s.depth() >= _minDepth){
			if(s.refDepth() < s.depth()){
				_bed.add(_window.from() + p);
			}
		}
		++p;
	}
};

void TCreateNonRefBedMask::createVariantMask(){
	_createMask("nonRefMask.bed");
};

}; // end namespace
