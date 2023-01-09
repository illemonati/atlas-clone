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
	_minDepthForMask = parameters().getParameter<uint32_t>("minDepthForMask");
};

void TCreateBedMask::_createMask(const std::string fileTag){
	_traverseBAMWindows();

	//write mask
	std::string filename = _outputName + "_minDepth"+ toString(_minDepthForMask) + "_" + fileTag + ".bed";
	logfile().listFlush("Writing mask to BED file '" + filename + "' ...");

	_bed.write(filename);
	logfile().done();
};

//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
TCreateDepthBedMask::TCreateDepthBedMask():TCreateBedMask(){
	_maxDepthForMask = parameters().getParameter<uint32_t>("maxDepthForMask");
	logfile().list("Will create a mask for all sites with depth outside the range [" + toString(_minDepthForMask) + ", " + toString(_maxDepthForMask) + "].");

	if(_maxDepthForMask < _minDepthForMask){
		UERROR("maxDepthForMask must be > minDepthForMask!");
	}

	if(parameters().parameterExists("maxDepth") || parameters().parameterExists("minDepth"))
		UERROR("Cannot mask sites for sequencing depth (parameters 'minDepth' and 'maxDepth') while creating the mask!");
};

void TCreateDepthBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto& s : _window){
		if(s.depth() < _minDepthForMask || s.depth() > _maxDepthForMask){
			_bed.add(_window.from() + p);
		}
		++p;
	}
};

void TCreateDepthBedMask::createDepthMask(){
	_createMask("maxDepth" + toString(_maxDepthForMask) + "_depthMask.bed");
};

//--------------------------------------
// TCreateInvariantBedMask
//--------------------------------------
TCreateInvariantBedMask::TCreateInvariantBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepthForMask) + " (parameter 'minDepthForMask') for which a single allele was observed (invariant).");

	if(_minDepthForMask < 2){
		UERROR("minDepthForMask must be >= 2 to assess variant / invariant status!");
	}
};

void TCreateInvariantBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto& s : _window){
		if(s.depth() >= _minDepthForMask){
			s.countAlleles(_baseCounts);
			if(coretools::numNonZero(_baseCounts) == 1){
				_bed.add(_window.from() + p);
			}
		}
		++p;
	}
};

void TCreateInvariantBedMask::createInvariantMask(){
	_createMask("_invariantMask.bed");
};

//--------------------------------------
// TCreateVariantBedMask
//--------------------------------------
TCreateVariantBedMask::TCreateVariantBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepthForMask) + " (parameter 'minDepthForMask') for which multiple alleles were observed (variant).");

	if(_minDepthForMask < 2){
		UERROR("minDepthForMask must be >= 2 to assess variant / invariant status!");
	}
};

void TCreateVariantBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto& s : _window){
		if(s.depth() >= _minDepthForMask){
			s.countAlleles(_baseCounts);
			if(coretools::numNonZero(_baseCounts) > 1){
				_bed.add(_window.from() + p);
			}
		}
		++p;
	}
};

void TCreateVariantBedMask::createVariantMask(){
	_createMask("_variantMask.bed");
};

//--------------------------------------
// TCreateNonRefBedMask
//--------------------------------------
TCreateNonRefBedMask::TCreateNonRefBedMask():TCreateBedMask(){
	logfile().list("Will create a mask of all sites with depth >= " + toString(_minDepthForMask) + " (parameter 'minDepthForMask') for which at least one non-ref allele was observed.");

	if(_minDepthForMask < 1){
		UERROR("maxDepthForMask must be > 1 to check for ref / non-ref status!");
	}

	_openReference(true);
};

void TCreateNonRefBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto& s : _window){
		if(s.depth() >= _minDepthForMask){
			if(s.refDepth() < s.depth()){
				_bed.add(_window.from() + p);
			}
		}
		++p;
	}
};

void TCreateNonRefBedMask::createVariantMask(){
	_createMask("_nonRefMask.bed");
};

}; // end namespace
