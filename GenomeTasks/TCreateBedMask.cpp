/*
 * TCreateBedMask.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#include "TCreateBedMask.h"

namespace GenomeTasks{

//--------------------------------------
// TCreateBedMask
//--------------------------------------
TCreateBedMask::TCreateBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Parameters, Logfile, RandomGenerator){
	_minDepthForMask = Parameters.getParameterInt("minDepthForMask");
};

void TCreateBedMask::_createMask(const std::string fileTag){
	_traverseBAMWindows();

	//write mask
	std::string filename = _outputName + "_minDepth"+ toString(_minDepthForMask) + "_" + fileTag + ".bed";
	_logfile->listFlush("Writing mask to BED file '" + filename + "' ...");

	_bed.write(filename);
	_logfile->done();
};


//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
TCreateDepthBedMask::TCreateDepthBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCreateBedMask(Parameters, Logfile, RandomGenerator){
	_maxDepthForMask = Parameters.getParameterInt("maxDepthForMask");
	_logfile->list("Will create a mask for all sites with depth outside the range [" + toString(_minDepthForMask) + ", " + toString(_maxDepthForMask) + "].");

	if(_maxDepthForMask < _minDepthForMask){
		throw "maxDepthForMask must be > minDepthForMask!";
	}

	if(Parameters.parameterExists("maxDepth") || Parameters.parameterExists("minDepth"))
		throw "Cannot mask sites for sequencing depth (parameters 'minDepth' and 'maxDepth') while creating the mask!";
};

void TCreateDepthBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto s = _window.begin(); s!=_window.end(); ++s){
		if(s->depth() < _minDepthForMask || s->depth() > _maxDepthForMask){
			_bed.addSite(_window.posInRef(p));
		}
	}
};

void TCreateDepthBedMask::createDepthMask(){
	_createMask("maxDepth" + toString(_maxDepthForMask) + "_depthMask.bed");
};

//--------------------------------------
// TCreateInvariantBedMask
//--------------------------------------
TCreateInvariantBedMask::TCreateInvariantBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCreateBedMask(Parameters, Logfile, RandomGenerator){
	_logfile->list("Will create a mask of all sites with depth >= " + toString(_minDepthForMask) + " (parameter 'minDepthForMask') for which a single allele was observed (invariant).");

	if(_minDepthForMask < 2){
		throw "minDepthForMask must be >= 2 to assess variant / invariant status!";
	}
};

void TCreateInvariantBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto s = _window.begin(); s!=_window.end(); ++s){
		if(s->depth() >= _minDepthForMask){
			s->countAlleles(_baseCounts);
			if(_baseCounts.numAlleles() == 1){
				_bed.addSite(_window.posInRef(p));
			}
		}
	}
};

void TCreateInvariantBedMask::createInvariantMask(){
	_createMask("_invariantMask.bed");
};

//--------------------------------------
// TCreateVariantBedMask
//--------------------------------------
TCreateVariantBedMask::TCreateVariantBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCreateBedMask(Parameters, Logfile, RandomGenerator){
	_logfile->list("Will create a mask of all sites with depth >= " + toString(_minDepthForMask) + " (parameter 'minDepthForMask') for which multiple alleles were observed (variant).");

	if(_minDepthForMask < 2){
		throw "minDepthForMask must be >= 2 to assess variant / invariant status!";
	}
};

void TCreateVariantBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto s = _window.begin(); s!=_window.end(); ++s){
		if(s->depth() >= _minDepthForMask){
			s->countAlleles(_baseCounts);
			if(_baseCounts.numAlleles() > 1){
				_bed.addSite(_window.posInRef(p));
			}
		}
	}
};

void TCreateVariantBedMask::createVariantMask(){
	_createMask("_variantMask.bed");
};

//--------------------------------------
// TCreateNonRefBedMask
//--------------------------------------
TCreateNonRefBedMask::TCreateNonRefBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TCreateBedMask(Parameters, Logfile, RandomGenerator){
	_logfile->list("Will create a mask of all sites with depth >= " + toString(_minDepthForMask) + " (parameter 'minDepthForMask') for which at least one non-ref allele was observed.");

	if(_minDepthForMask < 1){
		throw "maxDepthForMask must be > 1 to check for ref / non-ref status!";
	}

	_openReference(true);
};

void TCreateNonRefBedMask::_handleWindow(){
	uint32_t p = 0;
	for(auto s = _window.begin(); s!=_window.end(); ++s){
		if(s->depth() >= _minDepthForMask){
			if(s->refDepth() < s->depth()){
				_bed.addSite(_window.posInRef(p));
			}
		}
	}
};

void TCreateNonRefBedMask::createVariantMask(){
	_createMask("_nonRefMask.bed");
};

}; // end namespace
