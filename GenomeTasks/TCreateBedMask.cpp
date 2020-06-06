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
TCreateBedMask::TCreateBedMask(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Params, Logfile, RandomGenerator){
	_minDepthForMask = Params.getParameterInt("minDepthForMask");
};


//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
TCreateDepthBedMask::TCreateDepthBedMask(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TCreateBedMask(Params, Logfile, RandomGenerator){
	_maxDepthForMask = Params.getParameterInt("maxDepthForMask");
	_logfile->list("Will create a mask for all sites with outside the range [" + toString(_minDepthForMask) + ", " + toString(_maxDepthForMask) + "].");

	if(_maxDepthForMask < _minDepthForMask){
		throw "maxDepthForMask must be > minDepthForMask!";
	}

	if(Params.parameterExists("maxDepth") || Params.parameterExists("minDepth"))
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
	_traverseBAMWindows();

	//write mask
	std::string filename = _outputName + "_minDepth"+ toString(_minDepthForMask) + "_maxDepth" + toString(_maxDepthForMask) + "_depthMask.bed";
	_logfile->listFlush("Writing mask to BED file '" + filename + "' ...");

	_bed.write(filename);
	_logfile->done();
};

//--------------------------------------
// TCreateConservedBedMask
//--------------------------------------
TCreateConservedBedMask::TCreateConservedBedMask(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator):TCreateBedMask(Params, Logfile, RandomGenerator){
	_logfile->list("Will create a mask of all sites with depth >= " + toString(_minDepthForMask) + " (parameter 'minDepthForMask') for which a single allele is observed (invariant).");

	if(_minDepthForMask < 2){
		throw "maxDepthForMask must be >= 2!";
	}
};

void TCreateConservedBedMask::_handleWindow(){
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

void TCreateConservedBedMask::createConservedMask(){
	_traverseBAMWindows();

	//write mask
	std::string filename = _outputName + "_minDepth"+ toString(_minDepthForMask) + "_invariantMask.bed";
	_logfile->listFlush("Writing mask to BED file '" + filename + "' ...");

	_bed.write(filename);
	_logfile->done();
};

}; // end namespace
