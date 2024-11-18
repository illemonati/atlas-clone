#include "TBamReadMask.h"

#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;

void TBamReadMask::setMasks(const genometools::TChromosomes &Chromosomes){
	// normal mask
	if (parameters().exists("mask") || parameters().exists("regions")) {
		std::string filename;

		if (parameters().exists("mask")) {
			logfile().startIndent("Will apply BED mask:");
			if (parameters().exists("regions")) UERROR("Cannot use mask and regions at the same time.");

			filename = parameters().get<std::string>("mask");
			logfile().list("Will ignore reads overlapping with windows listed in BED file '" + filename + "'. (parameter 'mask')");
			_doMasking       = true;
			_considerRegions = false;			
		} else {
			filename = parameters().get<std::string>("regions");
			logfile().startIndent("Will limit to BED regions:");
			logfile().list("Will limit analysis to reads overlapping with windows listed in BED file '" + filename +
								  "' (parameter 'regions'):");
			_doMasking       = false;
			_considerRegions = true;
		}

		// read file
		logfile().listFlush("Reading file ...");
		_mask.parse(filename, Chromosomes);
		logfile().done();
		logfile().conclude("Read ", _mask.size(), " sites on ", _mask.NChrWindows(), " chromosomes.");		

		// read porosity		
		_porousProb = parameters().get("maskPorosity", coretools::Probability(0.0));
		if(_porousProb > 0.0){
			_maskIsPorous = true;
			if(_doMasking){			
				logfile().list("Mask will be porous: overlapping reads will be kept with probability ", _porousProb, ". (parameter 'maskPorosity')");
			} else {
				logfile().list("Regions will be porous: non-overlapping reads will be kept with probability ", _porousProb, ". (parameter 'maskPorosity')");
			}

		} else {
			_maskIsPorous = false;
			if(_doMasking){
				logfile().list("Mask will be strict: no overlapping reads will be kept. (use 'maskPorosity' to make it porous)");
			} else {
				logfile().list("Regions will be strict: only overlapping reads will be kept. (use 'maskPorosity' to make it porous)");
			}
		}
		logfile().endIndent();
	} else {
		logfile().list("Will use reads of the entire genome. (limit with 'mask' or 'regions')");
		_doMasking       = false;
		_considerRegions = false;
		_maskIsPorous = false;
		_porousProb = coretools::Probability(0.0);
	}
}

bool TBamReadMask::keepSingle(const genometools::TGenomeWindow& aln) const { 
	bool discard = (_doMasking && _mask.overlaps(aln)) || (_considerRegions && !_mask.overlaps(aln));
	if (discard) {
		return _maskIsPorous && coretools::instances::randomGenerator().getRand() < _porousProb;
	} 
	return true;
}

bool TBamReadMask::keepPaired(const genometools::TGenomeWindow& aln, const genometools::TGenomeWindow& mate) const {
	// if masking: neither of both can overlap
	// if regions: at least one overlaps
	bool discard = (_doMasking && (_mask.overlaps(aln) || _mask.overlaps(mate))) || (_considerRegions && !_mask.overlaps(aln) && !_mask.overlaps(mate));
	if (discard) {
		return _maskIsPorous && coretools::instances::randomGenerator().getRand() < _porousProb;
	} 
	return true;
}

}
