/*
 * TPMDSCalculator.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */


#include "TPMDSCalculator.h"

namespace GenomeTasks{

//----------------------------------------------
// TPMDSCalculator
//----------------------------------------------
//TODO: should that filter pairs as in TBamFilter?
TPMDSCalculator::TPMDSCalculator(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_parsed(Parameters, Logfile, RandomGenerator){
	//get parameters
	_pi = Parameters.getParameterDoubleWithDefault("pi", 0.001);
	_logfile->list("Running PMDS with rate of polymorphism (pi) = " + toString(_pi));
	_minPMDS = Parameters.getParameterDoubleWithDefault("minPMDS", -10000);
	_maxPMDS = Parameters.getParameterDoubleWithDefault("maxPMDS", 10000);
	_logfile->list("Filtering out reads with PMDS outside the range [" + toString(_minPMDS) + ", " + toString(_maxPMDS) + "].");
};

double TPMDSCalculator::_calculatePMDS(){
	//calculate PMDS (is in log)
	double PMDS = 0.0;
	uint16_t d = 0;
	for(auto& b : _alignment){
		//limit to aligned positions
		Base ref = _genoMap.getBase(_alignment.referenceAtInternalPos(d));
		if(b.isAligned() && b.base != N && ref != N){
			PMDS +=  _genotypeLikelihoodCalculator.calculateLogPMDS(b, ref, _pi);
		}
		++d;
	}

	return PMDS;
};

void TPMDSCalculator::_handleAlignment(){
	//calc PMD
	double PMDS = _calculatePMDS();

	//filter
	if(PMDS < _minPMDS || PMDS < _maxPMDS){
		_bamFile.curFilterOut();
	} else {
		//update and write
		//TODO: discuss if DS is the right tag. User-defined tags should have lower case letters, but we need to maintain consistency with other tools
		_bamFile.curAddSamField("DS", PMDS);
		_bamFile.writeCurAlignment(_outBam);
	}
};

void TPMDSCalculator::calculateLogPMDS(){
	//parse bam file and calculate PMDS for each read (seeSkoglund et al. 2014)
	//write new bam file with PMDS score added
	//parser.add_option("--writesamfield", action="store_true", dest="writesamfield",help="add 'DS:Z:<PMDS>' field to SAM output, will overwrite if already present",default=False)

	//open a bam file for writing
	_openBamForWriting(_outputName + "_PMDS.bam", _outBam);
	_bamFile.setExternalFilterReason("PMDS outside range [" + toString(_minPMDS) + ", " + toString(_maxPMDS) + "]");

	//traverse BAM
	_traverseBAMPassedQC();

	//report
	_outBam.close(_logfile);
};

}; // end namespace
