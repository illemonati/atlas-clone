/*
 * TPMDSCalculator.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */


#include "TPMDSCalculator.h"

#include <stdint.h>
#include <exception>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "TAlignment.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TLog.h"
#include "TSequencedBase.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks{

using coretools::str::toString;
using coretools::instances::logfile;
using coretools::instances::parameters;

//----------------------------------------------
// TPMDSCalculator
//----------------------------------------------
//TODO: should that filter pairs as in TBamFilter?
	TPMDSCalculator::TPMDSCalculator():TGenome_parsed(), _outBam(_genome.outputName() + "_PMDS.bam", _genome.bamFile()) {
	//get parameters
	_pi = parameters().get<coretools::Probability>("pi", coretools::Probability(0.001));
	logfile().list("Running PMDS with rate of polymorphism (pi) = " + toString(_pi));
	if(parameters().exists("filterPMDS")){
		_filterRange.set(parameters().get("filterPMDS"));
		_doFilter = true;
		logfile().list("Filtering out reads with PMDS outside the range " + _filterRange.rangeString() + ".");
	} else {
		_doFilter = true;
		logfile().list("Not applying any filter on PMDs when writing BAM file. (use 'filterPMDS' to filter)");
	}
	_parser.openReference(true);
};

double TPMDSCalculator::_calculatePMDS(BAM::TAlignment& alignment){
	//calculate PMDS (is in log)
	double PMDS = 0.0;
	for (size_t d = 0; d < alignment.size(); ++d) {
		if (alignment.isAlignedAtInternalPos(d)) {
			PMDS += _parser.errorModels().calculateLogPMDS(alignment[d], alignment.referenceAtInternalPos(d),
																   _pi);
		}
	}
	return PMDS;
};

void TPMDSCalculator::_handleAlignment(BAM::TAlignment& alignment){
	//calc PMD
	const auto PMDS = _calculatePMDS(alignment);

	//filter
	if(_doFilter && !_filterRange.within(PMDS)){
		_genome.bamFile().curFilterOut();
	} else {
		//update and write
		//TODO: discuss if DS is the right tag. User-defined tags should have lower case letters, but we need to maintain consistency with other tools
		_genome.bamFile().curAddSamField("DS", PMDS);
		_genome.bamFile().writeCurAlignment(_outBam);
	}
};

void TPMDSCalculator::run(){
	//parse bam file and calculate PMDS for each read (seeSkoglund et al. 2014)
	//write new bam file with PMDS score added
	//parser.add_option("--writesamfield", action="store_true", dest="writesamfield",help="add 'DS:Z:<PMDS>' field to SAM output, will overwrite if already present",default=False)

	//open a bam file for writing
	_genome.bamFile().setExternalFilterReason("PMDS outside range " + _filterRange.rangeString());

	//traverse BAM
	_traverseBAMPassedQC();
};

}; // end namespace
