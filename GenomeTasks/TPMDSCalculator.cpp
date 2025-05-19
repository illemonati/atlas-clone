/*
 * TPMDSCalculator.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#include "TPMDSCalculator.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TParameters.h"
#include <cmath>

namespace GenomeTasks{

using coretools::str::toString;
using coretools::instances::logfile;
using coretools::instances::parameters;

//----------------------------------------------
// TPMDSCalculator
//----------------------------------------------
//TODO: should that filter pairs as in TBamFilter?
TPMDSCalculator::TPMDSCalculator(): _outBam(_genome.outputName() + "_PMDS.bam", _genome.bamFile()) {
	//get parameters
	if (!_genome.errorModels().postMortemDamageModels().hasPMD()) {
		UERROR("Cannot estimate PMDS without pmd model!");
	}


	_pi = parameters().get<coretools::Probability>("pi", coretools::Probability(0.001));
	logfile().list("Running PMDS with rate of polymorphism (pi) = " + toString(_pi));
	if(parameters().exists("filterPMDS")){
		_filterRange.set(parameters().get("filterPMDS"));
		_doFilter = true;
		logfile().list("Filtering out reads with PMDS outside the range " + _filterRange.rangeString() + ".");
	} else {
		_doFilter = false;
		logfile().list("Not applying any filter on PMDs when writing BAM file. (use 'filterPMDS' to filter)");
	}
	_parser.openReference(true);
};

double TPMDSCalculator::_calculatePMDS(BAM::TAlignment& alignment){
	//calculate PMDS (is in log)
	double PMDS = 0.0;
	for (size_t d = 0; d < alignment.size(); ++d) {
		if (alignment.isAlignedAtInternalPos(d)) {
			const auto ref = alignment.referenceAtInternalPos(d);
			if (ref == genometools::Base::N) continue;

			PMDS += _genome.errorModels().calculateLogPMDS(alignment[d], alignment.referenceAtInternalPos(d), _pi);
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
		if (PMDS > 0) {
			const size_t idx = std::lround(PMDS*10);
			if (_dPMDs_pos.size() <= idx) _dPMDs_pos.resize(idx + 1, 0);
			++_dPMDs_pos[idx];
		} else {
			const size_t idx = std::lround(-PMDS*10);
			if (_dPMDs_neg.size() <= idx) _dPMDs_neg.resize(idx + 1, 0);
			++_dPMDs_neg[idx];
		}
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

	coretools::TOutputFile hist(_genome.outputName() + "_PMDS_hist.txt.gz", {"PMDS", "count"});
	const auto N = _dPMDs_neg.size();
	for (size_t i = 0; i < N; ++i) {
		hist.writeln((N - i)/-10., _dPMDs_neg[N - 1 - i]);
	}
	for (size_t i = 0; i < _dPMDs_pos.size(); ++i) {
		hist.writeln(i/10., _dPMDs_pos[i]);
	}
};

}; // end namespace
