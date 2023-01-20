/*
 * TBamDiagnoser.cpp
 *
 *  Created on: May 30, 2020
 *      Author: phaentu
 */


#include "TBamDiagnoser.h"
#include <stdint.h>
#include <ostream>
#include "TBamFile.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "TCigar.h"
#include "coretools/Files/TFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TReadGroups.h"

namespace GenomeTasks{

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::TCountDistributionVector;

TBamDiagnoser::TBamDiagnoser():TGenome_filtered(){
    // initialize TGenome_basic stuff

	//settings
	_bamFile.readGroups().fillVectorWithNames(_readGroupNames);

};

void TBamDiagnoser::_writeHistogram(const TCountDistributionVector<> & distVec, const std::string& header, const std::string& name){
	//displays distributions of type 'TCountDistributionVector' as a histogram
	std::string filename = _outputName + "_" + header + "Histogram.txt";
	logfile().listFlush("Writing " + name + " histogram to '" + filename + "' ...");
	coretools::TOutputFile out(filename, {"readGroup", header, "count"});

    // Should file contain read groups with 0 counts?
	if (parameters().parameterExists("writeZeroCounts")) {
		distVec.template writeCombined<true>(out, "allReadGroups");
		distVec.template write<true>(out, _readGroupNames);
	} else {
		distVec.template writeCombined<false>(out, "allReadGroups");
		distVec.template write<false>(out, _readGroupNames);
	}

	out.close();
	logfile().done();
};

void TBamDiagnoser::_handleAlignment() {
    //get read group
    uint16_t readGroup = _bamFile.curReadGroupID();

    //increments for each read that passed filters
    _passedQC.add(readGroup);

    //add to counters
    _readLength.add(readGroup, _bamFile.curCIGAR().lengthRead());
    _usableLength.add(readGroup, _bamFile.curCIGAR().lengthAligned());
    _softClippedLength.add(readGroup, _bamFile.curCIGAR().lengthSoftClipped());
    _mappingQuality.add(readGroup, _bamFile.curMappingQuality());

    //fragment length: only for proper pairs and only once
    if(_bamFile.curIsProperPair() && !_bamFile.curIsReverseStrand())
        _fragmentLength.add(readGroup, _bamFile.curFragmentLength());
}

void TBamDiagnoser::diagnose(){
    //calculate length of genome
    double totLengthOfGenome = _bamFile.chromosomes().referenceLength();

    //initialize counters
    uint32_t numRG = _bamFile.readGroups().size();

    // resize distributions
    _passedQC.resize(numRG);
    _readLength.resize(numRG);
    _usableLength.resize(numRG);
    _softClippedLength.resize(numRG);
    _mappingQuality.resize(numRG);
    _fragmentLength.resize(numRG);

	//now parse through bam file
    _traverseBAMPassedQC();
    if(!parameters().parameterExists("splitMergeInput")){
    	logfile().list("Will not create input file for splitMerge. (use 'splitMergeInput' to do so).");
    }
	logfile().list("Approximate sequencing depth was estimated at ", (double) _usableLength.sum() / (double) totLengthOfGenome, ".");

	//writing output files
	logfile().startIndent("Writing output files:");

	//writing read group summary
	std::string filename = _outputName + "_diagnostics.txt";
	logfile().listFlush("Writing general diagnostics to '" + filename + "' ...");
	coretools::TOutputFile out(filename, {"readGroup", "totalReads", "passedQC", "duplicates", "avgReadLength", "seqCycles", "properPairs", "avgFragmentLength", "softClipped", "avgSoftClippedLength", "avgUsableAlignedLength", "approximateDepth", "avgMappingQuality", "seqType"});

	//determine sequencing type of BAM file
	uint32_t paired_count = 0;
	uint32_t single_count = 0;
	for(uint32_t rg = 0; rg < numRG; ++rg){
		if (_fragmentLength[rg].counts() == 0){
			++single_count;
		} else {
			++paired_count;
		}
	}

	//write for combined
	out.write("allReadGroups", (_bamFile.numAlignmentsReadPerReadGroup()).counts(), _passedQC.counts(), _bamFile.getDuplicateFilter().getCombinedCounts(),
			  _readLength.mean(), _readLength.max(), _fragmentLength.counts(), _fragmentLength.mean(),
			  _softClippedLength.countsLargerZero(), _softClippedLength.mean(), _usableLength.mean(),
			  (double)_usableLength.sum() / (double)totLengthOfGenome, _mappingQuality.mean());
	if (numRG == single_count) {
		out.writeln("single");
	} else if (numRG == paired_count) {
		out.writeln("paired");
	} else {
		out.writeln("mixed");
	}

	//write per read group
	for(uint32_t rg = 0; rg < numRG; ++rg){
		out.write(_bamFile.readGroups().getName(rg), (_bamFile.numAlignmentsReadPerReadGroup())[rg], _passedQC[rg], _bamFile.getDuplicateFilter().getCounts(rg),
				  _readLength[rg].mean(), _readLength[rg].max(), _fragmentLength[rg].counts(), _fragmentLength.mean(),
				  _softClippedLength[rg].countsLargerZero(), _softClippedLength[rg].mean(), _usableLength[rg].mean(),
				  (double)_usableLength[rg].sum() / (double)totLengthOfGenome, _mappingQuality[rg].mean());
		if (_fragmentLength[rg].counts() == 0) {
			out.writeln("single");
		} else {
			out.writeln("paired");
		}
	}
	out.close();
	logfile().done();

	if(parameters().parameterExists("splitMergeInput")){
	//write file used by split merge
	std::string splitmergename = _outputName + "_splitMergeInput.txt";
	logfile().listFlush("Outputting input file for splitMerge to '" + splitmergename + "' ...");
	coretools::TOutputFile splitm (splitmergename, {"readGroup", "seqType", "seqCycles"});
	for(uint32_t rg = 0; rg < numRG; ++rg){
		splitm << _bamFile.readGroups().getName(rg);
		if (_fragmentLength[rg].counts() == 0){
			splitm << "single";
		} else {
			splitm << "paired";
		}
		splitm.writeln(_readLength.max());
	}
	splitm.close();
	logfile().done();
	}


	//writing distributions
	_writeHistogram(_readLength, "readLength", "read length");
	_writeHistogram(_usableLength, "alignedLength", "aligned length");
	_writeHistogram(_softClippedLength, "softClippedLength", "soft clipped length");
	_writeHistogram(_fragmentLength, "fragmentLength", "fragment length");
	_writeHistogram(_mappingQuality, "mappingQuality", "mapping quality");

    logfile().endIndent(); //end writing output files
};

}; // end namespace
