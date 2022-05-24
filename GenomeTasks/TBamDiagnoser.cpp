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
#include "TChromosomes.h"
#include "TCigar.h"
#include "TFile.h"
#include "TLog.h"
#include "TRandomGenerator.h"
#include "TReadGroups.h"

namespace GenomeTasks{

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::TCountDistributionVector;

TBamDiagnoser::TBamDiagnoser():TGenome_filtered(){
    // initialize TGenome_basic stuff

	//settings
	_qualFilter.set(parameters(), &logfile());
	_bamFile.readGroups().fillVectorWithNames(_readGroupNames);

	// should we filter -> only call setFilters if we want to filter!
	// By default, all filters are off
	if (parameters().parameterExists("filterBAM")){
	    _bamFile.setFilters(parameters(), &logfile());
	} else {
	    logfile().list("Will keep all reads. (use 'filterBAM' to limit)");
	}
};

void TBamDiagnoser::_writeHistogram(const TCountDistributionVector & distVec, const std::string& header, const std::string& name){
	// 1) read length
	std::string filename = _outputName + "_" + header + "Histogram.txt";
	logfile().listFlush("Writing " + name + " histogram to '" + filename + "' ...");
	coretools::TOutputFile out(filename, {"readGroup", header, "count"});

	distVec.writeCombined(out, "allReadGroups");
	distVec.write(out, _readGroupNames);

	out.close();
	logfile().done();
};

void TBamDiagnoser::_handleAlignment() {
    //get read group
    uint16_t readGroup = _bamFile.curReadGroupID();
    _totalReads.add(readGroup);

    //passed filters?
    _passedQC.add(readGroup);

    if(_bamFile.curIsDuplicate()) {
    	_numDuplicates.add(readGroup);
    	}

    //add to counters
    _readLength.add(readGroup, _bamFile.curCIGAR().lengthRead());
    _usableLength.add(readGroup, _bamFile.curUsableAlignedLength(_qualFilter));
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
    _totalReads.resize(numRG);
    _passedQC.resize(numRG);
    _numDuplicates.resize(numRG);
    _readLength.resize(numRG);
    _usableLength.resize(numRG);
    _softClippedLength.resize(numRG);
    _mappingQuality.resize(numRG);
    _fragmentLength.resize(numRG);

	//now parse through bam file
    _traverseBAMPassedQC();
	logfile().list("Approximate sequencing depth was estimated at ", (double) _usableLength.sum() / (double) totLengthOfGenome, ".");

	//writing output files
	logfile().startIndent("Writing output files:");

	//writing read group summary
	std::string filename = _outputName + "_diagnostics.txt";
	logfile().listFlush("Writing general diagnostics to '" + filename + "' ...");
	coretools::TOutputFile out(filename, {"readGroup", "reads", "passedQC", "fracPassedQC", "avgReadLength", "maxReadLength", "properPairs", "avgFragmentLength", "softClipped", "avgSoftClippedLength", "avgUsableAlignedLength", "approximateDepth", "avgMappingQuality", "numDuplicates"});

	//write for combined
	out << "allReadGroups";
	out << _totalReads.counts() << _passedQC.counts() << (double) _passedQC.counts() / (double) _totalReads.counts();
	out << _readLength.mean() << _readLength.max();
	out << _fragmentLength.counts() << _fragmentLength.mean();
	out << _softClippedLength.countsLargerZero() << _softClippedLength.mean();
	out << _usableLength.mean() << (double) _usableLength.sum() / (double) totLengthOfGenome;
	out << _mappingQuality.mean();
	out << _numDuplicates.counts() << std::endl;

	//write per read group
	for(uint32_t rg = 0; rg < numRG; ++rg){
		out << _bamFile.readGroups().getName(rg);
		out << _totalReads[rg] << _passedQC[rg] << (double) _passedQC[rg] / (double) _totalReads[rg];
		out << _readLength.mean() << _readLength.max();
		out << _fragmentLength[rg].counts() << _fragmentLength.mean();
		out << _softClippedLength[rg].countsLargerZero() << _softClippedLength[rg].mean();
		out << _usableLength[rg].mean() << (double) _usableLength[rg].sum() / (double) totLengthOfGenome;
		out << _mappingQuality[rg].mean();
		out << _numDuplicates[rg] << std::endl;
	}
	out.close();
	logfile().done();

	//TODO: write file used by split merge

	//writing distributions
	_writeHistogram(_readLength, "readLength", "read length");
	_writeHistogram(_usableLength, "usableLength", "usable length");
	_writeHistogram(_softClippedLength, "softClippedLength", "soft clipped length");
	_writeHistogram(_fragmentLength, "fragmentLength", "fragment length");
	_writeHistogram(_mappingQuality, "mappingQuality", "mapping quality");

    logfile().endIndent(); //end writing output files
};

}; // end namespace
