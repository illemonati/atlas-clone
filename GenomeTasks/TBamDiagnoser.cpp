/*
 * TBamDiagnoser.cpp
 *
 *  Created on: May 30, 2020
 *      Author: phaentu
 */


#include "TBamDiagnoser.h"

namespace GenomeTasks{

TBamDiagnoser::TBamDiagnoser(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_basic(Parameters, Logfile, RandomGenerator){
	//settings
	qualFilter.set(Parameters, _logfile);
	_bamFile.readGroups().fillVectorWithNames(_readGroupNames);
};

void TBamDiagnoser::_writeHistogram(const TCountDistributionVector & distVec, const std::string header, const std::string name){
	// 1) read length
	std::string filename = _outputName + "_" + header + "Histogram.txt";
	_logfile->listFlush("Writing " + name + " histogram to '" + filename + "' ...");
	TOutputFile out(filename, {"readGroup", header, "count"});

	distVec.writeCombined(out, "allReadGroups");
	distVec.write(out, _readGroupNames);

	out.close();
	_logfile->done();
};

void TBamDiagnoser::diagnose(){
    //calculate length of genome
    double totLengthOfGenome = (double) _bamFile.chromosomes().referenceLength();

    //initialize counters
    uint32_t numRG = _bamFile.readGroups().size();

    TCountDistribution totalReads(numRG);
    TCountDistribution passedQC(numRG);
    TCountDistributionVector readLength(numRG);
    TCountDistributionVector usableLength(numRG);
    TCountDistributionVector softClippedLength(numRG);
    TCountDistributionVector mappingQuality(numRG);
    TCountDistributionVector fragmentLength(numRG);

	//now parse through bam file
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignment()){
		//get read group
		uint16_t readGroup = _bamFile.curReadGroupID();
		totalReads.add(readGroup);

		//passed filters?
		if(_bamFile.curPassedQC()){
			passedQC.add(readGroup);

			//add to counters
			readLength.add(readGroup, _bamFile.curCIGAR().lengthRead());
			usableLength.add(readGroup, _bamFile.curUsableAlignedLength(qualFilter));
			softClippedLength.add(readGroup, _bamFile.curCIGAR().lengthSoftClipped());
			mappingQuality.add(readGroup, _bamFile.curMappingQuality());

			//fragment length: only for proper pairs and only once
			if(_bamFile.curIsProperPair()){
				if(!_bamFile.curIsReverseStrand()){
					fragmentLength.add(readGroup, _bamFile.curFragmentLength());
				}
			}
		}

        //report progress
        _bamFile.printProgress();
    }

	//report end
	_bamFile.printEndWithSummary();

	_logfile->list("Approximate sequencing depth was estimated at " + toString((double) usableLength.sum() / (double) totLengthOfGenome));

	//writing output files
	_logfile->startIndent("Writing output files:");

	//writing read group summary
	std::string filename = _outputName + "_diagnostics.txt";
	_logfile->listFlush("Writing general diagnostics to '" + filename + "' ...");
	TOutputFile out(filename, {"readGroup", "reads", "passedQC", "fracPassedQC", "avgReadLength", "maxReadLength", "properPairs", "avgFragmentLength", "softClipped", "avgSoftClippedLength", "avgUsableAligneLength", "approximatDepth", "avgMappingQuality"});

	//write for combined
	out << "allReadGroups";
	out << totalReads.counts() << passedQC.counts() << (double) passedQC.counts() / (double) totalReads.counts();
	out << readLength.mean() << readLength.max();
	out << fragmentLength.counts() << fragmentLength.mean();
	out << softClippedLength.countsLargerZero() << softClippedLength.mean();
	out << usableLength.mean() << (double) usableLength.sum() / (double) totLengthOfGenome;
	out << mappingQuality.mean() << std::endl;

	//write per read group
	for(uint16_t rg = 0; rg < numRG; ++rg){
		out << _bamFile.readGroups().getName(rg);
		out << totalReads[rg] << passedQC[rg] << (double) passedQC[rg] / (double) totalReads[rg];
		out << readLength.mean() << readLength.max();
		out << fragmentLength[rg].counts() << fragmentLength.mean();
		out << softClippedLength[rg].countsLargerZero() << softClippedLength[rg].mean();
		out << usableLength[rg].mean() << (double) usableLength[rg].sum() / (double) totLengthOfGenome;
		out << mappingQuality[rg].mean() << std::endl;
	}
	out.close();
	_logfile->done();

	//TODO: write file used by split merge

	//writing distributions
	_writeHistogram(readLength, "readLength", "read length");
	_writeHistogram(usableLength, "usableLength", "usable length");
	_writeHistogram(softClippedLength, "softClippedLength", "soft clipped length");
	_writeHistogram(fragmentLength, "fragmentLength", "fragment length");
	_writeHistogram(mappingQuality, "mappingQuality", "mapping quality");

    _logfile->endIndent(); //end writing output files
};

}; // end namespace
