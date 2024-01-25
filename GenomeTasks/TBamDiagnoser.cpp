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
#include "TBamFilters.h"
#include "coretools/Main/TParameters.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "TCigar.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TReadGroups.h"
#include "genometools/GenomePositions/TGenomePosition.h"

namespace GenomeTasks{

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::TCountDistributionVector;

;

void TBamDiagnoser::_writeHistogram(const std::vector<TCountDistributionVector<>> & distVec, const std::string& header, const std::string& name){
	//displays distributions of type 'TCountDistributionVector' as a histogram
	std::string filename = _genome.outputName() + "_" + header + "Histogram.txt";
	logfile().listFlush("Writing " + name + " histogram to '" + filename + "' ...");
	coretools::TOutputFile out(filename, {"readGroup", header, "count"});

	TCountDistributionVector<> distributionPerReadGroup;
	coretools::TCountDistribution<> temp;
	for (size_t i = 0; i < distVec.size(); i++){
		distVec[i].fillCombinedDistribution(temp);
		distributionPerReadGroup.add(i, temp);
	}

    // Should file contain read groups with 0 counts?
	if (parameters().exists("writeZeroCounts")) {
		distributionPerReadGroup.template writeCombined<true>(out, "allReadGroups");
		distributionPerReadGroup.template write<true>(out, _readGroupNames);
	} else {
		distributionPerReadGroup.template writeCombined<false>(out, "allReadGroups");
		distributionPerReadGroup.template write<false>(out, _readGroupNames);
	}

	out.close();
	logfile().done();
};

void TBamDiagnoser::_writeTable(const TCountDistributionVector<> & distVec, const std::string& header, const std::string &name){
	//displays distributions of type 'TCountDistributionVector' as a table
	std::string filename = _genome.outputName() + "_" + header + "Table.txt";
	logfile().listFlush("Writing " + name + " table to '" + filename + "' ...");
	std::vector<std::string> chromNames;
	auto chromIt = _genome.bamFile().chromosomes().cbegin();
	while(chromIt != _genome.bamFile().chromosomes().cend()){
		chromNames.push_back(chromIt->name());
		chromIt++;
	}
	std::vector<std::string> readGroupNames = _readGroupNames;
	readGroupNames.insert(readGroupNames.begin(), "allReadGroups");
	distVec.writeAsMatrixCombined(filename, header, chromNames, readGroupNames);
}

void TBamDiagnoser::_handleAlignment() {
    //get read group
	const size_t readGroup  = _genome.bamFile().curReadGroupID();
	if (readGroup == BAM::TReadGroups::noReadGroupId) return;

	const size_t chromosome = _genome.bamFile().refID();

	//increments for each read that passed filters
    _passedQC.add(readGroup, _genome.bamFile().curChromosome().refID());

    //add to counters
	const auto& curPosition = _genome.bamFile().curPosition();
	if (curPosition.refID() == _oldPositions[readGroup].refID()) {
		_readDist[readGroup].add(chromosome, _genome.bamFile().curPosition() - _oldPositions[readGroup]);
	}
	_oldPositions[readGroup] = curPosition;
	_readLength[readGroup].add(chromosome, _genome.bamFile().curCIGAR().lengthRead());
	_usableLength[readGroup].add(chromosome, _genome.bamFile().curCIGAR().lengthAligned());
	_softClippedLength[readGroup].add(chromosome, _genome.bamFile().curCIGAR().lengthSoftClipped());
	_mappingQuality[readGroup].add(chromosome, _genome.bamFile().curMappingQuality());

    //fragment length: only for proper pairs and only once
    if(_genome.bamFile().curIsProperPair() && !_genome.bamFile().curIsReverseStrand())
		_fragmentLength[readGroup].add(chromosome, _genome.bamFile().curFragmentLength());
}

void TBamDiagnoser::run(){
    //calculate length of genome
    double totLengthOfGenome = _genome.bamFile().chromosomes().referenceLength();

    //initialize counters
	const size_t numRG    = _genome.bamFile().readGroups().size();
	const size_t numChrom = _genome.bamFile().chromosomes().size();

	// resize distributions
    _passedQC.resize(numRG);
	_passedQC.resizeDistributions(numChrom);

	constexpr size_t BIG = -1;
	_oldPositions.resize(numRG, {BIG, BIG});
    _readLength.resize(numRG);
	_readDist.resize(numRG);
    _usableLength.resize(numRG);
    _softClippedLength.resize(numRG);
    _mappingQuality.resize(numRG);
    _fragmentLength.resize(numRG);
	for (size_t i = 0; i < numRG; i++) {
		_readLength[i].resize(numChrom);
		_readDist[i].resize(numChrom);
		_usableLength[i].resize(numChrom);
		_softClippedLength[i].resize(numChrom);
		_mappingQuality[i].resize(numChrom);
		_fragmentLength[i].resize(numChrom);
	}


	if(parameters().exists("perChromosome")){
		_chromStats = true;
		logfile().list("Will output data per chromosome into diagnostics-file. (parameter 'perChromosome')");
	} else {
		logfile().list("Will not output data per chromosome into diagnostics-file. (parameter 'perChromosome')");
	}

	//now parse through bam file
    _traverseBAMPassedQC();
    if(!parameters().exists("splitMergeInput")){
    	logfile().list("Will not create input file for splitMerge. (use 'splitMergeInput' to do so).");
    }
	if(!parameters().exists("printReferenceLength")){
		logfile().list("Will not print reference lengths of chromosomes to file. (use 'printReferenceLength' to do so).");
	}
	logfile().list("Approximate sequencing depth was estimated at ", (double) sumOverAllReadGroups(_usableLength) / (double) totLengthOfGenome, ".");

	//writing output files
	logfile().startIndent("Writing output files:");

	//writing read group summary
	std::string filename = _genome.outputName() + "_diagnostics.txt";
	logfile().listFlush("Writing general diagnostics to '" + filename + "' ...");
	coretools::TOutputFile out(filename);

	std::vector<std::string_view> header{"readGroup"};
	if (_chromStats) header.push_back("chromosome");
	header.insert(header.end(), {"totalReads", "passedQC", "duplicates", "avgReadLength", "seqCycles", "avgReadDist",
								 "properPairs", "avgFragmentLength", "softClipped", "avgSoftClippedLength",
								 "avgUsableAlignedLength", "approximateDepth", "avgMappingQuality", "seqType"});
	out.writeHeader(header);

	//determine sequencing type of BAM file
	size_t paired_count = 0;
	size_t single_count = 0;
	for(size_t rg = 0; rg < numRG; ++rg){
		if (_fragmentLength[rg].counts() == 0){
			++single_count;
		} else {
			++paired_count;
		}
	}
	//write for all read groups and all chromosomes
	out.write("allReadGroups");
	if (_chromStats) {
		out.write("allChromosomes");
	}
	out.write(_genome.bamFile().numAlignmentReadPerReadGroupPerChromosome().counts(), _passedQC.counts(),
	          _genome.bamFile().filter(BAM::FilterType::Duplicate).getCombinedCounts(),
	          meanOverAllReadGroups(_readLength), maxOverAllReadGroups(_readLength), meanOverAllReadGroups(_readDist),
	          countsOverAllReadGroups(_fragmentLength), meanOverAllReadGroups(_fragmentLength),
			  countsLargerZeroOverAllReadGroups(_softClippedLength), meanOverAllReadGroups(_softClippedLength),
			  meanOverAllReadGroups(_usableLength), (double)sumOverAllReadGroups(_usableLength) / totLengthOfGenome,
	          meanOverAllReadGroups(_mappingQuality));
	if (numRG == single_count) {
		out.writeln("single");
	} else if (numRG == paired_count) {
		out.writeln("paired");
	} else {
		out.writeln("mixed");
	}

	//write for all read groups per chromosome
	if(_chromStats){
		auto chromIt = _genome.bamFile().chromosomes().cbegin();
		while(chromIt != _genome.bamFile().chromosomes().cend()){
			size_t refID = chromIt->refID();
			out.write("allReadGroups", chromIt->name(),
					  (_genome.bamFile().numAlignmentReadPerReadGroupPerChromosome()).horizontalCounts(refID),
					  _passedQC.horizontalCounts(refID),
					  _genome.bamFile().filter(BAM::FilterType::Duplicate).getCountsPerChromosome(refID),
					  meanForChromosome(_readLength, refID), maxForChromosome(_readLength, refID), meanForChromosome(_readDist, refID),
					  countsForChromosome(_fragmentLength, refID), meanForChromosome(_fragmentLength, refID),
					  countsLargerZeroForChromosome(_softClippedLength, refID),
					  meanForChromosome(_softClippedLength, refID), meanForChromosome(_usableLength, refID),
					  (double)sumForChromosome(_usableLength, refID) / (double)chromIt->length(),
			          meanForChromosome(_mappingQuality, refID));
			if (numRG == single_count) {
				out.writeln("single");
			} else if (numRG == paired_count) {
				out.writeln("paired");
			} else {
				out.writeln("mixed");
			}
			chromIt++;
		}
	}

	//write per read group for all chromosomes
	for(size_t rg = 0; rg < numRG; ++rg){
		out.write(_genome.bamFile().readGroups().getName(rg));
		if (_chromStats) out.write("allChromosomes");
		out.write((_genome.bamFile().numAlignmentReadPerReadGroupPerChromosome())[rg].counts(), _passedQC[rg].counts(),
				  _genome.bamFile().filter(BAM::FilterType::Duplicate).getCounts(rg), _readLength[rg].mean(),
				  _readLength[rg].max(), _readDist[rg].mean(), _fragmentLength[rg].counts(), _fragmentLength[rg].mean(),
				  _softClippedLength[rg].countsLargerZero(), _softClippedLength[rg].mean(), _usableLength[rg].mean(),
				  (double)_usableLength[rg].sum() / (double)totLengthOfGenome, _mappingQuality[rg].mean());
		if (_fragmentLength[rg].counts() == 0) {
			out.writeln("single");
		} else {
			out.writeln("paired");
		}
		if(_chromStats){
			// write per read group per chromosome
			auto chromIt = _genome.bamFile().chromosomes().cbegin();
			while(chromIt != _genome.bamFile().chromosomes().cend()){
				size_t refID = chromIt->refID();
				out.write(
					_genome.bamFile().readGroups().getName(rg), chromIt->name(),
					(_genome.bamFile().numAlignmentReadPerReadGroupPerChromosome())[rg][refID], _passedQC[rg][refID],
					_genome.bamFile().filter(BAM::FilterType::Duplicate).getCountsAtReadGroupAndChromosome(rg, refID),
					_readLength[rg][refID].mean(), _readLength[rg][refID].max(), _readDist[rg][refID].mean(),
					_fragmentLength[rg][refID].counts(), _fragmentLength[rg][refID].mean(),
					_softClippedLength[rg][refID].countsLargerZero(), _softClippedLength[rg][refID].mean(),
					_usableLength[rg][refID].mean(), (double)_usableLength[rg][refID].sum() / (double)chromIt->length(),
					_mappingQuality[rg][refID].mean());
				if (_fragmentLength[rg][refID].counts() == 0) {
					out.writeln("single");
				} else {
					out.writeln("paired");
				}
				chromIt++;
			}
		}
	}
	out.close();
	logfile().done();

	if(parameters().exists("splitMergeInput")){
		//write file used by split merge
		std::string splitmergename = _genome.outputName() + "_splitMergeInput.txt";
		logfile().listFlush("Outputting input file for splitMerge to '" + splitmergename + "' ...");
		coretools::TOutputFile splitm (splitmergename, {"readGroup", "seqType", "seqCycles"});
		for(size_t rg = 0; rg < numRG; ++rg){
			splitm << _genome.bamFile().readGroups().getName(rg);
			if (_fragmentLength[rg].counts() == 0){
				splitm << "single";
			} else {
				splitm << "paired";
			}
			splitm.writeln(maxOverAllReadGroups(_readLength));
		}
		splitm.close();
		logfile().done();
	}

	if(parameters().exists("printReferenceLength")){
		//write file with length of all contigs
		std::string referenceLengthName = _genome.outputName() + "_referenceLengths.txt";
		logfile().listFlush("Outputting reference lengths of all contigs to '" + referenceLengthName + "' ...");
		coretools::TOutputFile refLen (referenceLengthName, {"chromosome", "length"});
		auto it = _genome.bamFile().chromosomes().cbegin();
		while(it != _genome.bamFile().chromosomes().cend()){
			refLen << it->name() << it->length() << coretools::endl;
			++it;
		}
		refLen.close();
		logfile().done();
	}

	//writing distributions
	_writeHistogram(_readLength, "readLength", "read length");
	_writeHistogram(_readDist, "readDist", "read length");
	_writeHistogram(_usableLength, "alignedLength", "aligned length");
	_writeHistogram(_softClippedLength, "softClippedLength", "soft clipped length");
	_writeHistogram(_fragmentLength, "fragmentLength", "fragment length");
	_writeHistogram(_mappingQuality, "mappingQuality", "mapping quality");

    logfile().endIndent(); //end writing output files
};

double TBamDiagnoser::meanOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec){
	size_t counts = countsOverAllReadGroups(vec);
	if (counts == 0) { return 0.0; }
	size_t sum = sumOverAllReadGroups(vec);
	return ((double) sum / (double) counts);
}

double TBamDiagnoser::meanForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID){
	size_t counts = countsForChromosome(vec, chromRefID);
	if (counts == 0) { return 0.0; }
	size_t sum = sumForChromosome(vec, chromRefID);

	return ((double) sum / (double) counts);
}

size_t TBamDiagnoser::maxOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec){
	size_t max = 0;
	for (const auto &s: vec){
		size_t curMax = s.max();
		if (curMax > max)
			max = curMax;
	}
	return max;
}

size_t TBamDiagnoser::maxForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID){
	size_t max = 0;
	for (const auto &s: vec){
		size_t curMax = s[chromRefID].max();
		if (curMax > max)
			max = curMax;
	}
	return max;
}

size_t TBamDiagnoser::countsOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec){
	size_t counts = 0;
	for (const auto &s: vec){
		counts += s.counts();
	}
	return counts;
}

size_t TBamDiagnoser::countsForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID){
	size_t counts = 0;
	for (const auto &s: vec){
		counts += s[chromRefID].counts();
	}
	return counts;
}


size_t TBamDiagnoser::countsLargerZeroOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec){
	size_t counts = 0;
	for (const auto &s: vec){
		counts += s.countsLargerZero();
	}
	return counts;
}

size_t TBamDiagnoser::countsLargerZeroForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID){
	size_t counts = 0;
	for (const auto &s: vec){
		counts += s[chromRefID].countsLargerZero();
	}
	return counts;
}

size_t TBamDiagnoser::sumOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec){
	size_t sum = 0;
	for (const auto &s: vec){
		sum += s.sum();
	}
	return sum;
}

size_t TBamDiagnoser::sumForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID){
	size_t sum = 0;
	for (const auto &s: vec){
		sum += s[chromRefID].sum();
	}
	return sum;
}

}; // end namespace
