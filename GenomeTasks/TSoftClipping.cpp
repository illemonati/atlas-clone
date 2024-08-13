/*
 * TSoftClipping.cpp
 *
 *  Created on: Dec 12, 2019
 *      Author: schulzi
 */

#include "TSoftClipping.h"
#include <vector>

#include "TAlignment.h"
#include "TBamTraverser.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "TCigar.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;

//------------------------------------------
// TSoftClippingStatsFile
//------------------------------------------
TSoftClippingStatsFile::TSoftClippingStatsFile(const std::string &Filename, bool PrintSequences) {
	open(Filename, PrintSequences);
};

void TSoftClippingStatsFile::open(const std::string &Filename, bool PrintSequences) {
	_printSequences = PrintSequences;

	// open output file
	_out.open(Filename);

	// write header
	std::vector<std::string> header = {"Read", "chromosome", "position", "nClippedLeft", "nNotClipped", "nClippedRight"};
	if (_printSequences) {
		header.emplace_back("sClippedLeft");
		header.emplace_back("sNotClipped");
		header.emplace_back("sClippedRight");
	}

	_out.writeHeader(header);
};

void TSoftClippingStatsFile::write(const BAM::TBamFile &bamFile) {
	_out << bamFile.curName() << bamFile.curChromosome().name() << bamFile.curPosition().position();
	const BAM::TCigar &cigar = bamFile.curCIGAR();

	_out << cigar.lengthSoftClippedLeft() << cigar.lengthSequenced() << cigar.lengthSoftClippedRight();
	if (_printSequences) {
		// left
		if (cigar.lengthSoftClippedLeft() > 0) {
			_out << bamFile.curQuerySequence(0, cigar.lengthSoftClippedLeft());
		} else {
			_out << "";
		}

		// middle
		_out << bamFile.curQuerySequence(cigar.lengthSoftClippedLeft(), cigar.lengthSequenced());

		// right
		if (cigar.lengthSoftClippedRight() > 0) {
			_out << bamFile.curQuerySequence(cigar.lengthSoftClippedLeft() + cigar.lengthSequenced(),
			                                 cigar.lengthSoftClippedRight());
		} else {
			_out << "";
		}
	}
	_out.endln();
};

//--------------------------------------------------------
// TAssessSoftClipping
//--------------------------------------------------------
TAssessSoftClipping::TAssessSoftClipping() : TBamReadTraverser<ReadType::Filtered>() {
	// limit input / output
	if (parameters().exists("writeReads")) {
		_writeAlignments     = true;
		const auto filename = _genome.outputName() + "_softClippingStats.txt.gz";
		logfile().list("Will write alignments with softclipping to file '" + filename + "'. (parameter 'writeReads')");

		// write all reads?
		if (parameters().exists("printAll")) {
			_printAll = true;
			logfile().list("Writing soft clipping stats for all reads to file. (parameter 'printAll')");
		} else {
			logfile().list(
			    "Writing soft clipping stats of soft clipped reads to file. (use 'printAll' to write for all reads)");
		}

		bool printSequences = false;
		if (parameters().exists("printSequences")) {
			printSequences = true;
			logfile().list("Writing soft clipped bases to file. (parameter 'printSequences')");
		} else {
			logfile().list(
			    "Writing only counts of soft clipped bases to file. (use 'printSequences' to also print sequences)");
		}

		// open file
		statFile.open(filename, printSequences);
	}
};

void TAssessSoftClipping::_handleAlignment() {
	// add to counters
	const BAM::TCigar &cigar = _genome.bamFile().curCIGAR();
	left.add(cigar.lengthRead(), cigar.lengthSoftClippedLeft());
	right.add(cigar.lengthRead(), cigar.lengthSoftClippedRight());
	total.add(cigar.lengthRead(), cigar.lengthSoftClippedLeft() + cigar.lengthSoftClippedRight());

	// write to file
	if (_writeAlignments && (cigar.lengthSoftClipped() > 0 || _printAll)) { statFile.write(_genome.bamFile()); }
};

void TAssessSoftClipping::run() {
	// now parse through bam file and write alignments
	_traverseBAMPassedQC();

	// write counts
	logfile().startIndent("Writing soft clipping distributions:");
	std::string filename = _genome.outputName() + "_softClippingMatrixLeft.txt";
	logfile().listFlush("Writing distribution of soft clipping on left to file '" + filename + "' ...");
	left.write(filename, "readLength", "softClippedLengthLeft");
	logfile().done();

	filename = _genome.outputName() + "_softClippingMatrixRight.txt";
	logfile().listFlush("Writing distribution of soft clipping on right to file '" + filename + "' ...");
	left.write(filename, "readLength", "softClippedLengthRight");
	logfile().done();

	filename = _genome.outputName() + "_softClippingMatrixBoth.txt";
	logfile().listFlush("Writing distribution of soft clipping on both combined to file '" + filename + "' ...");
	left.write(filename, "readLength", "softClippedLengthBoth");
	logfile().done();
	logfile().endIndent();
};

//--------------------------------------------------------
// TRemoveSoftClippedBases
//--------------------------------------------------------
	TSoftClipsTrimmer::TSoftClipsTrimmer() : TBamReadTraverser<ReadType::Parsed>(), _outBam(_genome.outputName() + "_softClippedBasesRemoved.bam", _genome.bamFile()){
		
	};

void TSoftClipsTrimmer::_handleAlignment(BAM::TAlignment& alignment) {
	alignment.trimSoftClips();
	_outBam.writeAlignment(alignment);
};

void TSoftClipsTrimmer::run() {
	logfile().list("Writing reads after soft-clip trimming to file '", _genome.outputName(),  "_softClippedBasesRemoved.bam'.");

	// traverse BAM
	_traverseBAMPassedQC();
};

}; // namespace GenomeTasks
