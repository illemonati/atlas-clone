/*
 * TSoftClipping.cpp
 *
 *  Created on: Dec 12, 2019
 *      Author: schulzi
 */

#include "TSoftClipping.h"
#include <algorithm>
#include <ostream>
#include <vector>

#include "TAlignment.h"
#include "TChromosomes.h"
#include "TCigar.h"
#include "TGenomePosition.h"
#include "TLog.h"
#include "TParameters.h"

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
	_out << bamFile.curName() << bamFile.curChromosome().name << bamFile.curPosition().position();
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
	_out << std::endl;
};

//--------------------------------------------------------
// TAssessSoftClipping
//--------------------------------------------------------
TAssessSoftClipping::TAssessSoftClipping() : TGenome_filtered() {
	// limit input / output
	if (parameters().parameterExists("writeReads")) {
		_writeAlignments     = true;
		std::string filename = _outputName + "_softClippingStats.txt.gz";
		logfile().list("Will write alignments with softclipping to file '" + filename + "'. (parameter 'writeReads')");

		// write all reads?
		if (parameters().parameterExists("printAll")) {
			_printAll = true;
			logfile().list("Writing soft clipping stats for all reads to file. (parameter 'printAll')");
		} else {
			logfile().list(
			    "Writing soft clipping stats of soft clipped reads to file. (use 'printAll' to write for all reads)");
		}

		bool printSequences = false;
		if (parameters().parameterExists("printSequences")) {
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
	const BAM::TCigar &cigar = _bamFile.curCIGAR();
	left.add(cigar.lengthRead(), cigar.lengthSoftClippedLeft());
	right.add(cigar.lengthRead(), cigar.lengthSoftClippedRight());
	total.add(cigar.lengthRead(), cigar.lengthSoftClippedLeft() + cigar.lengthSoftClippedRight());

	// write to file
	if (_writeAlignments && (cigar.lengthSoftClipped() > 0 || _printAll)) { statFile.write(_bamFile); }
};

void TAssessSoftClipping::assess() {
	// now parse through bam file and write alignments
	_traverseBAMPassedQC();

	// write counts
	logfile().startIndent("Writing soft clipping distributions:");
	std::string filename = _outputName + "_softClippingMatrixLeft.txt";
	logfile().listFlush("Writing distribution of soft clipping on left to file '" + filename + "' ...");
	left.write(filename, "readLength", "softClippedLengthLeft");
	logfile().done();

	filename = _outputName + "_softClippingMatrixRight.txt";
	logfile().listFlush("Writing distribution of soft clipping on right to file '" + filename + "' ...");
	left.write(filename, "readLength", "softClippedLengthRight");
	logfile().done();

	filename = _outputName + "_softClippingMatrixBoth.txt";
	logfile().listFlush("Writing distribution of soft clipping on both combined to file '" + filename + "' ...");
	left.write(filename, "readLength", "softClippedLengthBoth");
	logfile().done();
	logfile().endIndent();
};

//--------------------------------------------------------
// TRemoveSoftClippedBases
//--------------------------------------------------------
TRemoveSoftClippedBases::TRemoveSoftClippedBases() : TGenome_parsed(){};

void TRemoveSoftClippedBases::_handleAlignment() {
	if (_bamFile.curCIGAR().lengthSoftClipped() > 0) {
		// remove softclipped reads
		_alignment.removeSoftClippedBases();
		// write
		_outBam.writeAlignment(_alignment);
	} else {
		// just write current alignment
		_bamFile.writeCurAlignment(_outBam);
	}
};

void TRemoveSoftClippedBases::removeSoftClippedBases() {
	std::string filename = _outputName + "_softClippedBasesRemoved.bam";
	logfile().list("Writing reads after soft-clip trimming to file '" + filename + "'.");
	_openBamForWriting(filename, _outBam);

	// traverse BAM
	_traverseBAMPassedQC();

	// report
	_outBam.close(&logfile());
};

}; // namespace GenomeTasks
