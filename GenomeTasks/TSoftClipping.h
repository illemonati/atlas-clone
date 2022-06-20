/*
 * TSoftClipping.h
 *
 *  Created on: Dec 9, 2019
 *      Author: wegmannd
 */

#ifndef TSOFTCLIPPING_H_
#define TSOFTCLIPPING_H_

#include <string>

#include "TBamFile.h"
#include "TFile.h"
#include "TGenome.h"
#include "TTask.h"
#include "counters.h"

namespace GenomeTasks {

// TODO: does using "left" and "right" make sense? Should we not rather use 3' and 5' ends?

//--------------------------------------------------------
// TSoftClippingStatsFile
//--------------------------------------------------------
class TSoftClippingStatsFile {
private:
	coretools::TOutputFile _out;
	bool _printSequences = false;

public:
	TSoftClippingStatsFile() = default;
	TSoftClippingStatsFile(const std::string &Filename, bool PrintSequences);
	void open(const std::string &Filename, bool PrintSequences);
	void write(const BAM::TBamFile &bamFile);
};

//--------------------------------------------------------
// TAssessSoftClipping
//--------------------------------------------------------
class TAssessSoftClipping : public TGenome_filtered {
private:
	bool _writeAlignments = false;
	bool _printAll        = false;

	coretools::TCountDistributionVector<> left, right, total;
	TSoftClippingStatsFile statFile;

	void _handleAlignment() override;

public:
	TAssessSoftClipping();
	void assess();
};

//--------------------------------------------------------
// TRemoveSoftClippedBases
//--------------------------------------------------------
class TRemoveSoftClippedBases : public TGenome_parsed {
private:
	BAM::TOutputBamFile _outBam;
	void _handleAlignment() override;

public:
	TRemoveSoftClippedBases();
	void removeSoftClippedBases();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_assessSoftClipping : public coretools::TTask {
public:
	TTask_assessSoftClipping() { _explanation = "Assessing level of soft clipping in BAM file"; };

	void run() override {
		TAssessSoftClipping assessor;
		assessor.assess();
	};
};

class TTask_removeSoftClippedBasesFromReads : public coretools::TTask {
public:
	TTask_removeSoftClippedBasesFromReads() { _explanation = "Removing soft clipped bases from reads"; };

	void run() override {
		TRemoveSoftClippedBases remover;
		remover.removeSoftClippedBases();
	};
};

}; // namespace GenomeTasks

#endif /* TSOFTCLIPPING_H_ */
