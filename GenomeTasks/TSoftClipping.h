/*
 * TSoftClipping.h
 *
 *  Created on: Dec 9, 2019
 *      Author: wegmannd
 */

#ifndef TSOFTCLIPPING_H_
#define TSOFTCLIPPING_H_

#include <string>

#include "coretools/Files/TFile.h"
#include "coretools/Main/TTask.h"
#include "coretools/Math/counters.h"

#include "TBamFile.h"
#include "TGenome_OLD.h"
#include "TOutputBamFile.h"

namespace GenomeTasks {

// TODO: does using "left" and "right" make sense? Should we not rather use 3' and 5' ends?
// YES, we should use 3' and 5', please change.

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
class TAssessSoftClipping : public old::TGenome_filtered {
private:
	bool _writeAlignments = false;
	bool _printAll        = false;

	coretools::TCountDistributionVector<> left, right, total;
	TSoftClippingStatsFile statFile;

	void _handleAlignment() override;

public:
	TAssessSoftClipping();
	void run();
};

//--------------------------------------------------------
// TRemoveSoftClippedBases
//--------------------------------------------------------
class TRemoveSoftClippedBases : public old::TGenome_parsed {
private:
	BAM::TOutputBamFile _outBam;
	void _handleAlignment(BAM::TAlignment& alignment) override;

public:
	TRemoveSoftClippedBases();
	void run();
};

}; // namespace GenomeTasks

#endif /* TSOFTCLIPPING_H_ */
