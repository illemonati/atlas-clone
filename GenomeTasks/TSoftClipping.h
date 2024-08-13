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
#include "TBamTraverser.h"
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
class TAssessSoftClipping final : public TBamReadTraverser<ReadType::Filtered> {
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
class TSoftClipsTrimmer final : public TBamReadTraverser<ReadType::Parsed> {
private:
	BAM::TOutputBamFile _outBam;
	void _handleAlignment(BAM::TAlignment& alignment) override;

public:
	TSoftClipsTrimmer();
	void run();
};

}; // namespace GenomeTasks

#endif /* TSOFTCLIPPING_H_ */
