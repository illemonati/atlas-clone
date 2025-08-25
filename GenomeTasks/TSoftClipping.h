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
#include "TOutputBamFile.h"
#include "TReadTraverser.h"

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
	void open(const std::string &Filename, bool PrintSequences);
	void write(const BAM::TBamFile &bamFile);
};

//--------------------------------------------------------
// TAssessSoftClipping
//--------------------------------------------------------
class TAssessSoftClipping {
private:
	TReadTraverser _readTraverser;
	bool _writeAlignments = false;
	bool _printAll        = false;

	coretools::TCountDistributionVector<> left, right, total;
	TSoftClippingStatsFile statFile;

	void _traverseReads();

public:
	TAssessSoftClipping();
	void run();
};

}; // namespace GenomeTasks

#endif /* TSOFTCLIPPING_H_ */
