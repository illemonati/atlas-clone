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
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TTask.h"
#include "counters.h"

namespace GenomeTasks{

//TODO: does using "left" and "right" make sense? Should we not rather use 3' and 5' ends?

//--------------------------------------------------------
// TSoftClippingStatsFile
//--------------------------------------------------------
class TSoftClippingStatsFile{
private:
	coretools::TOutputFile _out;
	bool _printSequences;

public:
	TSoftClippingStatsFile(){ _printSequences = false; };
	TSoftClippingStatsFile(const std::string filename, const bool PrintSequences);
	void open(const std::string filename, const bool PrintSequences);
	void write(const BAM::TBamFile & bamFile);
};

//--------------------------------------------------------
// TAssessSoftClipping
//--------------------------------------------------------
class TAssessSoftClipping:public TGenome_filtered{
private:
	bool _writeAlignments;
	bool _printAll;

	coretools::TCountDistributionVector left, right, total;
	TSoftClippingStatsFile statFile;

	void _handleAlignment();

public:
	TAssessSoftClipping(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void assess();
};

//--------------------------------------------------------
// TRemoveSoftClippedBases
//--------------------------------------------------------
class TRemoveSoftClippedBases:public TGenome_parsed{
private:
	BAM::TOutputBamFile _outBam;
	void _handleAlignment();
public:
	TRemoveSoftClippedBases(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void removeSoftclippedBases();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_assessSoftClipping:public coretools::TTask{
public:
	TTask_assessSoftClipping(){ _explanation = "Assessing level of soft clipping in BAM file"; };

	void run(){
		using namespace coretools::instances;
		TAssessSoftClipping assessor(parameters(), &logfile(), &randomGenerator());
		assessor.assess();
	};
};

class TTask_removeSoftClippedBasesFromReads:public coretools::TTask{
public:
	TTask_removeSoftClippedBasesFromReads(){ _explanation = "Removing soft clipped bases from reads"; };

	void run(){
		using namespace coretools::instances;
		TRemoveSoftClippedBases remover(parameters(), &logfile(), &randomGenerator());
		remover.removeSoftclippedBases();
	};
};

}; // end namespace

#endif /* TSOFTCLIPPING_H_ */
