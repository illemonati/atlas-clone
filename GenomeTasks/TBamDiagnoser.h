/*
 * TBamDiagnoser.h
 *
 *  Created on: May 30, 2020
 *      Author: phaentu
 */

#ifndef TBAMDIAGNOSER_H_
#define TBAMDIAGNOSER_H_

#include <string>
#include <vector>

#include "TBamFilter.h"
#include "TGenome.h"
#include "TTask.h"
#include "counters.h"

namespace GenomeTasks{

//-------------------------------------------
// TBamDiagnoser
// A class to get basic characteristics of a BAM file
//-------------------------------------------

class TBamDiagnoser:public TGenome_filtered{
private:
	BAM::TQualityFilter _qualFilter;
	std::vector<std::string> _readGroupNames;

    // distributions
    coretools::TCountDistribution _totalReads;
    coretools::TCountDistribution _passedQC;
    coretools::TCountDistributionVector _readLength;
    coretools::TCountDistributionVector _usableLength;
    coretools::TCountDistributionVector _softClippedLength;
    coretools::TCountDistributionVector _mappingQuality;
    coretools::TCountDistributionVector _fragmentLength;

	void _writeHistogram(const coretools::TCountDistributionVector & distVec, const std::string& header, const std::string& name);
    void _handleAlignment() override;

public:
	TBamDiagnoser();

	void diagnose();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_BAMDiagnostics:public coretools::TTask{
public:
	TTask_BAMDiagnostics(){ _explanation = "Estimating approximate depth, read length frequencies and mapping quality frequencies"; };

	void run(){
		using namespace coretools::instances;
		TBamDiagnoser diagnoser;
		diagnoser.diagnose();
	};
};


}; // end namespace

#endif /* TBAMDIAGNOSER_H_ */
