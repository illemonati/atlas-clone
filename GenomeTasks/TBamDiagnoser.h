/*
 * TBamDiagnoser.h
 *
 *  Created on: May 30, 2020
 *      Author: phaentu
 */

#ifndef TBAMDIAGNOSER_H_
#define TBAMDIAGNOSER_H_

#include "TGenome.h"
#include "TFile.h"
#include "TTask.h"

namespace GenomeTasks{

using coretools::TCountDistribution;
using coretools::TCountDistributionVector;

//-------------------------------------------
// TBamDiagnoser
// A class to get basic characteristics of a BAM file
//-------------------------------------------

class TBamDiagnoser:public TGenome_filtered{
private:
	BAM::TQualityFilter _qualFilter;
	std::vector<std::string> _readGroupNames;

    // distributions
    TCountDistribution _totalReads;
    TCountDistribution _passedQC;
    TCountDistributionVector _readLength;
    TCountDistributionVector _usableLength;
    TCountDistributionVector _softClippedLength;
    TCountDistributionVector _mappingQuality;
    TCountDistributionVector _fragmentLength;

	void _writeHistogram(const TCountDistributionVector & distVec, const std::string& header, const std::string& name);
    void _handleAlignment() override;

public:
	TBamDiagnoser(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);

	void diagnose();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_BAMDiagnostics:public coretools::TTask{
public:
	TTask_BAMDiagnostics(){ _explanation = "Estimating approximate depth, read length frequencies and mapping quality frequencies"; };

	void run(coretools::TParameters & Parameters, coretools::TLog* Logfile){
		TBamDiagnoser diagnoser(Parameters, Logfile, _randomGenerator);
		diagnoser.diagnose();
	};
};


}; // end namespace

#endif /* TBAMDIAGNOSER_H_ */
