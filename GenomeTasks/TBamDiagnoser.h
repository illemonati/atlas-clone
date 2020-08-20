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

//-------------------------------------------
// TBamDiagnoser
// A class to get basic characteristics of a BAM file
//-------------------------------------------

class TBamDiagnoser:public TGenome_basic{
private:
	BAM::TQualityFilter qualFilter;
	std::vector<std::string> _readGroupNames;

	void _writeHistogram(const TCountDistributionVector & distVec, const std::string header, const std::string name);

public:
	TBamDiagnoser(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);

	void diagnose();

};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_BAMDiagnostics:public TTask{
public:
	TTask_BAMDiagnostics(){ _explanation = "Estimating approximate depth, read length frequencies and mapping quality frequencies"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TBamDiagnoser diagnoser(Parameters, Logfile, _randomGenerator);
		diagnoser.diagnose();
	};
};


}; // end namespace

#endif /* TBAMDIAGNOSER_H_ */
