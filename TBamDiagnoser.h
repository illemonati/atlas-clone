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

//-------------------------------------------
// TBamDiagnoser
// A class to get basic characteristics of a BAM file
//-------------------------------------------

class TBamDiagnoser:public TGenome_basic{
private:
	TQualityFilter qualFilter;
	std::vector<std::string> _readGroupNames;

	void _writeHistogram(const TCountDistributionVector & distVec, const std::string header, const std::string name);

public:
	TBamDiagnoser(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);

	void diagnose();

};



#endif /* TBAMDIAGNOSER_H_ */
