/*
 * TCreateBedMask.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCREATEBEDMASK_H_
#define GENOMETASKS_TCREATEBEDMASK_H_

#include "TGenome.h"

namespace GenomeTasks{

//--------------------------------------
// TCreateBedMask
//--------------------------------------
class TCreateBedMask:public TGenome_windows{
protected:
	BAM::TBed _bed;
	uint32_t _minDepthForMask;
public:
	TCreateBedMask(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
class TCreateDepthBedMask:public TCreateBedMask{
private:
	uint32_t _maxDepthForMask;

	void _handleWindow();
public:
	TCreateDepthBedMask(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void createDepthMask();
};

//--------------------------------------
// TCreateConservedBedMask
//--------------------------------------
class TCreateConservedBedMask:public TCreateBedMask{
private:
	void _handleWindow();

	//tmp variables
	GenotypeLikelihoods::TBaseCounts _baseCounts;
public:
	TCreateConservedBedMask(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void createConservedMask();
};

}; // end namespace

#endif /* GENOMETASKS_TCREATEBEDMASK_H_ */
