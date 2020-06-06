/*
 * TCreateBedMask.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCREATEBEDMASK_H_
#define GENOMETASKS_TCREATEBEDMASK_H_

#include "TGenome.h"
#include "TTask.h"

namespace GenomeTasks{

//--------------------------------------
// TCreateBedMask
//--------------------------------------
class TCreateBedMask:public TGenome_windows{
protected:
	BAM::TBed _bed;
	uint32_t _minDepthForMask;
public:
	TCreateBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};

//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
class TCreateDepthBedMask:public TCreateBedMask{
private:
	uint32_t _maxDepthForMask;

	void _handleWindow();
public:
	TCreateDepthBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
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
	TCreateConservedBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void createConservedMask();
};


//--------------------------------------
// Tasks
//--------------------------------------
class TTask_createDepthMask:public TTask{
public:
	TTask_createDepthMask(){ _explanation = "Creating depth mask BED file"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TCreateDepthBedMask depthMask(Parameters, Logfile, _randomGenerator);
		depthMask.createDepthMask();
	};
};

class TTask_createNonRefMask:public TTask{
public:
	TTask_createNonRefMask(){ _explanation = "Creating mask BED file with sites with non-reference alleles"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TCreateConservedBedMask conservedMask(Parameters, Logfile, _randomGenerator);
		conservedMask.createConservedMask();
	};
};


}; // end namespace

#endif /* GENOMETASKS_TCREATEBEDMASK_H_ */
