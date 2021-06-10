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

	void _createMask(const std::string fileTag);

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
class TCreateInvariantBedMask:public TCreateBedMask{
private:
	void _handleWindow();

	//tmp variables
	GenotypeLikelihoods::TBaseCounts _baseCounts;
public:
	TCreateInvariantBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void createInvariantMask();
};

//--------------------------------------
// TCreateVariantBedMask
//--------------------------------------
class TCreateVariantBedMask:public TCreateBedMask{
private:
	void _handleWindow();

	//tmp variables
	GenotypeLikelihoods::TBaseCounts _baseCounts;
public:
	TCreateVariantBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void createVariantMask();
};

//--------------------------------------
// TCreateNonRefBedMask
//--------------------------------------
class TCreateNonRefBedMask:public TCreateBedMask{
private:
	void _handleWindow();

	//tmp variables
	GenotypeLikelihoods::TBaseCounts _baseCounts;
public:
	TCreateNonRefBedMask(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void createVariantMask();
};


//--------------------------------------
// Tasks
//--------------------------------------
class TTask_createMask:public coretools::TTask{
public:
	TTask_createMask(){ _explanation = "Creating a mask BED file"; };

	void run(coretools::TParameters & Parameters, coretools::TLog* Logfile){

		//which mask?
		std::string mask = Parameters.getParameter<std::string>("mask");
		if(mask == "depth"){
			TCreateDepthBedMask depthMask(Parameters, Logfile, _randomGenerator);
			depthMask.createDepthMask();
		} else if(mask == "nonRef"){

		} else if(mask == "invariant"){
			TCreateInvariantBedMask conservedMask(Parameters, Logfile, _randomGenerator);
			conservedMask.createInvariantMask();
		} else if(mask == "variant"){

		} else {
			throw "Unknown mask '" + mask + "'! Valid types are 'depth', 'invariant', 'nonRef' and 'polymorphic'";
		}
	};
};


}; // end namespace

#endif /* GENOMETASKS_TCREATEBEDMASK_H_ */
