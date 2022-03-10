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
	genometools::TBed _bed;
	uint32_t _minDepthForMask;

	void _createMask(const std::string fileTag);

public:
	TCreateBedMask(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
};

//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
class TCreateDepthBedMask:public TCreateBedMask{
private:
	uint32_t _maxDepthForMask;

	void _handleWindow();
public:
	TCreateDepthBedMask(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
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
	TCreateInvariantBedMask(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
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
	TCreateVariantBedMask(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
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
	TCreateNonRefBedMask(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void createVariantMask();
};


//--------------------------------------
// Tasks
//--------------------------------------
class TTask_createMask:public coretools::TTask{
public:
	TTask_createMask(){ _explanation = "Creating a mask BED file"; };

	void run(){
		using namespace coretools::instances;

		//which mask?
		std::string mask = parameters().getParameter<std::string>("mask");
		if(mask == "depth"){
			TCreateDepthBedMask depthMask(parameters(), &logfile(), &randomGenerator());
			depthMask.createDepthMask();
		} else if(mask == "nonRef"){

		} else if(mask == "invariant"){
			TCreateInvariantBedMask conservedMask(parameters(), &logfile(), &randomGenerator());
			conservedMask.createInvariantMask();
		} else if(mask == "variant"){

		} else {
			throw "Unknown mask '" + mask + "'! Valid types are 'depth', 'invariant', 'nonRef' and 'polymorphic'";
		}
	};
};


}; // end namespace

#endif /* GENOMETASKS_TCREATEBEDMASK_H_ */
