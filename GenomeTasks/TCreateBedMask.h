/*
 * TCreateBedMask.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCREATEBEDMASK_H_
#define GENOMETASKS_TCREATEBEDMASK_H_

#include <stdint.h>
#include <exception>
#include <string>

#include "genometools/BED/TBed.h"
#include "TGenome.h"
#include "TGenotypeData.h"
#include "coretools/Main/TTask.h"

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
	TCreateBedMask();
};

//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
class TCreateDepthBedMask:public TCreateBedMask{
private:
	uint32_t _maxDepthForMask;

	void _handleWindow() override;
	void _handleAlignment() override {}
public:
	TCreateDepthBedMask();
	void createDepthMask();
};

//--------------------------------------
// TCreateConservedBedMask
//--------------------------------------
class TCreateInvariantBedMask:public TCreateBedMask{
private:
	void _handleWindow() override;
	void _handleAlignment() override {}

	//tmp variables
	GenotypeLikelihoods::TBaseCounts _baseCounts;
public:
	TCreateInvariantBedMask();
	void createInvariantMask();
};

//--------------------------------------
// TCreateVariantBedMask
//--------------------------------------
class TCreateVariantBedMask:public TCreateBedMask{
private:
	void _handleWindow() override;

	//tmp variables
	GenotypeLikelihoods::TBaseCounts _baseCounts;
public:
	TCreateVariantBedMask();
	void createVariantMask();
};

//--------------------------------------
// TCreateNonRefBedMask
//--------------------------------------
class TCreateNonRefBedMask:public TCreateBedMask{
private:
	void _handleWindow();
public:
	TCreateNonRefBedMask();
	void createVariantMask();
};

struct TMaskCreator {
	void run() {
		// which mask?
		std::string mask = coretools::instances::parameters().getParameter<std::string>("mask");
		if (mask == "depth") {
			TCreateDepthBedMask depthMask;
			depthMask.createDepthMask();
		} else if (mask == "nonRef") {

		} else if (mask == "invariant") {
			TCreateInvariantBedMask conservedMask;
			conservedMask.createInvariantMask();
		} else if (mask == "variant") {

		} else {
			UERROR("Unknown mask '", mask, "'! Valid types are 'depth', 'invariant', 'nonRef' and 'polymorphic'");
		}
	};
};


}; // end namespace

#endif /* GENOMETASKS_TCREATEBEDMASK_H_ */
