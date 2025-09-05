/*
 * TCreateBedMask.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCREATEBEDMASK_H_
#define GENOMETASKS_TCREATEBEDMASK_H_

#include <string>

#include "TSiteTraverser.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "genometools/TBed.h"

#include "TBamWindowTraverser.h"

namespace GenomeTasks{

//--------------------------------------
// TCreateBedMask
//--------------------------------------
class TCreateBedMask {
protected:
	BAM::TSiteTraverser _siteTraverser;

	genometools::TBed _bed;
	uint32_t _minDepth;

	void _createMask(std::string_view fileTag);
	virtual void _handleSite(const GenotypeLikelihoods::TSite& Site, genometools::TGenomePosition Position) = 0;
public:
	TCreateBedMask();
};

//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
class TCreateDepthBedMask:public TCreateBedMask{
private:
	uint32_t _maxDepth;

	void _handleSite(const GenotypeLikelihoods::TSite& Site, genometools::TGenomePosition Position) override;
public:
	TCreateDepthBedMask();
	void createDepthMask();
};

//--------------------------------------
// TCreateConservedBedMask
//--------------------------------------
class TCreateInvariantBedMask:public TCreateBedMask{
	void _handleSite(const GenotypeLikelihoods::TSite& Site, genometools::TGenomePosition Position) override;
public:
	TCreateInvariantBedMask();
	void createInvariantMask();
};

//--------------------------------------
// TCreateVariantBedMask
//--------------------------------------
class TCreateVariantBedMask:public TCreateBedMask{
	void _handleSite(const GenotypeLikelihoods::TSite& Site, genometools::TGenomePosition Position) override;
public:
	TCreateVariantBedMask();
	void createVariantMask();
};

//--------------------------------------
// TCreateNonRefBedMask
//--------------------------------------
class TCreateNonRefBedMask:public TCreateBedMask{
	void _handleSite(const GenotypeLikelihoods::TSite& Site, genometools::TGenomePosition Position) override;
public:
	TCreateNonRefBedMask();
	void createVariantMask();
};

struct TMaskCreator {
	void run() {
		// which mask?
		const auto mask = coretools::instances::parameters().get("type");
		if (mask == "depth") {
			TCreateDepthBedMask depthMask;
			depthMask.createDepthMask();
		} else if (mask == "nonRef") {
			TCreateNonRefBedMask mask;
			mask.createVariantMask();
		} else if (mask == "invariant") {
			TCreateInvariantBedMask conservedMask;
			conservedMask.createInvariantMask();
		} else if (mask == "variant") {
			TCreateVariantBedMask mask;
			mask.createVariantMask();
		} else {
			throw coretools::TUserError("Unknown mask '", mask, "'! Valid types are 'depth', 'invariant', 'variant' and 'nonRef'.");
		}
	};
};


}; // end namespace

#endif /* GENOMETASKS_TCREATEBEDMASK_H_ */
