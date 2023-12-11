/*
 * TCreateBedMask.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TCREATEBEDMASK_H_
#define GENOMETASKS_TCREATEBEDMASK_H_

#include <string>

#include "genometools/BED/TBed.h"

#include "TBamWindowTraverser.h"
#include "TGenotypeData.h"

namespace GenomeTasks{

//--------------------------------------
// TCreateBedMask
//--------------------------------------
class TCreateBedMask : public TBamWindowTraverser<WindowType::SingleBam> {
protected:
	genometools::TBed _bed;
	uint32_t _minDepth;

	void _createMask(const std::string fileTag);
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}
public:
	TCreateBedMask();
};

//--------------------------------------
// TCreateDepthBedMask
//--------------------------------------
class TCreateDepthBedMask:public TCreateBedMask{
private:
	uint32_t _maxDepth;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
public:
	TCreateDepthBedMask();
	void createDepthMask();
};

//--------------------------------------
// TCreateConservedBedMask
//--------------------------------------
class TCreateInvariantBedMask:public TCreateBedMask{
private:
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
public:
	TCreateInvariantBedMask();
	void createInvariantMask();
};

//--------------------------------------
// TCreateVariantBedMask
//--------------------------------------
class TCreateVariantBedMask:public TCreateBedMask{
private:
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
public:
	TCreateVariantBedMask();
	void createVariantMask();
};

//--------------------------------------
// TCreateNonRefBedMask
//--------------------------------------
class TCreateNonRefBedMask:public TCreateBedMask{
private:
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
public:
	TCreateNonRefBedMask();
	void createVariantMask();
};

struct TMaskCreator {
	void run() {
		// which mask?
		const std::string mask = coretools::instances::parameters().get<std::string>("type");
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
			UERROR("Unknown mask '", mask, "'! Valid types are 'depth', 'invariant', 'variant' and 'nonRef'.");
		}
	};
};


}; // end namespace

#endif /* GENOMETASKS_TCREATEBEDMASK_H_ */
