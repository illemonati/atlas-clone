/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include <array>
#include <set>
#include <string>

#include "coretools/Containers/TBitSet.h"
#include "coretools/Files/TFile.h"
#include "TGenome.h"
#include "TGenotypeData.h"
#include "coretools/Main/TTask.h"
#include "TAllelicDepthCounts.h"

namespace GenomeTasks {

//---------------------------------
// TPileup
//---------------------------------
class TPileup : public TGenome_windows {
private:
	coretools::TOutputFile _out;
	coretools::TOutputFile _outDepthHistogram;
	coretools::TOutputFile _outDepthPerChromosome;

	coretools::TCountDistribution<> _depthPerSite;
	coretools::TCountDistribution<> _depthPerSitePerChromosome;
	coretools::TCountDistributionVector<> _qualDist;
	coretools::TCountDistributionVector<> _contextDist;

	// what to print?
	coretools::TBitSet<8> _printSettings;
	enum {OnlySitesWithData, Depth, Bases, Qualities, Alleles, Mates, Strand, Likelihoods};

	coretools::TBitSet<8> _histSettings;
	enum {Depths, Quality, Contexts, AllelicDepth};

	TAllelicDepthCounts _counts;
	bool _writeEmpty;

	void _handleWindow() override;
	void _handleAlignment() override {}
public:
	TPileup();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPILEUP_H_ */
