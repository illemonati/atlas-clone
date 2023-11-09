/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include "coretools/Containers/TBitSet.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Math/counters.h"

#include "TGenome.h"
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

	enum class Print: size_t {min, OnlySitesWithData=min, Depth, Bases, Qualities, Alleles, Mates, Strand, Likelihoods, max};
	coretools::TStrongBitSet<Print> _printSettings;

	enum class Hist: size_t {min, Depths, Quality, Contexts, AllelicDepth, max};
	coretools::TStrongBitSet<Hist> _histSettings;

	TAllelicDepthCounts _counts;
	bool _writeEmpty;
	bool _onlySummary;

	void _handleWindow() override;
	void _handleAlignment() override {}
public:
	TPileup();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPILEUP_H_ */
