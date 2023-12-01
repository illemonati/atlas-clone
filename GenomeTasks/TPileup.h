/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include "coretools/Files/TOutputFile.h"
#include "coretools/Math/counters.h"

#include "TBamWindowTraverser.h"
#include "TAllelicDepthCounts.h"

namespace GenomeTasks {

//---------------------------------
// TPileup
//---------------------------------
class TPileup final : public TBamWindowTraverser {
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

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _onChrChange(const genometools::TChromosome&) override {}
public:
	TPileup();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPILEUP_H_ */
