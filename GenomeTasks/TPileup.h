/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include "TSequencedBase.h"
#include "coretools/Files/TOutputFile.h"

#include "TBamWindowTraverser.h"
#include "TAllelicDepthCounts.h"

namespace GenomeTasks {

//---------------------------------
// TPileup
//---------------------------------
class TPileup final : public TBamWindowTraverser<WindowType::MultiBam> {
private:
	coretools::TOutputFile _out;
	coretools::TOutputFile _outDepthHistogram;
	coretools::TOutputFile _outDepthPerChromosome;

	coretools::TCountDistribution<> _depthPerSite;
	coretools::TCountDistribution<> _depthPerSitePerChromosome;
	coretools::TCountDistributionVector<> _qualDist;
	coretools::TCountDistributionVector<> _contextDist;
	coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<std::vector<size_t>, genometools::Base, 5>, genometools::Base, 5>, BAM::End> _transitionsPerPos;

	enum class Print: size_t {min, OnlySitesWithData=min, Depth, Bases, SampleBases, Qualities, Alleles, Mates, Strand, Likelihoods, HML, max};
	coretools::TStrongBitSet<Print> _printSettings;

	enum class Hist: size_t {min, Depths, Quality, Contexts, AllelicDepth, Transitions, max};
	coretools::TStrongBitSet<Hist> _histSettings;

	TAllelicDepthCounts _counts;
	bool _writeEmpty;
	bool _onlySummary;

	void _handleWindow(GenotypeLikelihoods::TWindow& Window) override;
	void _startChromosome(const genometools::TChromosome& ) override {}
	void _endChromosome(const genometools::TChromosome& Chr) override;
public:
	TPileup();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPILEUP_H_ */
