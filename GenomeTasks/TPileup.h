/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include "TSequencedData.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/TOutputFile.h"

#include "TBamWindowTraverser.h"
#include "TAllelicDepthCounts.h"

namespace GenomeTasks {

//---------------------------------
// TPileup
//---------------------------------
class TPileup final : public TBamWindowTraverser<WindowType::MultiBam> {
	using Rho         = coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base, 5>, genometools::Base, 5>;
	using Transitions = coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<std::vector<Rho>, BAM::End>, BAM::Strand>, BAM::Mate>;
	using PrevBases   = coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base, 5>, genometools::Base, 5>, genometools::Base, 5>, BAM::Strand>, BAM::Mate>;

	enum class Print: size_t {min, OnlySitesWithData=min, Depth, Bases, SampleBases, Qualities, Alleles, Mates, Strand, Likelihoods, HML, max};
	enum class Hist: size_t {min, Depths, Quality, Contexts, AllelicDepth, Transitions, PrevBases, StrandMate, max};

	coretools::TOutputFile _out;
	coretools::TOutputFile _outDepthHistogram;
	coretools::TOutputFile _outDepthPerChromosome;
	coretools::TOutputFile _outTransitions;
	coretools::TOutputFile _outTransitionsRel;
	coretools::TOutputFile _outTransitionsPsi;
	coretools::TOutputFile _outTransitionsRho;

	coretools::TCountDistribution<> _depthPerSite;
	coretools::TCountDistribution<> _depthPerSitePerChromosome;
	coretools::TCountDistributionVector<> _qualDist;
	coretools::TCountDistributionVector<> _contextDist;

	Transitions _transitionsChr;
	Transitions _transitionsTot;

	PrevBases _prevBases;

	std::vector<coretools::TStrongArray<coretools::TStrongArray<size_t, BAM::Strand>, BAM::Mate>> _strandMate;

	coretools::TStrongBitSet<Print> _printSettings;
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
