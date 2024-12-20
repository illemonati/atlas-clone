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
	coretools::TOutputFile _outTransitions;
	coretools::TOutputFile _outTransitionsPsi;
	coretools::TOutputFile _outTransitionsRho;
	coretools::TOutputFile _outRho;

	coretools::TCountDistribution<> _depthPerSite;
	coretools::TCountDistribution<> _depthPerSitePerChromosome;
	coretools::TCountDistributionVector<> _qualDist;
	coretools::TCountDistributionVector<> _contextDist;

	using Rho         = coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base>, genometools::Base>;
	using Transitions = coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<std::vector<Rho>, BAM::End>, BAM::Strand>, BAM::Mate>;

	Transitions _transitionsChr;
	Transitions _transitionsTot;

	enum class Print: size_t {min, OnlySitesWithData=min, Depth, Bases, SampleBases, Qualities, Alleles, Mates, Strand, Likelihoods, HML, max};
	coretools::TStrongBitSet<Print> _printSettings;

	enum class Hist: size_t {min, Depths, Quality, Contexts, AllelicDepth, Transitions, max};
	coretools::TStrongBitSet<Hist> _histSettings;

	TAllelicDepthCounts _counts;
	bool _writeEmpty;
	bool _onlySummary;

	void _writeTransitions(const Transitions& transitions, std::string_view Chr);
	void _writeRho(const Rho& rho, std::string_view mate, std::string_view strand, std::string_view end, std::string_view Chr);
	void _handleWindow(GenotypeLikelihoods::TWindow& Window) override;
	void _startChromosome(const genometools::TChromosome& ) override {}
	void _endChromosome(const genometools::TChromosome& Chr) override;
public:
	TPileup();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPILEUP_H_ */
