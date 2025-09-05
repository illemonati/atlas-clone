/*
 * TPileup.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TPILEUP_H_
#define GENOMETASKS_TPILEUP_H_

#include "TSequencedData.h"
#include "TSiteTraverser.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/TOutputFile.h"

#include "TBamWindowTraverser.h"
#include "TAllelicDepthCounts.h"
#include "coretools/enum.h"

namespace GenomeTasks {

//---------------------------------
// TPileup
//---------------------------------
class TPileup {
	using Rho         = coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base, 5>, genometools::Base, 5>;
	using Transitions = coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<std::vector<Rho>, BAM::End>, BAM::Strand>, BAM::Mate>;
	using PrevBases   = coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base, 5>, genometools::Base, 5>, genometools::Base, 5>, BAM::Strand>, BAM::Mate>;

	enum class Print: size_t {min, Depth=min, Bases, SampleBases, Qualities, Alleles, Mates, Strands, Likelihoods, HML, PrintAll, max = PrintAll};
	enum class Hist : size_t {min, Depth=min, Qualities, Contexts, AllelicDepth, Transitions, PrevBases, max };

	BAM::TSiteTraverser _siteTraverser;

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

	Transitions _transitionsChr{};
	Transitions _transitionsTot{};

	PrevBases _prevBases{};

	coretools::TStrongBitSet<Print, coretools::index(Print::max) + 1> _printSettings{0};
	coretools::TStrongBitSet<Hist> _histSettings{0};

	TAllelicDepthCounts _counts;
	bool _writeEmpty;
	bool _onlyHistograms;

	void _traverseSites();
	void _addHist(const GenotypeLikelihoods::TSite& Site);
	void _printSite(const GenotypeLikelihoods::TSite& Site);
	void _endChromosome(const genometools::TChromosome& Chr);

public:
	TPileup();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TPILEUP_H_ */
