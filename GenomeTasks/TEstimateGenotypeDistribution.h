/*
 * TEstimateGenotypeDistribution.h
 */

#ifndef GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_
#define GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_

#include "coretools/Files/TOutputFile.h"

#include "TGenotypeDistribution.h"
#include "TBamWindowTraverser.h"

namespace GenomeTasks {
class TEstimateGenotypeDistribution final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	std::unique_ptr<GenotypeLikelihoods::TGenotypeDistribution> _genoDist;
	std::vector<coretools::Probability> _probs;

	// genomeWide
	std::vector<std::vector<GenotypeLikelihoods::TSite>> _sites_P;
	struct TStats {
		double NData    = 0;
		size_t NMissing = 0;
	};
	std::vector<TStats> _stats;

	size_t _numEMIterations;
	double _minDeltaLL;
	bool _genomeWide       = false;
	size_t _totMaskedSites = 0;
	size_t _totSites       = 0;
	coretools::TOutputFile _out;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}

	double _runEM(const std::vector<GenotypeLikelihoods::TSite>& Sites);
	double _LL(const std::vector<GenotypeLikelihoods::TSite>& Sites);

public:
	TEstimateGenotypeDistribution();
	void run();
};

} // namespace GenomeTasks

#endif
