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
	std::vector<GenotypeLikelihoods::TSite> _sites;
	size_t _numEMIterations;
	double _minDeltaLL;
	coretools::TOutputFile _out;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}

	double _runEM();
	double _LL();

public:
	TEstimateGenotypeDistribution();
	void run();
};

} // namespace GenomeTasks

#endif
