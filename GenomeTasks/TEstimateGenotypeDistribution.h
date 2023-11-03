/*
 * TEstimateGenotypeDistribution.h
 *
 */

#ifndef GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_
#define GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_

#include "TGenome.h"

namespace GenomeTasks {
class TEstimateGenotypeDistribution final : public TGenome_windows {
private:
	std::unique_ptr<GenotypeLikelihoods::TGenotypeDistribution> _genoDist;
	std::vector<GenotypeLikelihoods::TSite> _sites;
	size_t _numEMIterations;
	double _minDeltaLL;

	void _handleWindow() override;
	void _handleAlignment() override {}

	void _runEM();
	double _LL();

public:
	TEstimateGenotypeDistribution();
	void run();
};

} // namespace GenomeTasks

#endif
