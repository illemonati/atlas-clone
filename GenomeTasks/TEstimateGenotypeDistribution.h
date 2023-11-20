/*
 * TEstimateGenotypeDistribution.h
 *
 */

#ifndef GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_
#define GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_

#include "TGenome_OLD.h"
#include "coretools/Files/TOutputFile.h"

namespace GenomeTasks {
class TEstimateGenotypeDistribution final : public old::TGenome_windows {
private:
	std::unique_ptr<GenotypeLikelihoods::TGenotypeDistribution> _genoDist;
	std::vector<GenotypeLikelihoods::TSite> _sites;
	size_t _numEMIterations;
	double _minDeltaLL;
	coretools::TOutputFile _out;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _handleAlignment(BAM::TAlignment&) override {}

	double _runEM();
	double _LL();

public:
	TEstimateGenotypeDistribution();
	void run();
};

} // namespace GenomeTasks

#endif
