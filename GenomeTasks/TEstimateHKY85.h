/*
 * TEstimateGenotypeDistribution.h
 */

#ifndef GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_
#define GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_

#include "PMD/TModel.h"
#include "coretools/Containers/TMultiVector.h"
#include "coretools/Files/TOutputFile.h"

#include "TGenotypeDistribution.h"
#include "TBamWindowTraverser.h"
#include <memory>

namespace GenomeTasks {
class TEstimateHKY85 final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	enum class Sample : size_t {min = 0, reads=min, sites, upToDepth, max};

	std::unique_ptr<GenotypeLikelihoods::TGenotypeDistribution> _genoDist;
	std::unique_ptr<GenotypeLikelihoods::PMD::TWithPMD> _pmd;

	// EM
	size_t _numEMIterations;
	double _minDeltaLL;
	size_t _totMaskedSites = 0;
	size_t _totSites       = 0;

	coretools::TOutputFile _out;
	std::vector<double> _depthOrProbs;
	bool _genomeWide = false;
	Sample _sample;

	// genomeWide data
	struct TStats {
		double NData    = 0;
		size_t NMissing = 0;
	};
	TStats _stats_full;
	std::vector<GenotypeLikelihoods::TSite> _sites_full;
	coretools::TMultiVector<std::vector<GenotypeLikelihoods::TSite>> _sites_P;
	coretools::TMultiVector<TStats> _stats_P;


	void _handleGenomeWide(GenotypeLikelihoods::TWindow& window);
	void _handlePerWindow(GenotypeLikelihoods::TWindow& window);

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}

	double _runEM(const std::vector<GenotypeLikelihoods::TSite>& Sites);
	double _LL(const std::vector<GenotypeLikelihoods::TSite>& Sites);

	void _openFile();

public:
	TEstimateHKY85();
	void run();
};

} // namespace GenomeTasks

#endif
