/*
 * TEstimateGenotypeDistribution.h
 */

#ifndef GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_
#define GENOMETASKS_TESTIMATEGENOTYPEDISTRIBUTION_H_

#include "coretools/Containers/TNestedVector.h"
#include "coretools/Containers/TStandarray.h"
#include "coretools/Files/TOutputFile.h"

#include "TGenotypeDistribution.h"
#include "TBamWindowTraverser.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/Containers.h"
#include <memory>

namespace GenomeTasks {

class TEstimateHKY85 final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	enum class Sample : size_t {min = 0, readsProb=min, readsDepth, sitesProb, sitesDepth, upToDepth, max};

	std::unique_ptr<GenotypeLikelihoods::TGenotypeDistribution> _genoDist;

	// EM
	size_t _numEMIterations;
	double _minDeltaLL;
	size_t _totMaskedSites = 0;
	size_t _totSites       = 0;
	size_t _nRounds        = 1;

	coretools::TOutputFile _out;
	std::vector<double> _depthOrProbs;
	bool _genomeWide = false;
	Sample _sample;

	std::vector<genometools::Base> _refBases;
	std::vector<genometools::TGenotypeLikelihoods> _P_g_I_ds;

	// genomeWide data
	struct TStats {
		double NData    = 0;
		size_t NMissing = 0;
	};
	TStats _stats;
	coretools::TNestedVector<size_t> _readIDs;
	using Standarray = coretools::TStrongStandarray<float, genometools::Base>;
	coretools::TNestedVector<Standarray> _P_d_I_b;
	size_t _lastReadID = 0;

	bool _downSample() const noexcept {return !_depthOrProbs.empty();}
	bool _downsampeDepth() {return _sample == Sample::readsDepth || _sample == Sample::sitesDepth;}

	void _handleGenomeWide(GenotypeLikelihoods::TWindow& window);
	void _handlePerWindow(GenotypeLikelihoods::TWindow& window);
	void _handleGenomeWide();

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}

	void _addSites(const GenotypeLikelihoods::TWindow &Window);
	template<typename Container>
	double _addToPg(genometools::TGenotypeLikelihoods &P_g_I_di, const Container &P_dij_I_b, double sum) {
		using genometools::Genotype;
		using genometools::Base;
		using coretools::P;
		const auto isInvariant = _genoDist->ploidy() == genometools::Ploidy::haploid;
		const double sum_inv   = 1. / sum;
		sum                    = 0.;
		for (auto k = Base::min; k < Base::max; ++k) {
			const auto kk = genometools::genotype(k, k);
			P_g_I_di[kk] *= P(P_dij_I_b[k] * sum_inv);
			sum += P_g_I_di[kk];
		}
		if (!isInvariant) {
			for (const auto kl : {Genotype::AC, Genotype::AG, Genotype::AT, Genotype::CG, Genotype::CT, Genotype::GT}) {
				const auto k = genometools::first(kl);
				const auto l = genometools::second(kl);
				P_g_I_di[kl] *= P(0.5 * (P_dij_I_b[k] + P_dij_I_b[l]) * sum_inv);
				sum += P_g_I_di[kl];
			}
		}

		return sum;
	}
	std::pair<size_t, size_t> _downsampeSites(double ProbOrDepth, double depthInit);
	double _runEM();
	double _LL();

	void _openFile();

public:
	TEstimateHKY85();
	void run();
};

} // namespace GenomeTasks

#endif
