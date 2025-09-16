/*
 * TEstimateMutationLoad.h
 */

#ifndef GENOMETASKS_TESTIMATEMUTATIONLOAD_H_
#define GENOMETASKS_TESTIMATEMUTATIONLOAD_H_

#include <vector>

#include "TSiteTraverser.h"
#include "coretools/Containers/TView.h"
#include "genometools/Genotypes/Base.h"
#include "stattools/EM/TEMPriorIndependent.h"
#include "stattools/EM/TLatentVariable.h"

#include "TBamWindowTraverser.h"
#include "genometools/Genotypes/Containers.h"

namespace GenomeTasks{
namespace MutationLoad {

using PrecisionType = double;
using NumStatesType = int;
using LengthType    = size_t;

//------------------------------------------------
// TGenotypeProbabilities
//------------------------------------------------
class TGenotypeProbabilities{
private:
	std::array<PrecisionType, 4> _pi;
	coretools::TStrongArray<coretools::TStrongArray<PrecisionType, genometools::Genotype>, genometools::Base> _genotypeProbs;

	void _calculateGenotypeProbs();
public:
	TGenotypeProbabilities();

	void setPi(const std::array<PrecisionType, 4>& Pi);
	const std::array<PrecisionType, 4>& getPi() const { return _pi; };
	PrecisionType operator()(genometools::Base PreferredBase, genometools::Genotype Geno) const {
		return _genotypeProbs[PreferredBase][Geno];
	}
};

//------------------------------------------------
// TMutationLoadEMPrior
//------------------------------------------------
class TMutationLoadEMPrior : public stattools::TEMPriorIndependent_base<PrecisionType, NumStatesType, LengthType>{
private:
	using TPiIndex = coretools::TStrongArray<coretools::TStrongArray<NumStatesType, genometools::Genotype>, genometools::Base>;
	coretools::TConstView<genometools::TGenotypeLikelihoods> _siteLikelihoods;
	coretools::TConstView<genometools::Base> _refBases;

	TGenotypeProbabilities _genoProbs;
	TPiIndex _piIndex{initIndex()};
	std::array<double, 4> _tmpPiForEstimation;

	static TPiIndex initIndex();

public:
	TMutationLoadEMPrior(coretools::TConstView<genometools::TGenotypeLikelihoods> SiteLikelihoods, coretools::TConstView<genometools::Base> RefBases) : 
		TEMPriorIndependent_base(genometools::TGenotypeLikelihoods::capacity),
		_siteLikelihoods(SiteLikelihoods), _refBases(RefBases) {};

	PrecisionType operator()(LengthType Index, NumStatesType State) const override;

	// EM functions
	void prepareEMParameterEstimationOneIteration() override;
	void handleEMParameterEstimationOneIteration(size_t Index, const stattools::TDataVector<PrecisionType, NumStatesType> &Weights) override;
	void finalizeEMParameterEstimationOneIteration() override;
	const std::array<double, 4>& getPi() const { return _genoProbs.getPi(); }
	void reportEMParameters() override;
};

//------------------------------------------------
// TMutationLoadLatentVariable
//------------------------------------------------
class TMutationLoadLatentVariable final : public stattools::TLatentVariable<PrecisionType, NumStatesType, LengthType>{
private:
	coretools::TConstView<genometools::TGenotypeLikelihoods> _siteLikelihoods;
	coretools::TConstView<genometools::Base> _refBases;

public:
	TMutationLoadLatentVariable(coretools::TConstView<genometools::TGenotypeLikelihoods> SiteLikelihoods,
								coretools::TConstView<genometools::Base> RefBases)
		: _siteLikelihoods(SiteLikelihoods), _refBases(RefBases){};

	// EM functions
	void calculateEmissionProbabilities(size_t Index, stattools::TDataVector<PrecisionType, NumStatesType> &Emission) const override;
};

} // namespace MutationLoad
//-----------------------------------
// TEstimateMutationLoad
//-----------------------------------
class TEstimateMutationLoad  {
private:
	BAM::TSiteTraverser _siteTraverser{genometools::Morphic::Mono};
	std::vector<genometools::TGenotypeLikelihoods> _siteLikelihoods;
	std::vector<genometools::Base> _refBases;
	std::string _fileName;

	void _traverseSites();
	void _addSite(const GenotypeLikelihoods::TSite& site);
public:
	TEstimateMutationLoad();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TESTIMATEMUTATIONLOAD_H_ */
