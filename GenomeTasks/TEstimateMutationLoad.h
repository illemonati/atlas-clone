/*
 * TEstimateMutationLoad.h
 */

#ifndef GENOMETASKS_TESTIMATEMUTATIONLOAD_H_
#define GENOMETASKS_TESTIMATEMUTATIONLOAD_H_

#include <vector>

#include "genometools/BED/TBed.h"
#include "genometools/GenotypeTypes.h"
#include "stattools/EM/TEMPriorIndependent.h"
#include "stattools/EM/TLatentVariable.h"


#include "TBamWindowTraverser.h"
#include "TGenotypeData.h"
#include "TWindow.h"

namespace GenomeTasks{
namespace MutationLoad {


using PrecisionType = double;
using NumStatesType = int;
using LengthType = size_t;


//------------------------------------------------
// TSiteData
//------------------------------------------------
class TSiteData {
public:
	GenotypeLikelihoods::TGenotypeLikelihoods likelihoods;
	genometools::Base preferredBase;

	TSiteData(const GenotypeLikelihoods::TGenotypeLikelihoods &Likelihoods, const genometools::Base PreferredBase)
		: likelihoods(Likelihoods), preferredBase(PreferredBase){};

	TSiteData(const TSiteData & other) = delete;
	TSiteData(TSiteData && other){
		likelihoods = std::move(other.likelihoods);
		preferredBase = std::move(other.preferredBase);
	};
};

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
	~TGenotypeProbabilities() = default;
	void setPi(std::array<double, 4> Pi);
	const std::array<double, 4>& getPi() const { return _pi; };
	double operator()(genometools::Base PreferredBase, genometools::Genotype Geno) const {
		return _genotypeProbs[PreferredBase][Geno];
	}
};

//-------------------------------------
// TPiIndex
//-------------------------------------
class TPiIndex{
private:
	coretools::TStrongArray<coretools::TStrongArray<NumStatesType, genometools::Genotype>, genometools::Base> _index;

public:
	TPiIndex();
	~TPiIndex() = default;

	int operator()(genometools::Base PreferredBase, genometools::Genotype Geno){
		return _index[PreferredBase][Geno];
	}
};

//------------------------------------------------
// TMutationLoadEMPrior
//------------------------------------------------
class TMutationLoadEMPrior : public stattools::TEMPriorIndependent_base<PrecisionType, NumStatesType, LengthType>{
private:
	std::vector<MutationLoad::TSiteData>& _sites;
	TGenotypeProbabilities _genoProbs;
	TPiIndex _piIndex;
	std::array<double, 4> _tmpPiForEstimation;

public:
	TMutationLoadEMPrior(std::vector<MutationLoad::TSiteData>& Sites) : 
		TEMPriorIndependent_base(GenotypeLikelihoods::TGenotypeLikelihoods::capacity),
		_sites(Sites) {};

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
	std::vector<MutationLoad::TSiteData>& _sites;

public:
	TMutationLoadLatentVariable(std::vector<MutationLoad::TSiteData>& Sites) : _sites(Sites) {};

	// EM functions
	void calculateEmissionProbabilities(size_t Index, stattools::TDataVector<PrecisionType, NumStatesType> &Emission) const override;
};

} // namespace MutationLoad
//-----------------------------------
// TEstimateMutationLoad
//-----------------------------------
class TEstimateMutationLoad final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	std::vector<MutationLoad::TSiteData> _sites;
	bool _parseFromBed;
	std::string _bedFileName;
	genometools::TBed _bedFile;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}

	void _addSite(const GenotypeLikelihoods::TSite& site, const genometools::Base PreferredBase);
public:
	TEstimateMutationLoad();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TESTIMATEMUTATIONLOAD_H_ */
