/*
 * TEstimateMutationLoad.h
 */

#ifndef GENOMETASKS_TESTIMATEMUTATIONLOAD_H_
#define GENOMETASKS_TESTIMATEMUTATIONLOAD_H_

#include <stdint.h>
#include <vector>

#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "genometools/BED/TBed.h"
#include "TGenome.h"
#include "coretools/Main/TTask.h"
#include "TWindow.h"
#include "coretools/Types/probability.h"
#include "stattools/EM/TLatentVariable.h"
#include "stattools/EM/TEMPriorIndependent.h"
#include "TGenotypeData.h"
#include "genometools/GenotypeTypes.h"

namespace GenomeTasks {

namespace MutationLoad {

using genometools::Genotype;
using genometools::Base;
using GenotypeLikelihoods::TGenotypeLikelihoods;

using PrecisionType = double;
using NumStatesType = int;
using LengthType = size_t;


//------------------------------------------------
// TSiteData
//------------------------------------------------
class TSiteData {
public:
	TGenotypeLikelihoods likelihoods;
	Base preferredBase;

	TSiteData(const TGenotypeLikelihoods & Likelihoods, const Base PreferredBase):
	preferredBase(PreferredBase),
	likelihoods(Likelihoods) {};

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
	coretools::TStrongArray<coretools::TStrongArray<PrecisionType, Genotype>, Base> _genotypeProbs;

	void _calculateGenotypeProbs();
public:
	TGenotypeProbabilities();
	~TGenotypeProbabilities() = default;
	void setPi(std::array<double, 4> Pi);
	const std::array<double, 4>& getPi() const { return _pi; };
	double operator()(Base PreferredBase, Genotype Geno) const {
		return _genotypeProbs[PreferredBase][Geno];
	}
};

//-------------------------------------
// TPiIndex
//-------------------------------------
class TPiIndex{
private:
	coretools::TStrongArray<coretools::TStrongArray<NumStatesType, Genotype>, Base> _index;

public:
	TPiIndex();
	~TPiIndex() = default;

	int operator()(Base PreferredBase, Genotype Geno){
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
		TEMPriorIndependent_base(TGenotypeLikelihoods::capacity),
		_sites(Sites) {};

	PrecisionType operator()(LengthType Index, NumStatesType State) const override;

	// EM functions
	void prepareEMParameterEstimationOneIteration() override;
	void handleEMParameterEstimationOneIteration(size_t Index, const stattools::TDataVector<PrecisionType, NumStatesType> &Weights) override;
	void finalizeEMParameterEstimationOneIteration() override;
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
class TEstimateMutationLoad : public TGenome_windows {
private:
	std::vector<MutationLoad::TSiteData> _sites;

	void _handleWindow() override;
	void _handleAlignment() override {}
public:
	TEstimateMutationLoad();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TESTIMATEMUTATIONLOAD_H_ */
