/*
 * TEstimateMutationLoad.h
 */

#ifndef GENOMETASKS_TESTIMATETHETA_H_
#define GENOMETASKS_TESTIMATETHETA_H_

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
#include "TGenotypeData.h"
#include "genometools/GenotypeTypes.h"

namespace GenomeTasks {

namespace MutationLoad {

//------------------------------------------------
// TSiteData
//------------------------------------------------
class TSiteData {
public:
	GenotypeLikelihoods::TGenotypeLikelihoods likelihoods;
	genometools::Base preferredBase;

	TSiteData(const GenotypeLikelihoods::TGenotypeLikelihoods & Likelihoods, const genometools::Base PreferredBase): likelihoods(Likelihoods), preferredBase(PreferredBase){};
};

//------------------------------------------------
// TGenotypeProbabilities
//------------------------------------------------
class TGenotypeProbabilities{
private:
	std::array<double, 4> _pi;
	coretools::TStrongArray<coretools::TStrongArray<double, genometools::Genotype>, genometools::Base> _genotypeProbs;

	void _calculateGenotypeProbs();
public:
	TGenotypeProbabilities();
	~TGenotypeProbabilities() = default;
	void setPi(std::array<double, 4> Pi);
	const std::array<double, 4>& getPi() const { return _pi; };
	double operator()(genometools::Base PreferredBase, genometools::Genotype Geno){
		return _genotypeProbs[PreferredBase][Geno];
	}
};

//-------------------------------------
// TPiIndex
//-------------------------------------
class TPiIndex{
private:
	coretools::TStrongArray<coretools::TStrongArray<int, genometools::Genotype>, genometools::Base> _index;

public:
	TPiIndex();
	~TPiIndex() = default;

	int operator()(genometools::Base PreferredBase, genometools::Genotype Geno){
		return _index[PreferredBase][Geno];
	}
};


coretools::TStrongArray<coretools::TStrongArray<int, genometools::Genotype>, genometools::Base> _piIndex; //maps genotype to pi

//------------------------------------------------
// TMutationLoadLatentVariable
//------------------------------------------------


class TMutationLoadLatentVariable final : public stattools::TLatentVariable<double, uint8_t, size_t>{
private:
	std::vector<TSiteData> _sites;
	TGenotypeProbabilities _genoProbs;
	TPiIndex _piIndex;
	std::array<double, 4> _tmpPiForEstimation;

public:
	void addData(const GenotypeLikelihoods::TGenotypeLikelihoods & Likelihoods, const genometools::Base PreferredBase){
		_sites.emplace_back(Likelihoods, PreferredBase);
	}

	// EM functions
	void calculateEmissionProbabilities(size_t Index, stattools::TDataVector<double, uint8_t> &Emission) const override;
	void prepareEMParameterEstimationOneIteration() override;
	void handleEMParameterEstimationOneIteration(size_t Index, const stattools::TDataVector<double, uint8_t> &Weights) override;
	void finalizeEMParameterEstimationOneIteration() override;
	void reportEMParameters() override;


};

}
//-----------------------------------
// TEstimateMutationLoad
//-----------------------------------
class TEstimateMutationLoad : public TGenome_windows {
private:
	MutationLoad::TMutationLoadLatentVariable _latentVar;

	void _handleWindow() override;
	void _handleAlignment() override {}
public:
	TEstimateMutationLoad();
	void run();
};

}; // namespace GenomeTasks

#endif /* GENOMETASKS_TESTIMATETHETA_H_ */
