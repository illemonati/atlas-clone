/*
 * TEstimateMutationLoad.cpp
*/

#include "TEstimateMutationLoad.h"
#include "coretools/traits.h"
#include "stattools/EM/TEM.h"

using coretools::instances::logfile;

namespace GenomeTasks {

namespace MutationLoad {


//------------------------------------------------
// TGenotypeProbabilities
//------------------------------------------------
TGenotypeProbabilities::TGenotypeProbabilities(){
	_pi = {0.25, 0.25, 0.25, 0.25};
	_calculateGenotypeProbs();
}

void TGenotypeProbabilities::setPi(std::array<double, 4> Pi){
	_pi = Pi;
	_calculateGenotypeProbs();
}

void TGenotypeProbabilities::_calculateGenotypeProbs(){
	// calculate the genotype probabilities according to the mutation load model
	// calculate it for each possible preferred base
	using genometools::Base;
	using genometools::genotype;

	for(Base r = Base::min; r < Base::max; ++r){
		// homozygous preferred
		_genotypeProbs[r][genotype(r,r)] = _pi[0];

		for(Base a = Base::min; a < Base::max; ++a){
			if(a != r){
				// heterozygous preferred
				_genotypeProbs[r][genotype(r,a)] = _pi[1] / 3.0;

				// homozygous alternative
				_genotypeProbs[r][genotype(a,a)] = _pi[2] / 3.0;

				for(Base b = Base::min; b < Base::max; ++b){
					if(b != a && b != r){
						_genotypeProbs[r][genotype(a,b)] = _pi[3] / 3.0;
					}
				}
			}
		}
	}
}

//-------------------------------------
// TPiIndex
//-------------------------------------
TPiIndex::TPiIndex(){
	// maps genotype to pi
	using genometools::Base;
	using genometools::genotype;

	for(Base r = Base::min; r < Base::max; ++r){
		// homozygous preferred
		_index[r][genotype(r,r)] = 0;

		for(Base a = Base::min; a < Base::max; ++a){
			if(a != r){
				// heterozygous preferred
				_index[r][genotype(r,a)] = 1;

				// homozygous alternative
				_index[r][genotype(a,a)] = 2;

				for(Base b = Base::min; b < Base::max; ++b){
					if(b != a && b != r){
						_index[r][genotype(a,b)] = 3;
					}
				}
			}
		}
	}
}

//------------------------------------------------
// TMutationLoadEMPrior
//------------------------------------------------

PrecisionType TMutationLoadEMPrior::operator()(size_t Index, NumStatesType State) const{
	return _genoProbs(_sites[Index].preferredBase, genometools::Genotype(State));
}

void TMutationLoadEMPrior::prepareEMParameterEstimationOneIteration(){
	_tmpPiForEstimation = {0.0, 0.0, 0.0, 0.0};
}

void TMutationLoadEMPrior::handleEMParameterEstimationOneIteration(size_t Index, const stattools::TDataVector<PrecisionType, NumStatesType> &Weights){
	for(genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
		_tmpPiForEstimation[ _piIndex(_sites[Index].preferredBase, g) ] += Weights[coretools::underlying(g)];
	}
}

void TMutationLoadEMPrior::finalizeEMParameterEstimationOneIteration(){
	coretools::normalize(_tmpPiForEstimation);
	_genoProbs.setPi(_tmpPiForEstimation);
}

void TMutationLoadEMPrior::reportEMParameters(){
	logfile().list("Pi = ", _genoProbs.getPi());
}


//------------------------------------------------
// TMutationLoadLatentVariable
//------------------------------------------------
void TMutationLoadLatentVariable::calculateEmissionProbabilities(size_t Index, stattools::TDataVector<PrecisionType, NumStatesType> &Emission) const {
	Emission = _sites[Index].likelihoods;
}

} // end namespace MutationLoad

//------------------------------------------------
// TEstimateMutationLoad
//------------------------------------------------
void TEstimateMutationLoad::_handleWindow() {
	// adding sites to estimator
	logfile().listFlushTime("Calculating genotype likelihoods and storing data ...");
	try {
		auto thesePositions = _subset->getPositionInWindow(_window);
		for(auto& it : thesePositions){
			uint32_t internalPos = it - _window.from();
			GenotypeLikelihoods::TSite& site = _window[internalPos];
			GenotypeLikelihoods::TGenotypeLikelihoods genoLik = _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(site);
			_sites.emplace_back(genoLik, it.ref());
		}
	} catch (...) {
		UERROR("Failed to allocate sufficient memory to store the data for so many sites. Consider using fewer sites.");
	}
	logfile().doneTime();
};


void TEstimateMutationLoad::run(){
	// traverse BAM file and store data
	_traverseBAMWindows();

	// now run estimation
	MutationLoad::TMutationLoadEMPrior prior(_sites);
	MutationLoad::TMutationLoadLatentVariable latentVar(_sites);

	stattools::TEM<MutationLoad::PrecisionType, MutationLoad::NumStatesType, MutationLoad::LengthType> EM(prior, latentVar);

}

} // end namespace GenomeTasks
