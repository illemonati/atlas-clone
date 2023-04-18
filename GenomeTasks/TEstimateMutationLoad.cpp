/*
 * TEstimateMutationLoad.cpp
*/

#include "TEstimateMutationLoad.h"
#include "coretools/traits.h"

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
	using genometools::Genotype;

	for(Base r = Base::min; r < Base::max; ++r){
		// homozygous preferred
		_genotypeProbs[r][Genotype(r,r)] = _pi[0];

		for(Base a = Base::min; a < Base::max; ++a){
			if(a != r){
				// heterozygous preferred
				_genotypeProbs[r][Genotype(r,a)] = _pi[1] / 3.0;

				// homozygous alternative
				_genotypeProbs[r][Genotype(a,a)] = _pi[2] / 3.0;

				for(Base b = Base::min; b < Base::max; ++b){
					if(b != a && b != r){
						_genotypeProbs[r][Genotype(a,b)] = _pi[3] / 3.0;
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
	using genometools::Genotype;

	for(Base r = Base::min; r < Base::max; ++r){
		// homozygous preferred
		_index[r][Genotype(r,r)] = 0;

		for(Base a = Base::min; a < Base::max; ++a){
			if(a != r){
				// heterozygous preferred
				_index[r][Genotype(r,a)] = 1;

				// homozygous alternative
				_index[r][Genotype(a,a)] = 2;

				for(Base b = Base::min; b < Base::max; ++b){
					if(b != a && b != r){
						_index[r][Genotype(a,b)] = 3;
					}
				}
			}
		}
	}
}

//------------------------------------------------
// TMutationLoadLatentVariable
//------------------------------------------------


void TMutationLoadLatentVariable::calculateEmissionProbabilities(size_t Index, stattools::TDataVector<double, uint8_t> &Emission) const{
	// calculate raw emissions
	double sum = 0.0;
	for(genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
		Emission[coretools::underlying(g)] = _sites[Index].likelihoods[g] * _genoProbs(_sites[Index].preferredBase, g);
		sum += Emission[coretools::underlying(g)];
	}

	// normalize emissions
	Emission.normalize(sum);
}

void TMutationLoadLatentVariable::prepareEMParameterEstimationOneIteration(){
	_tmpPiForEstimation = {0.0, 0.0, 0.0, 0.0};
}

void TMutationLoadLatentVariable::handleEMParameterEstimationOneIteration(size_t Index, const stattools::TDataVector<double, uint8_t> &Weights){
	for(genometools::Genotype g = genometools::Genotype::min; g < genometools::Genotype::max; ++g){
		_tmpPiForEstimation[ _piIndex(_sites[Index].preferredBase, g) ] += Weights[coretools::underlying(g)];
	}
}

void TMutationLoadLatentVariable::finalizeEMParameterEstimationOneIteration(){
	double sum = std::accumulate(_tmpPiForEstimation.begin(), _tmpPiForEstimation.end(), 0.0);
	std::transform(_tmpPiForEstimation.begin(), _tmpPiForEstimation.end(), _tmpPiForEstimation.begin(), [sum](auto x){return x / sum;});
	_genoProbs.setPi(_tmpPiForEstimation);
}

void TMutationLoadLatentVariable::reportEMParameters(){
	logfile().list("Pi = ", _genoProbs.getPi());
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
			_latentVar.addData(genoLik, it.ref());
		}
	} catch (...) {
		UERROR("Failed to allocate sufficient memory to store the data for so many sites. Consider using fewer sites.");
	}
	logfile().doneTime();
};


TEstimateMutationLoad::run(){
	// traverse BAM file and store data
	_traverseBAMWindows();

	// now run estimation

}

} // end namespace GenomeTasks
