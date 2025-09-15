/*
 * TEstimateMutationLoad.cpp
 */

#include "TEstimateMutationLoad.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TParameters.h"
#include "stattools/EM/TEM.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using genometools::Genotype;
using genometools::TGenotypeLikelihoods;


namespace MutationLoad {
using genometools::Base;

//------------------------------------------------
// TGenotypeProbabilities
//------------------------------------------------
TGenotypeProbabilities::TGenotypeProbabilities() {
	_pi = {0.25, 0.25, 0.25, 0.25};
	_calculateGenotypeProbs();
}

void TGenotypeProbabilities::setPi(std::array<double, 4> Pi) {
	_pi = Pi;
	_calculateGenotypeProbs();
}

void TGenotypeProbabilities::_calculateGenotypeProbs() {
	// calculate the genotype probabilities according to the mutation load model
	// calculate it for each possible preferred base
	using genometools::Base;
	using genometools::genotype;

	for (Base r = Base::min; r < Base::max; ++r) {
		// homozygous preferred
		_genotypeProbs[r][genotype(r, r)] = _pi[0];

		for (Base a = Base::min; a < Base::max; ++a) {
			if (a != r) {
				// heterozygous preferred
				_genotypeProbs[r][genotype(r, a)] = _pi[1] / 3.0;

				// homozygous alternative
				_genotypeProbs[r][genotype(a, a)] = _pi[2] / 3.0;

				for (Base b = Base::min; b < Base::max; ++b) {
					if (b != a && b != r) { _genotypeProbs[r][genotype(a, b)] = _pi[3] / 3.0; }
				}
			}
		}
	}
}

//-------------------------------------
// TPiIndex
//-------------------------------------
TPiIndex::TPiIndex() {
	// maps genotype to pi
	using genometools::Base;
	using genometools::genotype;

	for (Base r = Base::min; r < Base::max; ++r) {
		// homozygous preferred
		_index[r][genotype(r, r)] = 0;

		for (Base a = Base::min; a < Base::max; ++a) {
			if (a != r) {
				// heterozygous preferred
				_index[r][genotype(r, a)] = 1;

				// homozygous alternative
				_index[r][genotype(a, a)] = 2;

				for (Base b = Base::min; b < Base::max; ++b) {
					if (b != a && b != r) { _index[r][genotype(a, b)] = 3; }
				}
			}
		}
	}
}

//------------------------------------------------
// TMutationLoadEMPrior
//------------------------------------------------

PrecisionType TMutationLoadEMPrior::operator()(size_t Index, NumStatesType State) const {
	return _genoProbs(_sites[Index].preferredBase, Genotype(State));
}

void TMutationLoadEMPrior::prepareEMParameterEstimationOneIteration() { _tmpPiForEstimation = {0.0, 0.0, 0.0, 0.0}; }

void TMutationLoadEMPrior::handleEMParameterEstimationOneIteration(
    size_t Index, const stattools::TDataVector<PrecisionType, NumStatesType> &Weights) {
	for (Genotype g = Genotype::min; g < Genotype::max; ++g) {
		_tmpPiForEstimation[_piIndex(_sites[Index].preferredBase, g)] += Weights[coretools::index(g)];
	}
}

void TMutationLoadEMPrior::finalizeEMParameterEstimationOneIteration() {
	coretools::normalize(_tmpPiForEstimation);
	_genoProbs.setPi(_tmpPiForEstimation);
}

void TMutationLoadEMPrior::reportEMParameters() { logfile().list("Pi = ", _genoProbs.getPi()); }

//------------------------------------------------
// TMutationLoadLatentVariable
//------------------------------------------------
void TMutationLoadLatentVariable::calculateEmissionProbabilities(
    size_t Index, stattools::TDataVector<PrecisionType, NumStatesType> &Emission) const {
	Emission.copyToCurrent(_sites[Index].likelihoods.data());
}

} // end namespace MutationLoad

//------------------------------------------------
// TEstimateMutationLoad
//------------------------------------------------
void TEstimateMutationLoad::_addSite(const GenotypeLikelihoods::TSite &site) {
	//if (site.empty() || PreferredBase == genometools::Base::N) return;

	genometools::TGenotypeLikelihoods genoLik = _siteTraverser.errorModels().calculateGenotypeLikelihoods(site);
	_sites.emplace_back(genoLik, site.refBase);
}

void TEstimateMutationLoad::_traverseSites() {
	logfile().startIndent("Calculating genotype likelihoods and storing data ...");
	for (;!_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		for (; !_siteTraverser.endOfCurChr(); _siteTraverser.nextSite()) {
			_addSite(_siteTraverser.site());
		}
	}
	logfile().endIndent();
}

TEstimateMutationLoad::TEstimateMutationLoad()  {
	using coretools::instances::logfile;
	using coretools::instances::parameters;
	_siteTraverser.requireSingleBAM();
	_siteTraverser.skipEmpty();

	if (_siteTraverser.alleles()) {
		logfile().list("Limiting analysis to sites listed in alleles file.");
		_fileName = parameters().get("alleles");
	} else if (_siteTraverser.regions()) {
		logfile().list("Limiting analysis to sites listed in BED file.");
		logfile().list("Will assume that the reference allele is the preferred allele.");
		_siteTraverser.requireReference();
		_siteTraverser.filterRefN();
		_fileName = parameters().get("regions");
	} else {
		throw coretools::TUserError("Sites and preferred allele must be specified either using 'alleles' or 'bed'!");
	}
}

void TEstimateMutationLoad::run() {
	// traverse BAM file and store data
	//_traverseBAMWindows();
	_traverseSites();

	// check if sufficient sites
	coretools::user_assert(_sites.size() != 0, "No sites were kept after traversing BAM file!");

	// now run estimation
	MutationLoad::TMutationLoadEMPrior prior(_sites);
	MutationLoad::TMutationLoadLatentVariable latentVar(_sites);

	stattools::TEM<MutationLoad::PrecisionType, MutationLoad::NumStatesType, MutationLoad::LengthType> EM(prior,
	                                                                                                      latentVar);

	std::vector<MutationLoad::LengthType> chunkEnds = {_sites.size()};
	EM.runEM(chunkEnds);

	// write output file
	std::string filename = _siteTraverser.outputName() + "_mutationLoad.txt";
	coretools::TOutputFile out(filename, {"BAM", "Alleles", "Pi_rr", "Pi_ra", "Pi_aa", "Pi_ab"});
	out.writeln(_siteTraverser.bamFile().filename(), _fileName, prior.getPi());
}

} // end namespace GenomeTasks
