/*
 * TEstimateMutationLoad.cpp
 */

#include "TEstimateMutationLoad.h"
#include "coretools/traits.h"
#include "stattools/EM/TEM.h"

using coretools::instances::logfile;
using genometools::Genotype;
using GenotypeLikelihoods::TGenotypeLikelihoods;

namespace GenomeTasks {

namespace MutationLoad {

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
void TEstimateMutationLoad::_addSite(const GenotypeLikelihoods::TSite &site, const genometools::Base PreferredBase) {
	if (!site.empty()) {
		GenotypeLikelihoods::TGenotypeLikelihoods genoLik =
		    _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(site);
		_sites.emplace_back(genoLik, PreferredBase);
	}
}

void TEstimateMutationLoad::_handleWindow() {
	// adding sites to estimator
	logfile().listFlushTime("Calculating genotype likelihoods and storing data ...");
	try {
		if (_parseFromBed) {
			// get sites from bed file and alleles from reference
			auto it = _bedFile.lower_bound(_window);
			while (it != _bedFile.end() && _window.overlaps(*it)) {
				for (genometools::TGenomePosition s = std::max(it->from(), _window.from());
				     s < it->to() && s < _window.to(); ++s) {
					const GenotypeLikelihoods::TSite &site = _window[s - _window.from()];
					_addSite(site, site.refBase);
				}
				++it;
			}
		} else {
			// get sites and alleles from site subset
			auto thesePositions = _subsetMonomorphic->getPositionInWindow(_window);
			for (auto &it : thesePositions) { _addSite(_window[it - _window.from()], it.ref()); }
		}
	} catch (...) {
		UERROR("Failed to allocate sufficient memory to store the data for so many sites. Consider using fewer sites.");
	}
	logfile().doneTime();
};

TEstimateMutationLoad::TEstimateMutationLoad() : TGenome_windows() {
	using coretools::instances::logfile;
	using coretools::instances::parameters;
	// Two ways to read positions and preferred alleles:
	//  1) from an alleles file (chr, pos, allele)
	//  2) from a BED file and the reference
	if (parameters().exists("alleles")) {
		_openSiteSubset("alleles", false);
		_parseFromBed = false;
	} else if (parameters().exists("bed")) {
		logfile().startIndent("Limiting analysis to sites listed in BED file:");
		// open reference
		logfile().list("Will assume that the reference allele is the preferred allele.");
		_openReference(true);
		// parse BED
		_bedFileName = parameters().get("bed");
		logfile().listFlush("Reading BED file '", _bedFileName, "' (parameter 'bed') ...");
		_bedFile.add(_bedFileName, _bamFile.chromosomes());
		logfile().done();
		logfile().conclude("Read ", _bedFile.size(), " sites on ", _bedFile.numChromosomesWithWindows(),
		                   " chromosomes.");
		_parseFromBed = true;
		logfile().endIndent();
	} else {
		UERROR("Sites and preferred allele must be specified either using 'alleles' or 'bed'!");
	}
};

void TEstimateMutationLoad::run() {
	// traverse BAM file and store data
	_traverseBAMWindows();

	// check if sufficient sites
	if (_sites.size() == 0) { UERROR("No sites were kept after traversing BAM file!"); }

	// now run estimation
	MutationLoad::TMutationLoadEMPrior prior(_sites);
	MutationLoad::TMutationLoadLatentVariable latentVar(_sites);

	stattools::TEM<MutationLoad::PrecisionType, MutationLoad::NumStatesType, MutationLoad::LengthType> EM(prior,
	                                                                                                      latentVar);

	std::vector<MutationLoad::LengthType> chunkEnds = {_sites.size()};
	EM.runEM(chunkEnds);

	// write output file
	std::string filename = _outputName + "_mutationLoad.txt";
	coretools::TOutputFile out(filename, {"BAM", "Alleles", "Pi_rr", "Pi_ra", "Pi_aa", "Pi_ab"});
	if (_parseFromBed) {
		out.writeln(_bamFile.filename(), _bedFileName, prior.getPi());
	} else {
		out.writeln(_bamFile.filename(), _subsetMonomorphic->filename(), prior.getPi());
	}
}

} // end namespace GenomeTasks
