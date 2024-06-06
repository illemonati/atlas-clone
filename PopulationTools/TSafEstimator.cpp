#include "TSafEstimator.h"
#include "TGlfMultiReader.h"
#include "TSafFile.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/devtools.h"
#include <algorithm>

namespace PopulationTools {

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::LogProbability;
using coretools::Probability;
using genometools::Base;
using genometools::genotype;
using coretools::TStrongArray;
using coretools::P;

TSafEstimator::TSafEstimator() {
	_glfReader.addReference(parameters().get<std::string>("fasta"));
	_glfReader.openGLFs();
	_glfReader.setAllActive();

	const auto minSamples = parameters().get<size_t>("minSamplesWithData", 1);
	_glfReader.minSamplesWithData(minSamples);

	logfile().list("Will only print sites for which at least ", minSamples,
				   " samples have data. (parameter minSamplesWithData)");


	_logProbs.assign(2*_glfReader.numActiveSamples() + 1, LogProbability::lowest());
}

void TSafEstimator::_iterate(const TGenotypeLikelihoodsAllCombinationsVector &data, Base major) {
	assert(major == Base::N);
	constexpr TStrongArray<std::array<Base, 3>, Base> minors = []() {
		TStrongArray<std::array<Base, 3>, Base> minors{};
		minors[Base::A] = {Base::C, Base::G, Base::T};
		minors[Base::C] = {Base::A, Base::G, Base::T};
		minors[Base::G] = {Base::A, Base::C, Base::T};
		minors[Base::T] = {Base::A, Base::C, Base::G};
		return minors;
	}();
	const size_t nSamples = _glfReader.numActiveSamples();

	std::vector<double> probs(2*nSamples + 1, P(0.));
	const auto mama = genotype(major, major);
	for (auto minor: minors[major]) {
		std::vector<double> hs(2*nSamples + 1, P(0.));
		const auto mami = genotype(major, minor);
		const auto mimi = genotype(minor, minor);

		const auto& front = data.front();
		hs[0] = Probability(front[mama]);
		hs[1] = 2*Probability(front[mami]);
		hs[2] = Probability(front[mimi]);

		// TODO: more than one sample here

		for (size_t j = 0; j < probs.size(); ++j) {
			probs[j] += hs[j]/coretools::choose(2*nSamples, j);
		}
	}


	const auto max = *std::max_element(probs.begin(), probs.end());
	for (size_t i = 0; i < probs.size(); ++i) {
		_logProbs[i] = coretools::logP(std::log(probs[i]/max));
	}
}

void TSafEstimator::run() {
	const auto nSamples   = _glfReader.numActiveSamples();
	const auto oFile      = parameters().get("out", "saf");

	TSafFile safFile(oFile, nSamples);
	logfile().list("Will write saf-file with prefix ", oFile, ".");

	for (auto ids = _glfReader.readWindow(); !ids.empty(); ids = _glfReader.readWindow()) {
		const auto refs = _glfReader.refView();
		for (const auto iW : ids) {
			if (refs[iW] == Base::N) continue;

			_iterate(_glfReader.data(iW), refs[iW]);
			safFile.write(_glfReader.curChrName(), _glfReader.position(iW).position(), _logProbs);
		}
	}
}
}
