#include "TSafEstimator.h"
#include "TGlfMultiReader.h"
#include "TSafFile.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/mathFunctions.h"
#include <algorithm>
#include <cmath>

namespace PopulationTools {

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::LogProbability;
using coretools::Probability;
using genometools::Base;
using genometools::genotype;
using coretools::TStrongArray;
using coretools::P;

namespace impl {
		
double sumInLog(double logA, double logB) {
	if (!std::isfinite(logA)) return logB;
	if (!std::isfinite(logB)) return logA;
	return logA > logB
		? std::log(1 + std::exp(logB - logA)) + logA
		: std::log(1 + std::exp(logA - logB)) + logB;
}

auto logChoose(size_t N) {
	std::vector<double> res(N + 1);
	for (size_t i = 0; i < N; ++i) { res[i] = coretools::chooseLog(N, i); }
	return res;
}
}


TSafEstimator::TSafEstimator() {
	_glfReader.addReference(parameters().get<std::string>("fasta"));
	_glfReader.openGLFs();
	_glfReader.setAllActive();

	const auto minSamples = parameters().get<size_t>("minSamplesWithData", 1);
	_glfReader.minSamplesWithData(minSamples);

	logfile().list("Will only print sites for which at least ", minSamples,
				   " samples have data. (parameter minSamplesWithData)");

	const size_t nSamples = _glfReader.numActiveSamples();
	_logProbs.assign(2*nSamples + 1, LogProbability::lowest());

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
	const static std::vector<double> logChoose = impl::logChoose(2*nSamples);

	_logProbs.clear();

	std::vector<double> probs(2*nSamples + 1, 0.);
	const auto mama = genotype(major, major);
	for (auto minor : minors[major]) {
		std::vector<double> hs(2 * nSamples + 1, P(0.));
		const auto mami = genotype(major, minor);
		const auto mimi = genotype(minor, minor);

		const auto &front = data.front();
		hs[0]             = Probability(front[mama]);
		hs[1]             = 2 * Probability(front[mami]);
		hs[2]             = Probability(front[mimi]);

		for (size_t d = 1; d < data.size(); ++d) {
			for (size_t j = 2 * (d + 1); j > 1; --j) {
				hs[j] = Probability(data[d][mama]) * hs[j] + 2 * Probability(data[d][mami]) * hs[j - 1] +
				        Probability(data[d][mimi]) * hs[j - 2];
			}
			hs[1] = Probability(data[d][mama]) * hs[1] + 2 * Probability(data[d][mami]) * hs[0];
			hs[0] = Probability(data[d][mama]) * hs[0];
		}
		// for non-data
		for (size_t d = data.size(); d < nSamples; ++d) {
			for (size_t j = 2 * d; j > 1; --j) { hs[j] = hs[j] + 2 * hs[j - 1] + hs[j - 2]; }
			hs[1] = hs[1] + 2 * hs[0];
			hs[0] = hs[0];
		}

		if (_logProbs.empty()) {
			for (size_t j = 0; j < hs.size(); ++j) { _logProbs.push_back(std::log(hs[j]) - logChoose[j]); }
		} else {
			for (size_t j = 0; j < _logProbs.size(); ++j) {
				_logProbs[j] = impl::sumInLog(_logProbs[j], std::log(hs[j]) - logChoose[j]);
			}
		}
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
