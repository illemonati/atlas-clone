#include "TSafEstimator.h"
#include "TGlfMultiReader.h"
#include "TSafFile.h"
#include "coretools/Main/TParameters.h"
#include <algorithm>
#include <cmath>
#include <limits>

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

constexpr bool multiModal(double mama, double mami, double mimi) {
  return mami < mama && mami < mimi;
}
} // namespace impl

TSafEstimator::TSafEstimator() {
	_glfReader.addReference(parameters().get<std::string>("fasta"));
	_glfReader.openGLFs();
	_glfReader.setAllActive();
	const size_t nSamples = _glfReader.numActiveSamples();

	logfile().list("Will estimate Site Allele Frequency for ", nSamples, " samples.");

	const auto minSamples = parameters().get<size_t>("minSamplesWithData", 1);
	_glfReader.minSamplesWithData(minSamples);

	logfile().list("Will print sites for which at least ", minSamples,
				   " samples have data. (parameter minSamplesWithData)");

	_logProbs.assign(2*nSamples + 1, LogProbability::lowest());

}

void TSafEstimator::_iterate(const TGenotypeLikelihoodsAllCombinationsVector &data, Base major) {
	assert(major != Base::N);
	constexpr TStrongArray<std::array<Base, 3>, Base> minors = []() {
		TStrongArray<std::array<Base, 3>, Base> minors{};
		minors[Base::A] = {Base::C, Base::G, Base::T};
		minors[Base::C] = {Base::A, Base::G, Base::T};
		minors[Base::G] = {Base::A, Base::C, Base::T};
		minors[Base::T] = {Base::A, Base::C, Base::G};
		return minors;
	}();
	const size_t nSamples = _glfReader.numActiveSamples();

	_logProbs.clear();

	const auto mama = genotype(major, major);


	constexpr auto tol    = 1e-9;
	constexpr auto logTol = -20.72326583694641;

	size_t upper = 2;
	double lPmax = logTol;

	for (auto minor : minors[major]) {
		std::vector<double> hs(2 * nSamples + 1, P(0.));
		const auto mami = genotype(major, minor);
		const auto mimi = genotype(minor, minor);

		const auto &front  = data.front();
		const double pMama = Probability(front[mama]);
		const double pMami = Probability(front[mami]);
		const double pMimi = Probability(front[mimi]);
		hs[0]              = pMama;
		hs[1]              = pMami; // not 2*pMami, see eq (6) Han et al
		hs[2]              = pMimi;

		double max      = std::max({hs[0], hs[1], hs[2]});
		double sumLog   = 0.;
		bool multiModal = impl::multiModal(pMama, pMami, pMimi);

		for (size_t d = 1; d < nSamples; ++d) {
			double pMama = 1.;
			double pMami = 1.;
			double pMimi = 1.;
			if (d < data.size()) {
				pMama = Probability(data[d][mama]);
				pMami = Probability(data[d][mami]);
				pMimi = Probability(data[d][mimi]);
			}
			if (impl::multiModal(pMama, pMami, pMimi)) multiModal = true; // once true, always true

			double nextMax = 0.;
			size_t N = 2 * (d + 1);

			if (multiModal) {
				upper = N + 1;
			} else {
				upper = std::min(N, upper + 2);
				for (size_t j = upper; j > 1; --j) {
					hs[j] = ((N-j)*(N-j-1)*pMama*hs[j]
							 + 2*j*(N-j)*pMami*hs[j - 1]
							 + j*(j-1)*pMimi*hs[j - 2])/max;
					if (hs[j] < tol) {
						hs[j] = 0;
						--upper;
					} else {
						break;
					}
				}
			}
			for (size_t j = upper - 1; j > 1; --j) {
				hs[j] = ((N-j)*(N-j-1)*pMama*hs[j]
					  + 2*j*(N-j)*pMami*hs[j - 1]
						 + j*(j-1)*pMimi*hs[j - 2])/max;
				nextMax = std::max(hs[j], nextMax);
			}
			hs[1] = ((N-1)*(N-2)*pMama * hs[1] + 2*(N-1) * pMami * hs[0])/max;
			hs[0] = ((N)*(N-1)*pMama * hs[0])/max;
			sumLog += std::log(max);
			max = std::max({hs[1], hs[0], nextMax});
		}

		if (_logProbs.empty()) {
			for (size_t j = 0; j < hs.size(); ++j) {
				if (hs[j] > 0) {
					_logProbs.push_back(std::log(hs[j]) + sumLog);
					lPmax = std::max(_logProbs.back(), lPmax);
				} else {
					_logProbs.push_back(-std::numeric_limits<double>::infinity());
				}
			}
		} else {
			for (size_t j = 0; j <hs.size(); ++j) {
				if (hs[j] > 0) {
					_logProbs[j] = impl::sumInLog(_logProbs[j], std::log(hs[j]) + sumLog);
					lPmax = std::max(_logProbs[j], lPmax);
				}
			}
		}
	}
	// remove upper values
	while(!_logProbs.empty() && (!std::isfinite(_logProbs.back()) || _logProbs.back() - lPmax < logTol)) _logProbs.pop_back();
	for (auto& lP: _logProbs) {
		lP -= lPmax;
	}
	// remove lower values
	for (_lower = 0; _lower < _logProbs.size() - 1; ++_lower) { // ensure _lower remains < size
		if (_logProbs[_lower] > logTol) break;
	}
}

void TSafEstimator::run() {
	const auto nSamples   = _glfReader.numActiveSamples();
	const auto oFile      = parameters().get("out", "saf");

	TSafFile safFile(oFile, nSamples);
	logfile().list("Will write saf-file with prefix ", oFile, ".");


	coretools::TTimer timer;
	const auto dCounter = std::min<size_t>(10000000 / nSamples, 1000000);
	size_t counter        = 0;
	size_t counterF       = 0;
	size_t nextPrint      = dCounter;

	logfile().startIndent("Parsing through glf files:");
	for (auto ids = _glfReader.readWindow(); !ids.empty(); ids = _glfReader.readWindow()) {
		const auto refs = _glfReader.refView();
		for (const auto iW : ids) {
			if (refs[iW] == Base::N) {
				++counterF;
				continue;
			}

			_iterate(_glfReader.data(iW), refs[iW]);
			if (!_logProbs.empty()) safFile.write(_glfReader.curChrName(), _glfReader.position(iW).position(), _logProbs, _lower);
		}

		// report progress
		counter  += ids.size();
		if (counter >= nextPrint) {
			logfile().list("Parsed ", nextPrint, " positions in ", timer.formattedTime(), ".");
			while (nextPrint <= counter) nextPrint += dCounter;
		}
	}
	logfile().list("Reached end of glf files!");
	logfile().list("Parsed a total of ", counter, " positions, filtered: ", counterF, " (", (100.*counterF)/counter, "%).");
}
}
