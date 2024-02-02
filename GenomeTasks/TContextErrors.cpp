#include "TContextErrors.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Strings/toString.h"
#include "coretools/enum.h"
#include <algorithm>

namespace GenomeTasks{
using genometools::Base;
using coretools::str::toString;
using coretools::instances::logfile;
using coretools::instances::parameters;

TContextErrors::TContextErrors()
	: _out(_genome.outputName() + "_contextErrors.txt.gz"), _depth(parameters().get("depth", 10)) {
}

void TContextErrors::run() {
	logfile().list("Running Context Statistics for sites with depth ", _depth, ", counting the occurance of ", _depth, ":0 and ", _depth - 1, ":1.");
	_traverseBAMWindows();
	logfile().list("The count of major Bases [A, C, G, T] are", _major, ".");;
	logfile().list("The count of minor Bases [A, C, G, T, N] are", _minor, ".");;

	_out.writeHeader({"Major", "Previous", toString(_depth, ":0"), toString(_depth - 1, ":1"), "Fraction"});
	logfile().list("Writing Context Staticstics to ", _out.name(), ".");;
	for (auto maj = Base::min; maj < Base::max; ++maj) {
		for (auto prev = Base::min; prev <= Base::max; ++prev) {
			const auto nE = _counts[maj][prev].noError;
			const auto oE = _counts[maj][prev].oneError;
			const double f = nE > 0 ? double(oE)/nE : 0.;
			_out.writeln(maj, prev, nE, oE, f);
		}
	}
}


void TContextErrors::_handleWindow(GenotypeLikelihoods::TWindow& window) {
	using coretools::TStrongArray;
	using coretools::instances::randomGenerator;
	for (size_t i = 0; i < window.size(); ++i) {
		const auto & site = window[i];
		if (site.depth() != _depth) continue;

		TStrongArray<size_t, Base> counts_now{};
		TStrongArray<size_t, Base, coretools::index(Base::N) + 1> counts_prev{};

		for (const auto& b: site) {
			if (b.base == Base::N) continue;
			++counts_now[b.base];
			++counts_prev[b.previousBase];

		}
		auto bMajor = Base::N;
		auto bMinor = Base::N;
		auto bPrev  = Base::N;
		size_t nMajor = 0;
		size_t nMinor = 0;
		size_t nPrev  = 0;
		for (auto b = Base::min; b < Base::max; ++b) {
			const auto c_p = counts_prev[b];
			if (c_p > nPrev) {
				nPrev = c_p;
				bPrev = b;
			}

			const auto c_i = counts_now[b];
			if (c_i > nMajor) {
				nMinor = nMajor;
				bMinor = bMajor;
				nMajor = c_i;
				bMajor = b;
			} else if (c_i > nMinor) {
				nMinor = c_i;
				bMinor = b;
			}
		}
		++_major[bMajor];
		++_minor[bMinor];
		//const auto it = std::max_element(counts_now.begin(), counts_now.end());

		if (nMajor == _depth && nMinor == 0) {
			++_counts[bMajor][bPrev].noError;
		} else if (nMajor == _depth - 1 && nMinor == 1) {
			++_counts[bMajor][bPrev].oneError;
		}
	}
}
} // namespace GenomeTasks
