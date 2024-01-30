#include "TContextErrors.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/toString.h"
#include "coretools/enum.h"

namespace GenomeTasks{
using genometools::Base;
using coretools::str::toString;

TContextErrors::TContextErrors()
	: _out(_genome.outputName() + "_fromto.txt.gz"), _depth(coretools::instances::parameters().get("depth", 10)) {
	_windows.requireReference();
}

void TContextErrors::run() {
	_out.writeHeader({"Previous", "Major", toString(_depth, ":0"), toString(_depth - 1, ":1")});
	_traverseBAMWindows();
}


void TContextErrors::_handleWindow(GenotypeLikelihoods::TWindow& window) {
	using coretools::TStrongArray;
	TStrongArray<TStrongArray<size_t, Base>, Base, coretools::index(Base::N) + 1> counts;
	for (size_t i = 0; i < window.size(); ++i) {
		const auto & site = window[i];
		if (site.depth() != _depth) continue;

		TStrongArray<size_t, Base> counts_i{};

		for (const auto& b: site) {
			++counts_i[b.base];

		}
		auto bMajor = Base::N;
		auto bMinor = Base::N;
		size_t nMajor = 0;
		size_t nMinor = 0;
		for (auto b = Base::min; b < Base::max; ++b) {
			const auto c = counts_i[b];
			if (c > nMajor) {
				nMajor = c;
				bMajor = b;
			} else if (c > nMinor) {
				nMinor = c;
				bMinor = b;
			}
		}
		if (nMajor == _depth) {
			
		} else if (nMinor == 1) {
			
		}
	}
}
} // namespace GenomeTasks
