#include "TFromTo.h"
#include "TSequencedBase.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "genometools/GenotypeTypes.h"

namespace GenomeTasks{
using coretools::instances::randomGenerator;
using genometools::Base;

TFromTo::TFromTo() : _out(_genome.outputName() + "_fromto.txt.gz") {
	_windows.requireReference();
	_out.writeHeader({"Chr", "Pos", "Depth", "dist5_1", "dist3_1", "dist5_2", "dist3_2", "reveresed1", "reversed2"});
}

void TFromTo::run() { _traverseBAMWindows(); }


void TFromTo::_handleWindow(GenotypeLikelihoods::TWindow& window) {
	for (size_t i = 0; i < window.size(); ++i) {
		const auto & site = window[i];
		if (site.empty()) continue;

		const auto r = site.refBase;
		if (r == Base::C || r == Base::G) continue;

		const auto pos = window.position(i);
		std::vector<BAM::TSequencedBase> bases;
		for (const auto& b: site) {
			if (b.base != r) bases.push_back(b);
		}
		if (bases.size() < 2 || bases.size() > site.depth()/2) continue;

		const auto i1 = randomGenerator().getRand(size_t{}, bases.size());
		auto i2       = i1;
		while (i2 == i1) i2 = randomGenerator().getRand(size_t{}, bases.size());

		_out.writeln(window.chrName(), pos.position(), site.depth(), bases[i1].distFrom5Prime, bases[i1].distFrom3Prime,
					 bases[i2].distFrom5Prime, bases[i2].distFrom3Prime, bases[i1].isReverseStrand(), bases[i2].isReverseStrand());
	}
}
} // namespace GenomeTasks
