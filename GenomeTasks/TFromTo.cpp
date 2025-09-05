#include "TFromTo.h"
#include "TSequencedData.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "genometools/Genotypes/Base.h"

namespace GenomeTasks{
using coretools::instances::randomGenerator;
using genometools::Base;

TFromTo::TFromTo() : _out(_siteTraverser.outputName() + "_fromto.txt.gz") {
	_siteTraverser.requireReference();
	_siteTraverser.requireSingleBAM();
	_out.writeHeader({"Chr", "Pos", "Depth", "dist5_1", "dist3_1", "dist5_2", "dist3_2", "reveresed1", "reversed2"});
}

void TFromTo::run() { _traverseSites(); }

void TFromTo::_traverseSites() {
	using BAM::Flags;
	for (; !_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		for (; !_siteTraverser.endOfCurChr(); _siteTraverser.nextSite()) {
			const auto &site = _siteTraverser.site();
			const auto r     = site.refBase;
			if (r == Base::C || r == Base::G) continue;

			std::vector<BAM::TSequencedData> data;
			for (const auto &d : site) {
				if (d.base != r) data.push_back(d);
			}
			if (data.size() < 2 || data.size() > site.depth() / 2) continue;

			const auto i1 = randomGenerator().getRand(size_t{}, data.size());
			auto i2       = i1;
			while (i2 == i1) i2 = randomGenerator().getRand(size_t{}, data.size());

			_out.writeln(_siteTraverser.curChr().name(), _siteTraverser.position().position(), site.depth(), data[i1].distFrom5.linear(),
			             data[i1].distFrom3.linear(), data[i2].distFrom5.linear(), data[i2].distFrom3.linear(),
			             data[i1].get<Flags::ReversedStrand>(), data[i2].get<Flags::ReversedStrand>());
		}
	}
}
} // namespace GenomeTasks
