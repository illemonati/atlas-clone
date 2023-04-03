#include "TModel.h"
#include "TPerReadGroup.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"
#include <utility>

namespace GenotypeLikelihoods::PMD {
using genometools::Base;
using namespace coretools::str;

TBaseProbabilities TNoPMD::massFunction(genometools::Genotype g, const TBaseLikelihoods &baseLikelihoodsNoPMD) {
	using namespace genometools;
	const Base a = first(g);
	const Base b = second(g);
	TBaseLikelihoods mf{0};
	mf[a] = baseLikelihoodsNoPMD[a];
	mf[b] = baseLikelihoodsNoPMD[b];
	return TBaseProbabilities::normalize(mf);
}

Strand getStrand(std::string_view s) {
		constexpr coretools::TStrongArray<std::string_view, Strand> strands{{"singleStrand", "doubleStrand"}};
		for (auto strand = Strand::min; strand < Strand::max; ++strand) {
			if (s == strands[strand]) return strand;
		}
		UERROR("Wrong strand defined. Use either ", strands.front(), " or ", strands.back(), "!");
}


TModel *makeType(std::string_view pmdString) {
	// split by ':'

	TSplitter spl(pmdString, ':');

	if (spl.front() == TNoPMD::name) return new TNoPMD;

	const auto strand = getStrand(spl.front());
	spl.popFront();

	if (spl.front().front() == '[') {
		// per fragment length
		// example: doubleStrand:[30,120]:Exponential
		TSplitter from_to(spl.front().substr(1, spl.front().size() - 1), ',');
		const auto from = fromStringCheck<double>(from_to.front());
		from_to.popFront();
		const auto to   = fromStringCheck<double>(from_to.front());

		spl.popFront();
		const auto function5 = spl.front();
		spl.popFront();
		const auto function3 = spl.empty() ? function5 : spl.front();
		if (strand == Strand::Single) return new TPerReadGroup<Strand::Single>(function5, function3);
		/*else*/
		return new TPerReadGroup<Strand::Double>(function5, function3);
	} else {
		const auto function5 = spl.front();
		spl.popFront();
		const auto function3 = spl.empty() ? function5 : spl.front();
		if (strand == Strand::Single) return new TPerReadGroup<Strand::Single>(function5, function3);
		/*else*/
		return new TPerReadGroup<Strand::Double>(function5, function3);
	}
}
} // namespace GenotypeLikelihoods::PMD
