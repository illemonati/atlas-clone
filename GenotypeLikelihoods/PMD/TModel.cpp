#include "TModel.h"
#include "TWithPMD.h"
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

		// default values
		size_t from = 30;
		size_t to   = 120;

		std::string_view interval = spl.front();
		interval.remove_prefix(1);
		interval.remove_suffix(1);
		TSplitter from_to(interval, ',');
		if (!from_to.empty()) {
			from = fromStringCheck<size_t>(from_to.front());
			from_to.popFront();
			if (!from_to.empty()) to = fromStringCheck<size_t>(from_to.front());
		}

		spl.popFront();
		const auto function5 = spl.front();
		spl.popFront();
		const auto function3 = spl.empty() ? function5 : spl.front();
		if (strand == Strand::Single) return new TWithPMD<Strand::Single, true>(function5, function3, from, to);
		/*else*/
		return new TWithPMD<Strand::Single, true>(function5, function3, from, to);
	} else {
		const auto function5 = spl.front();
		spl.popFront();
		const auto function3 = spl.empty() ? function5 : spl.front();
		if (strand == Strand::Single) return new TWithPMD<Strand::Single, false>(function5, function3);
		/*else*/
		return new TWithPMD<Strand::Double, false>(function5, function3);
	}
}
} // namespace GenotypeLikelihoods::PMD
