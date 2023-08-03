#include "TModel.h"
#include "TWithPMD.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Strings/stringFunctions.h"
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
	if (pmdString.empty() || pmdString == TNoPMD::name) return new TNoPMD;

	// Possibilities:
	// only Function: Exponential | Empiric
	// only Function perLength: Exponential[30] | Empiric[30]
	// Everything: 
	// singleStrand:Empiric[...]:Empiric[...]
	// doubleStrand:Exponential[30,0.1,0.1,0.05]:Exponential[40,0.2,0.3,0.07]
	// doubleStrand[30]:Exponential[30,0.1,0.1,0.05]:Exponential[40,0.2,0.3,0.07]

	TSplitter spl(pmdString, ':');
	const auto front = spl.front(); // either function-name or strand-name
	spl.popFront();

	if (spl.empty()) {  // only function
		if (front.back() == ']') {
			const auto function = coretools::str::readBefore(front, '[');
			auto from           = coretools::str::readAfter(front, '[');
			from.remove_suffix(1);

			if (from.empty()) return new TWithPMD<true>(function, size_t{30});
			/*else*/ return new TWithPMD<true>(function, fromString<size_t, true>(from));
		} else {
			return new TWithPMD<false>(front);
		}
	}
	// else

	const auto function5 = spl.front();
	spl.popFront();
	if (spl.empty()) { UERROR("You need to specify two function, 5' and 3'!"); }
	const auto function3 = spl.front();

	if (front.back() == ']') {
		const auto strand = getStrand(coretools::str::readBefore(front, '['));
		auto from         = coretools::str::readAfter(front, '[');
		from.remove_suffix(1);
		if (from.empty()) return new TWithPMD<true>(function5, function3, strand, size_t{30});
		/*else*/ return new TWithPMD<true>(function5, function3, strand, fromString<size_t, true>(from));
	} else {
		const auto strand    = getStrand(front);
		return new TWithPMD<false>(function5, function3, strand);
	}
}
} // namespace GenotypeLikelihoods::PMD
