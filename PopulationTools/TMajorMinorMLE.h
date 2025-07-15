#ifndef TMAJORMINORMLE_H_
#define TMAJORMINORMLE_H_

#include <cstddef>

#include "TMajorMinor.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/TFrequencies.h"
#include "genometools/VCF/VCF.h"


namespace PopulationTools::MajorMinor {

class TMLE {
private:
	template<size_t N>
	static TMMData _estimate(coretools::TConstView<genometools::TGLFEntry> data, double maxF,
							 const std::array<genometools::AllelicCombination, N>& usedAllelicCombinations, genometools::Base ref,
							coretools::Probability minMAF, coretools::PhredInt minVariantQuality) {
		using coretools::Log10Probability;
		// calculate L10L for each allelic combination
		genometools::TGenotypeFrequencies bestFreqs;
		auto bestL  = coretools::Log10Probability::lowest();
		auto bestAC = genometools::AllelicCombination::min;

		for (const auto ac : usedAllelicCombinations) {
			auto GLs = fill(data, ac);
			genometools::TGenotypeFrequencies freqs;
			freqs.estimate<false>(GLs, GLs.size(), maxF);
			auto L = freqs.calculateLog10Likelihood(GLs, GLs.size());
			if ((L > bestL) || (L == bestL && coretools::instances::randomGenerator().getRand() > 0.5)) {
				bestL     = L;
				bestFreqs = freqs;
				bestAC    = ac;
			}
		}
		genometools::Base major, minor;

		// which one is major?
		if (bestFreqs.alleleFrequency() < 0.5) {
			major = first(bestAC);
			minor = second(bestAC);
		} else {
			major = second(bestAC);
			minor = first(bestAC);
		}

		// determine variant quality
		Log10Probability LL_fixed = LLFixedAllele(data, major);
		const coretools::PhredInt variantQuality{LL_fixed > bestL ? Log10Probability(0.0)
																			   : Log10Probability(LL_fixed - bestL)};

		if (bestFreqs.MAF() < minMAF || variantQuality < minVariantQuality) {
			return failedTMMData;
		}

		// ensure return value uses major = ref if ref was provided
		if (minor == ref) { // cannot happen if base == ref
			minor = major;
			major = ref;
		}

		return {true, bestFreqs.MAF(), variantQuality, major, minor};
	}

public:

	static TMMData estimate(coretools::TConstView<genometools::TGLFEntry> data, double maxF,
							coretools::Probability minMAF, coretools::PhredInt minVariantQuality) {
		// calculate L10L for each allelic combination
		using AC = genometools::AllelicCombination;
		constexpr std::array usedAllelicCombinations = {AC::AC, AC::AG, AC::AT, AC::CG, AC::CT, AC::GT};
		return _estimate(data, maxF, usedAllelicCombinations, genometools::Base::N, minMAF, minVariantQuality);
	}

	static TMMData estimate(coretools::TConstView<genometools::TGLFEntry> data, double maxF,
							genometools::Base ref, coretools::Probability minMAF,
							coretools::PhredInt minVariantQuality) {
		// calculate L10L for each allelic combination
		if (ref == genometools::Base::N) return estimate(data, maxF, minMAF, minVariantQuality);

		const auto usedAllelicCombinations = useAllelicCombinationsThatContain(ref);
		return _estimate(data, maxF, usedAllelicCombinations, ref, minMAF, minVariantQuality);
	}

	static TMMData estimate(coretools::TConstView<genometools::TGLFEntry> data, double maxF,
							genometools::Base ref, genometools::Base alt, coretools::Probability minMAF,
							coretools::PhredInt minVariantQuality) {
		const auto usedAllelicCombinations = std::array{allelicCombination(ref, alt)};
		return _estimate(data, maxF, usedAllelicCombinations, ref, minMAF, minVariantQuality);
	}
};

} // namespace PopulationTools::MajorMinor

#endif
