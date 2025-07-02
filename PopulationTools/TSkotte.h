#ifndef TSKOTTE_H_
#define TSKOTTE_H_

#include "TMajorMinor.h"
#include "coretools/Containers/TDualStrongArray.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Types/probability.h"
#include "coretools/Math/TSumLog.h"

#include "genometools/GLF/GLF.h"
#include "genometools/Genotypes/AllelicCombination.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/Genotype.h"
#include <cstddef>
namespace PopulationTools::MajorMinor {

class TSkotte {
	enum class HaploDiplo : size_t { min, first = min, second, homoFirst, het, homoSecond, max };

	static HaploDiplo _haploIndex(coretools::Probability pFirst, coretools::Probability pSecond) { return HaploDiplo(pSecond > pFirst); }

	static HaploDiplo _diploIndex(coretools::Probability pHomoFirst, coretools::Probability pHet, coretools::Probability pHomoSecond) {
		std::array<coretools::Probability, 3> es = {pHomoFirst, pHet, pHomoSecond};
		return HaploDiplo(2 + std::distance(es.begin(), std::max_element(es.begin(), es.end())));
	}

	static double _haploWeights(coretools::Probability pFirst, coretools::Probability pSecond,
								const coretools::TStrongArray<coretools::Probability, HaploDiplo> freqs) {
		const double weights0 = pFirst * freqs[HaploDiplo::first];
		const double weights1 = pSecond * freqs[HaploDiplo::second];
		const double sum      = weights0 + weights1;
		return weights0 / sum;
	}

	static std::pair<double, double> _diploWeights(coretools::Probability pHomoFirst, coretools::Probability pHet, coretools::Probability pHomoSecond,
												   const coretools::TStrongArray<coretools::Probability, HaploDiplo> freqs) {
		const double weights0 = pHomoFirst * freqs[HaploDiplo::homoFirst];
		const double weights1 = pHet * freqs[HaploDiplo::het];
		const double weights2 = pHomoSecond * freqs[HaploDiplo::homoSecond];
		const double sum_i    = 1. / (weights0 + weights1 + weights2);
		return {weights0 * sum_i, weights2 * sum_i};
	}

	template<bool hasHaploid, bool hasDiploid>
	static auto _iterate(coretools::TConstView<coretools::TDualStrongArray<coretools::Probability, genometools::Base, genometools::Genotype>> glfs,
						 genometools::AllelicCombination ac, double maxF, coretools::Probability minMAF) {
		using coretools::Probability;
		using coretools::P;
		using coretools::TStrongArray;

		constexpr bool hasBoth = hasHaploid && hasDiploid;

		const auto first      = genometools::first(ac);
		const auto second     = genometools::second(ac);
		const auto homoFirst  = genometools::homoFirst(ac);
		const auto het        = genometools::het(ac);
		const auto homoSecond = genometools::homoSecond(ac);

		// Initial guess
		coretools::TStrongArray<size_t, HaploDiplo> counts{0};
		for (const auto &g : glfs) {
			if constexpr (hasBoth) {
				if (g.isType<coretools::ABType::A>()) { // Haploid
					++counts[_haploIndex(g[first], g[second])];
				} else { // Diploid
					++counts[_diploIndex(g[homoFirst], g[het], g[homoSecond])];
				}
			} else if (hasHaploid) {
				++counts[_haploIndex(g[first], g[second])];
			} else if (hasDiploid) {
				++counts[_diploIndex(g[homoFirst], g[het], g[homoSecond])];
			}
		}

		constexpr Probability fMin{0.000001};
		TStrongArray<Probability, HaploDiplo> freqs{P(0.)};
		const auto nHaplo     = counts[HaploDiplo::first] + counts[HaploDiplo::second];
		const double nHaplo_1 = 1. / nHaplo;
		if constexpr (hasHaploid) {
			freqs[HaploDiplo::first]  = std::max<Probability>(P(counts[HaploDiplo::first] * nHaplo_1), fMin);
			freqs[HaploDiplo::second] = std::max<Probability>(P(counts[HaploDiplo::second] * nHaplo_1), fMin);
			const auto sum            = freqs[HaploDiplo::first] + freqs[HaploDiplo::second];
			freqs[HaploDiplo::first].scale(sum);
			freqs[HaploDiplo::second].scale(sum);
		}
		const auto nDiplo = counts[HaploDiplo::homoFirst] + counts[HaploDiplo::het] + counts[HaploDiplo::homoSecond];
		const double nDiplo_1 = 1. / nDiplo;
		if constexpr (hasDiploid) {
			freqs[HaploDiplo::homoFirst]  = std::max<Probability>(P(counts[HaploDiplo::homoFirst] * nDiplo_1), fMin);
			freqs[HaploDiplo::het]        = std::max<Probability>(P(counts[HaploDiplo::het] * nDiplo_1), fMin);
			freqs[HaploDiplo::homoSecond] = std::max<Probability>(P(counts[HaploDiplo::homoSecond] * nDiplo_1), fMin);
			const auto sum = freqs[HaploDiplo::homoFirst] + freqs[HaploDiplo::het] + freqs[HaploDiplo::homoSecond];
			freqs[HaploDiplo::homoFirst].scale(sum);
			freqs[HaploDiplo::het].scale(sum);
			freqs[HaploDiplo::homoSecond].scale(sum);
		}

		const double nHaplo2nDiplo_1 = 1. / (nHaplo + 2 * nDiplo);

		// iterate
		constexpr size_t maxIter = 1000;

		Probability aF  = P((nDiplo * (freqs[HaploDiplo::het] + 2.0 * freqs[HaploDiplo::homoSecond]) +
							nHaplo * (double)freqs[HaploDiplo::second]) *
							nHaplo2nDiplo_1);
		Probability MAF = (aF < 0.5 ? aF : aF.complement());

		if (MAF < minMAF / 2) {
			// will be automaticall filtered out
			return std::make_tuple(Probability::lowest(), coretools::Log10Probability::lowest(), first, second);
		}

		for (size_t _ = 0; _ < maxIter; ++_) {
			// set genofreq
			double hplF0 = 0;
			double dplF0 = 0;
			double dplF2 = 0;

			for (const auto &g : glfs) {
				if constexpr (hasBoth) {
					if (g.isType<coretools::ABType::A>()) { // haploid
						hplF0 += _haploWeights(g[first], g[second], freqs);
					} else {
						const auto [w0, w2] = _diploWeights(g[homoFirst], g[het], g[homoSecond], freqs);
						dplF0 += w0;
						dplF2 += w2;
					}
				} else if (hasHaploid) {
					hplF0 += _haploWeights(g[first], g[second], freqs);
				} else if (hasDiploid) {
					const auto [w0, w2] = _diploWeights(g[homoFirst], g[het], g[homoSecond], freqs);
					dplF0 += w0;
					dplF2 += w2;
				}
			}
			double maxF_i{};
			if constexpr (hasHaploid) {
				hplF0 *= nHaplo_1;
				const auto hplF1 = 1.0 - std::min(1.0, hplF0);

				// check if we stop
				maxF_i = std::max(fabs(hplF0 - freqs[HaploDiplo::first]), fabs(hplF1 - freqs[HaploDiplo::second]));

				freqs[HaploDiplo::first]  = P(hplF0);
				freqs[HaploDiplo::second] = P(hplF1);
			}
			if constexpr (hasDiploid) {
				dplF0 *= nDiplo_1;
				dplF2 *= nDiplo_1;
				const auto dplF1 = 1.0 - std::min(1.0, (dplF0 + dplF2));
				// 1 - sum ensures range despite numeric inaccuracies

				// check if we stop
				if constexpr (hasHaploid) {
					maxF_i =
						std::max({maxF_i, fabs(dplF0 - freqs[HaploDiplo::homoFirst]),
								  fabs(dplF1 - freqs[HaploDiplo::het]), fabs(dplF2 - freqs[HaploDiplo::homoSecond])});
				} else {
					maxF_i = std::max({fabs(dplF0 - freqs[HaploDiplo::homoFirst]), fabs(dplF1 - freqs[HaploDiplo::het]),
									   fabs(dplF2 - freqs[HaploDiplo::homoSecond])});
				}

				freqs[HaploDiplo::homoFirst]  = P(dplF0);
				freqs[HaploDiplo::het]        = P(dplF1);
				freqs[HaploDiplo::homoSecond] = P(dplF2);
			}

			aF = P((nDiplo * (freqs[HaploDiplo::het] + 2.0 * freqs[HaploDiplo::homoSecond]) +
					nHaplo * (double)freqs[HaploDiplo::second]) *
				   nHaplo2nDiplo_1);

			MAF = aF < 0.5 ? aF : aF.complement();

			if (MAF < minMAF - maxF_i) {
				// will be automaticall filtered out
				return std::make_tuple(Probability::lowest(), coretools::Log10Probability::lowest(), first, second);
			}

			if (maxF_i < maxF) { break; }
		}

		const auto [major, minor] = [aF, first, second]() {
			if (aF < 0.5) {
				return std::make_tuple(first, second);
			} else {
				return std::make_tuple(second, first);
			}
		}();

		coretools::TSumLogProbability L{};
		for (const auto &g : glfs) {
			if constexpr (hasBoth) {
				if (g.isType<coretools::ABType::A>()) { // haploid
					L.add(g[first] * freqs[HaploDiplo::first] + g[second] * freqs[HaploDiplo::second]);
				} else {
					L.add(g[homoFirst] * freqs[HaploDiplo::homoFirst] + g[het] * freqs[HaploDiplo::het] +
						  g[homoSecond] * freqs[HaploDiplo::homoSecond]);
				}
			} else if (hasHaploid) {
				L.add(g[first] * freqs[HaploDiplo::first] + g[second] * freqs[HaploDiplo::second]);
			} else if (hasDiploid) {
				L.add(g[homoFirst] * freqs[HaploDiplo::homoFirst] + g[het] * freqs[HaploDiplo::het] +
					  g[homoSecond] * freqs[HaploDiplo::homoSecond]);
			}
		}
		constexpr double l10_1 = 0.43429448190325176;
		const coretools::Log10Probability bestL{L.getSum() * l10_1};

		return std::make_tuple(MAF, bestL, major, minor);
	}

	template<size_t N>
	static TMMData _estimate(coretools::TConstView<genometools::TGLFEntry> data, double maxF,
							 const std::array<genometools::AllelicCombination, N>& usedAllelicCombinations, genometools::Base ref,
							coretools::Probability minMAF, coretools::PhredInt minVariantQuality) {
		using coretools::Probability;
		using coretools::Log10Probability;
		using coretools::P;
		using coretools::TStrongArray;
		using genometools::Base;
		using genometools::Genotype;

		static std::vector<coretools::TDualStrongArray<coretools::Probability, genometools::Base, genometools::Genotype>> glfs;
		glfs.clear();

		coretools::TStrongArray<coretools::TSumLogProbability, genometools::AllelicCombination> Ls{};

		bool hasHaploid = false;
		bool hasDiploid = false;

		for (const auto &d : data) {
			if (!d.depth) continue;
			if (d.isHaploid()) {
				hasHaploid = true;
				glfs.emplace_back(TStrongArray<coretools::Probability, Base>{{(coretools::Probability)d[Base::A], (coretools::Probability)d[Base::G],
																   (coretools::Probability)d[Base::C], (coretools::Probability)d[Base::T]}});
				for (auto ac : usedAllelicCombinations) {
					Ls[ac].add(0.5 * (glfs.back()[genometools::first(ac)] + glfs.back()[genometools::second(ac)]));
				}
			} else {
				hasDiploid = true;
				glfs.emplace_back(TStrongArray<coretools::Probability, Genotype>{
					{(coretools::Probability)d[Genotype::AA], (coretools::Probability)d[Genotype::AC], (coretools::Probability)d[Genotype::AG],
					 (coretools::Probability)d[Genotype::AT], (coretools::Probability)d[Genotype::CC], (coretools::Probability)d[Genotype::CG],
					 (coretools::Probability)d[Genotype::CT], (coretools::Probability)d[Genotype::GG], (coretools::Probability)d[Genotype::GT],
					 (coretools::Probability)d[Genotype::TT]}});
				for (auto ac : usedAllelicCombinations) {
					Ls[ac].add(
						0.25 * (glfs.back()[genometools::homoFirst(ac)] + glfs.back()[genometools::homoSecond(ac)]) +
						0.5 * glfs.back()[genometools::het(ac)]);
				}
			}
		}

		TStrongArray<double, genometools::AllelicCombination> LLs{std::numeric_limits<double>::lowest()};
		for (auto ac : usedAllelicCombinations) { LLs[ac] = Ls[ac].getSum(); }
		const genometools::AllelicCombination bestAC = chooseBestAllelicCombination(LLs);

		auto [MAF, bestL, major, minor] = [bestAC, maxF, hasHaploid, hasDiploid, minMAF]() {
			if (hasHaploid) {
				if (hasDiploid) {
					return _iterate<true, true>(glfs, bestAC, maxF, minMAF);
				} else {
					return _iterate<true, false>(glfs, bestAC, maxF, minMAF);
				}
			} else {
				if (hasDiploid) {
					return _iterate<false, true>(glfs, bestAC, maxF, minMAF);
				} else {
					throw coretools::TDevError("No Data!");
				}
			}
		}();

		// determine variant quality
		Log10Probability LL_fixed = LLFixedAllele(data, major);
		const coretools::PhredInt variantQuality{LL_fixed > bestL ? Log10Probability(0.0)
																			   : Log10Probability(LL_fixed - bestL)};

		if (MAF < minMAF || variantQuality < minVariantQuality) {
			return failedTMMData;
		}

		// ensure return value uses major = ref if ref was provided
		if (minor == ref) { // cannot happen if base == ref
			minor = major;
			major = ref;
		}

		return {true, MAF, variantQuality, major, minor};
	}
public:

	static TMMData estimate(coretools::TConstView<genometools::TGLFEntry> data, double maxF,
							coretools::Probability minMAF, coretools::PhredInt minVariantQuality) {
		using AC = genometools::AllelicCombination;
		constexpr std::array usedAllelicCombinations = {AC::AC, AC::AG, AC::AT, AC::CG, AC::CT, AC::GT};
		return _estimate(data, maxF, usedAllelicCombinations, genometools::Base::N, minMAF, minVariantQuality);
	}

	static TMMData estimate(coretools::TConstView<genometools::TGLFEntry> data, double maxF,
							genometools::Base ref, coretools::Probability minMAF,
							coretools::PhredInt minVariantQuality) {
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

}

#endif
