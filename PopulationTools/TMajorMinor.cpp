/*
 * TMajorMinor.cpp
 *
 *  Created on: Nov 5, 2018
 *      Author: phaentu
 */

#include "TMajorMinor.h"

#include <algorithm>
#include <array>
#include <exception>
#include <iterator>
#include <memory>
#include <numeric>
#include <set>
#include <stddef.h>
#include <stdint.h>
#include <string>
#include <tuple>
#include <vector>

#include "TGlfMultiReader.h"
#include "coretools/Containers/TDualArray.h"
#include "coretools/Containers/TDualStrongArray.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Main/TTask.h"
#include "coretools/Math/TSumLog.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/TTimer.h"
#include "coretools/Types/probability.h"
#include "coretools/Types/strongTypes.h"
#include "coretools/Types/weakTypes.h"
#include "coretools/devtools.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/TGenotypeFrequencies.h"

#ifdef _OPENMP
#include "omp.h"
#endif


namespace PopulationTools {


using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using genometools::AllelicCombination;
using genometools::Base;
using genometools::Genotype;
using coretools::index;
using coretools::Probability;
using coretools::TConstView;
using coretools::Log10Probability;
using coretools::TStrongArray;
using coretools::TDualStrongArray;
using GLF::TMultiGLFData;
//---------------------------------------------------
// TMajorMinorEstimatorBase
//---------------------------------------------------

namespace impl {

constexpr coretools::TDualArray<AllelicCombination, 3, index(AllelicCombination::max)>
useAllelicCombinationsThatContain(Base base) {
	using AC = AllelicCombination;
	switch (base) {
	case Base::A: return std::array{AC::AC, AC::AG, AC::AT};
	case Base::C: return std::array{AC::AC, AC::CG, AC::CT};
	case Base::G: return std::array{AC::AG, AC::CG, AC::GT};
	case Base::T: return std::array{AC::AT, AC::CT, AC::GT};
	case Base::N: return std::array{AC::AC, AC::AG, AC::AT, AC::CG, AC::CT, AC::GT};
	}
};

template <typename Container>
AllelicCombination chooseBestAllelicCombination(const Container &acd) {
	return randomGenerator()
		.sampleIndexOfMaxima<Container, AllelicCombination>(acd);
};

} // namespace impl

struct TMMData {
	bool pass;
	Probability MAF;
	genometools::PhredIntProbability variantQuality;
	genometools::Base major;
	genometools::Base minor;
};

class TSkotte {
	enum class HaploDiplo : size_t { min, first = min, second, homoFirst, het, homoSecond, max };

	static HaploDiplo _haploIndex(Probability pFirst, Probability pSecond) {
		return HaploDiplo(pSecond.isHigher(pFirst));
	}

	static HaploDiplo _diploIndex(Probability pHomoFirst, Probability pHet, Probability pHomoSecond) {
		std::array<Probability, 3> es = {pHomoFirst, pHet, pHomoSecond};
		return HaploDiplo(2 + std::distance(es.begin(), std::max_element(es.begin(), es.end())));
	}

	static double _haploWeights(Probability pFirst, Probability pSecond,
							  const TStrongArray<Probability, HaploDiplo> freqs) {
		const double weights0 = pFirst * freqs[HaploDiplo::first];
		const double weights1 = pSecond * freqs[HaploDiplo::second];
		const double sum      = weights0 + weights1;
		return weights0 / sum;
	}

	static std::pair<double, double> _diploWeights(Probability pHomoFirst, Probability pHet, Probability pHomoSecond,
												   const TStrongArray<Probability, HaploDiplo> freqs) {
		const double weights0 = pHomoFirst * freqs[HaploDiplo::homoFirst];
		const double weights1 = pHet * freqs[HaploDiplo::het];
		const double weights2 = pHomoSecond * freqs[HaploDiplo::homoSecond];
		const double sum_i    = 1./(weights0 + weights1 + weights2);
		return {weights0 * sum_i, weights2 * sum_i};
	}

	template<bool hasHaploid, bool hasDiploid>
	static auto _iterate(TConstView<coretools::TDualStrongArray<Probability, Base, Genotype>> glfs, AllelicCombination ac, double maxF, Probability minMAF) {
		constexpr bool hasBoth = hasHaploid && hasDiploid;

		const auto first      = genometools::first(ac);
		const auto second     = genometools::second(ac);
		const auto homoFirst  = genometools::homoFirst(ac);
		const auto het        = genometools::het(ac);
		const auto homoSecond = genometools::homoSecond(ac);

		// Initial guess
		TStrongArray<size_t, HaploDiplo> counts{0};
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

		constexpr Probability fMin = 0.000001;
		TStrongArray<Probability, HaploDiplo> freqs{0};
		const auto nHaplo     = counts[HaploDiplo::first] + counts[HaploDiplo::second];
		const double nHaplo_1 = 1. / nHaplo;
		if constexpr (hasHaploid) {
			freqs[HaploDiplo::first]  = std::max<Probability>(double(counts[HaploDiplo::first]) * nHaplo_1, fMin);
			freqs[HaploDiplo::second] = std::max<Probability>(double(counts[HaploDiplo::second]) * nHaplo_1, fMin);
			const auto sum            = freqs[HaploDiplo::first] + freqs[HaploDiplo::second];
			freqs[HaploDiplo::first].scale(sum);
			freqs[HaploDiplo::second].scale(sum);
		}
		const auto nDiplo     = counts[HaploDiplo::homoFirst] + counts[HaploDiplo::het] + counts[HaploDiplo::homoSecond];
		const double nDiplo_1 = 1. / nDiplo;
		if constexpr (hasDiploid) {
			freqs[HaploDiplo::homoFirst]  = std::max<Probability>(double(counts[HaploDiplo::homoFirst]) * nDiplo_1, fMin);
			freqs[HaploDiplo::het]        = std::max<Probability>(double(counts[HaploDiplo::het]) * nDiplo_1, fMin);
			freqs[HaploDiplo::homoSecond] = std::max<Probability>(double(counts[HaploDiplo::homoSecond]) * nDiplo_1, fMin);
			const auto sum = freqs[HaploDiplo::homoFirst] + freqs[HaploDiplo::het] + freqs[HaploDiplo::homoSecond];
			freqs[HaploDiplo::homoFirst].scale(sum);
			freqs[HaploDiplo::het].scale(sum);
			freqs[HaploDiplo::homoSecond].scale(sum);
		}

		const double nHaplo2nDiplo_1 = 1. / (nHaplo + 2*nDiplo);

		// iterate
		constexpr size_t maxIter = 1000;

		Probability aF = (nDiplo * (freqs[HaploDiplo::het] + 2.0 * freqs[HaploDiplo::homoSecond]) +
						  nHaplo * (double)freqs[HaploDiplo::second]) *
						 nHaplo2nDiplo_1;
		Probability MAF = (aF < 0.5 ? aF : aF.complement());

		if (MAF < minMAF/2) {
			// will be automaticall filtered out
			return std::make_tuple(Probability::lowest(), Log10Probability::lowest(), first, second);
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

				freqs[HaploDiplo::first] = hplF0;
				freqs[HaploDiplo::second]= hplF1;
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

				freqs[HaploDiplo::homoFirst]  = dplF0;
				freqs[HaploDiplo::het]        = dplF1;
				freqs[HaploDiplo::homoSecond] = dplF2;
			}

			aF = (nDiplo * (freqs[HaploDiplo::het] + 2.0 * freqs[HaploDiplo::homoSecond]) +
				  nHaplo * (double)freqs[HaploDiplo::second]) *
				 nHaplo2nDiplo_1;

			MAF = aF < 0.5 ? aF : aF.complement();

			if (MAF < minMAF - maxF_i) {
				// will be automaticall filtered out
				return std::make_tuple(Probability::lowest(), Log10Probability::lowest(), first, second);
			}

			if (maxF_i < maxF)  {
				break;
			}
		}

		const auto [major, minor] = [aF, first ,second]() {
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
		constexpr double l10_1       = 0.43429448190325176;
		const Log10Probability bestL = L.getSum()*l10_1;

		return std::make_tuple(MAF, bestL, major, minor);
	}

public:
	static TMMData estimate(coretools::TConstView<GLF::TMultiGLFDataSample> data, double maxF, genometools::Base base, Probability minMAF, genometools::PhredIntProbability minVariantQuality) {
		std::vector<coretools::TDualStrongArray<Probability, Base, Genotype>> glfs;
		glfs.reserve(data.size());

		const auto used = impl::useAllelicCombinationsThatContain(base);
		coretools::TStrongArray<coretools::TSumLogProbability, AllelicCombination> Ls{};

		bool hasHaploid = false;
		bool hasDiploid = false;

		for (const auto &d : data) {
			if (!d.hasData()) continue;
			if (d.isHaploid()) {
				hasHaploid = true;
				glfs.emplace_back(TStrongArray<Probability, Base>{{(Probability)d[Base::A], (Probability)d[Base::G],
																   (Probability)d[Base::C], (Probability)d[Base::T]}});
				for (auto ac : used) {
					Ls[ac].add(0.5 * (glfs.back()[genometools::first(ac)] + glfs.back()[genometools::second(ac)]));
				}
			} else {
				hasDiploid = true;
				glfs.emplace_back(TStrongArray<Probability, Genotype>{
				    {(Probability)d[Genotype::AA], (Probability)d[Genotype::AC],
				     (Probability)d[Genotype::AG], (Probability)d[Genotype::AT],
				     (Probability)d[Genotype::CC], (Probability)d[Genotype::CG],
				     (Probability)d[Genotype::CT], (Probability)d[Genotype::GG],
				     (Probability)d[Genotype::GT], (Probability)d[Genotype::TT]}});
				for (auto ac : used) {
					Ls[ac].add(0.25 * (glfs.back()[genometools::homoFirst(ac)] + glfs.back()[genometools::homoSecond(ac)]) +
							   0.5 * glfs.back()[genometools::het(ac)]);
				}
			}
		}

		TStrongArray<double, AllelicCombination> LLs{std::numeric_limits<double>::lowest()};
		for (auto ac : used) {
			LLs[ac] = Ls[ac].getSum();
		}
		const auto bestAC     = impl::chooseBestAllelicCombination(LLs);

		const auto [MAF, bestL, major, minor] = [&glfs, bestAC, maxF,hasHaploid, hasDiploid, minMAF]() {
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
				DEVERROR("No Data!");
			}
		}}();

		if (MAF < minMAF) return {false, Probability::lowest(), genometools::PhredIntProbability::lowest(), major, minor};

		const auto refHom         = genometools::genotype(major, major);
		Log10Probability LL_fixed = 0.0;
		for (size_t i = 0; i < data.size(); ++i) {
			if (data[i].hasData()) {
				if (data[i].isHaploid())
					LL_fixed += (Log10Probability)data[i][major];
				else
					LL_fixed += (Log10Probability)data[i][refHom];
			}
		}
		const genometools::PhredIntProbability variantQuality{
			LL_fixed > bestL ? Log10Probability(0.0) : Log10Probability(LL_fixed - bestL)};

		if (variantQuality < minVariantQuality) {
			return {false, Probability::lowest(), genometools::PhredIntProbability::lowest(), major, minor};
		}

		return {true, MAF, variantQuality, major, minor};
	}
};

struct TMLE {
	static TMMData estimate(coretools::TConstView<GLF::TMultiGLFDataSample> data, double maxF, genometools::Base base, Probability minMAF, genometools::PhredIntProbability minVariantQuality) {
		// calculate L10L for each allelic combination
		const auto used = impl::useAllelicCombinationsThatContain(base);

		genometools::TGenotypeFrequencies bestFreqs;
		coretools::Log10Probability bestL = coretools::Log10Probability::lowest();
		AllelicCombination bestAC         = AllelicCombination::min;

		for (const auto ac : used) {
			auto GLs = fill(data, ac);
			genometools::TGenotypeFrequencies freqs;
			freqs.estimate<false>(GLs, GLs.size(), maxF);
			auto L = freqs.calculateLog10Likelihood(GLs, GLs.size());
			if ((L > bestL) || (L == bestL && randomGenerator().getRand() > 0.5)) {
				bestL     = L;
				bestFreqs = freqs;
				bestAC    = ac;
			}
		}
		Base _major, _minor;

		// which one is major?
		if (bestFreqs.alleleFrequency() < 0.5) {
			_major = first(bestAC);
			_minor = second(bestAC);
		} else {
			_major = second(bestAC);
			_minor = first(bestAC);
		}

		if (_minor == base) { // cannot happen if base == ref
			_minor = _major;
			_major = base;
		}

		// calculate variant quality
		const auto refHom                  = genometools::genotype(_major, _major);
		Log10Probability LL_fixed_glfPhred = 0.0;
		for (size_t i = 0; i < data.size(); ++i) {
			if (data[i].hasData()) {
				if (data[i].isHaploid())
					LL_fixed_glfPhred += (Log10Probability)data[i][_major];
				else
					LL_fixed_glfPhred += (Log10Probability)data[i][refHom];
			}
		}
		genometools::PhredIntProbability _variantQuality{
			LL_fixed_glfPhred > bestL ? Log10Probability(0.0) : Log10Probability(LL_fixed_glfPhred - bestL)};

		if (bestFreqs.MAF() < minMAF || _variantQuality < minVariantQuality) {
			return {false, Probability::lowest(), genometools::PhredIntProbability::lowest(), _major, _minor};
		}
		return {true, bestFreqs.MAF(), _variantQuality, _major, _minor};
	}
};

template<typename Estimator> void iterate(double maxF) {
	// open GLF files
	GLF::TGlfMultiReader glfReader;
	glfReader.openGLFs();
	glfReader.setAllActive();

	// add reference, if provided
	if (parameters().parameterExists("fasta")) {
		logfile().list("Will only identify the most likely alternative allele (argument: fasta)");
		const std::string fastaFile = parameters().getParameter<std::string>("fasta");
		logfile().list("Reading reference sequence from '" + fastaFile + "'");
		glfReader.addReference(fastaFile);
	}

	// estimation method
	const std::string method = parameters().getParameterWithDefault<std::string>("method", "MLE");

	const bool usePhredLikelihoods = parameters().parameterExists("phredLik");
	if (usePhredLikelihoods) {
		logfile().list("Will write phred-scaled likelihoods. (parameter phredLik)");
	} else {
		logfile().list("Will write raw likelihoods. (use phredLik to phred-scale)");
	}

	// read filters
	size_t minSamplesWithData = 1;
	genometools::PhredIntProbability minVariantQuality{0};
	if (parameters().parameterExists("printAll")) {
		logfile().list("Will all sites and samples. (parameter printAll)");
		minSamplesWithData = 0;
		minVariantQuality  = 0;
	} else {
		minSamplesWithData = parameters().getParameterWithDefault<size_t>("minSamplesWithData", 1);
		if (minSamplesWithData > 0) {
			logfile().list("Will only print sites for which at least ", minSamplesWithData,
						   " samples have data. (parameter minSamplesWithData)");
		}

		minVariantQuality = parameters().getParameterWithDefault<genometools::PhredIntProbability>(
			"minVariantQual", genometools::PhredIntProbability::highest());
		if (minVariantQuality > genometools::PhredIntProbability::highest()) {
			logfile().list("Will only print sites with variant quality >= ", minVariantQuality,
						   " samples have data. (parameter minVariantQual)");
		}
	}
	glfReader.minSamplesWithData(minSamplesWithData);

	coretools::Probability minMAF = parameters().getParameterWithDefault("minMAF", 0.0);
	if (minMAF > 0.0) {
		logfile().list("Will filter on a minor allele frequency of ", minMAF, ". (parameter 'minMAF')");
	} else {
		logfile().list("Will keep sites regardless of their minor allele frequency. (use 'minMAF' to filter)");
	}

	// limit input
	const size_t limitSites = parameters().getParameterWithDefault("limitSites", 0);
	if (limitSites > 0) logfile().list("Will stop at input position ", limitSites, ". (parameter 'limitSites')");
	if (limitSites < 0) UERROR("maxPos cannot be negative!");

	// filename tag
	const std::string outname = parameters().getParameterWithDefault<std::string>("out", "ATLAS_majorMinor");
	logfile().list("Will write output files with tag '" + outname + "'. (parameter 'out')");

	// get sample names. Make sure they are unique in the vcf
	std::vector<std::string> sampleNames;
	if (parameters().parameterExists("sampleNames")) {
		logfile().startIndent("Using the following sample names (parameter 'sampleNames'):");
		parameters().fillParameterIntoContainer("sampleNames", sampleNames, ',');

		if (sampleNames.size() != glfReader.numActiveSamples()) {
			UERROR("Number of provided same names does not match number of GLF files!");
		}

		// report
		auto glfNames = glfReader.namesOfActiveFiles();
		for (size_t i = 0; i < glfReader.numActiveSamples(); ++i) {
			logfile().list(glfNames[i], " -> ", sampleNames[i]);
		}
		logfile().endIndent();
	} else {
		logfile().list(
			"Will deduce sample names from GLF file names. (use 'sampleNames' to provide alternative names)");
		sampleNames = glfReader.sampleNamesOfActiveFiles();

		// if there are duplicates, add suffix
		bool foundDuplicates = false;
		for (size_t i = 0; i < sampleNames.size(); ++i) {
			for (size_t j = i + 1; j < sampleNames.size(); ++j) {
				int counter = 1;
				if (sampleNames[i] == sampleNames[j]) {
					sampleNames[j] += "." + coretools::str::toString(counter++);
					if (!foundDuplicates) {
						logfile().startIndent("Duplicate samples will be rename as follows:");
						foundDuplicates = true;
					}
					logfile().list(sampleNames[i], " -> ", sampleNames[j]);
				}
			}
		}
		if (foundDuplicates) { logfile().endIndent(); }
	}

#ifdef _OPENMP
	size_t maxThreads = coretools::instances::parameters().getParameterWithDefault("maxThreads", omp_get_max_threads());
	coretools::instances::logfile().list("Running in parallel with a maximum of ", maxThreads,
										 " threads (argument 'maxThreads')");
#else
	coretools::instances::logfile().list("Not running in parallel");
#endif

	// open vcf file
	GLF::TGlfMultiReaderVcf vcf(outname + ".vcf.gz", "ATLAS_GLF_Caller", sampleNames, usePhredLikelihoods);

	// vars
	logfile().startIndent("Parsing through glf files:");
	coretools::TTimer timer;
	constexpr size_t dCounter = 1000000;
	size_t counter            = 0;
	size_t nextPrint          = dCounter;

	for (auto ids = glfReader.readWindow(); !ids.empty(); ids = glfReader.readWindow()) {
		std::vector<TMMData> data(ids.size());

#pragma omp parallel for num_threads(maxThreads)
		for (size_t i = 0; i < ids.size(); ++i) {
			const auto iW  = ids[i];
			const Base ref = glfReader.refBase(iW); // can be N
			data[i]        = Estimator::estimate(glfReader.data(iW), maxF, ref, minMAF, minVariantQuality);
		}

		// pass filter?
		for (size_t i = 0; i < ids.size(); ++i) {
			const auto iW = ids[i];
			const auto &di = data[i];
			if (!di.pass) continue;

			// write to VCF
			vcf.writeSite(glfReader.chr(), glfReader.position(iW), di.variantQuality, glfReader.data(iW), di.major,
						  di.minor);
		}

		counter += ids.back() + 1;

		// report progress
		if (counter >= nextPrint) {
			logfile().list("Parsed ", nextPrint, " positions in ", timer.formattedTime(), ".");
			while (nextPrint <= counter) nextPrint += dCounter;
		}

		if (limitSites > 0 && counter == limitSites) break;
	}

	logfile().list("Reached end of glf files!");
	logfile().removeIndent();
};

//---------------------------------------------------
// TMajorMinor
//---------------------------------------------------
void TMajorMinor::run() {
	const std::string method = parameters().getParameterWithDefault<std::string>("method", "MLE");

	const double maxF = parameters().getParameterWithDefault("maxF", 0.0000001);
	if (method == "Skotte") {
		logfile().list("Will estimate major / minor alleles using the Skotte method with maxF ", maxF,
					   ". (parameters method and maxF)");
		iterate<TSkotte>(maxF);
	} else if (method == "MLE") {
		logfile().list("Will estimate major / minor alleles using the MLE method with maxF ", maxF,
					   ". (parameters method and maxF)");
		iterate<TMLE>(maxF);
	} else {
		UERROR("Unknown MajorMinor method '", method, "'!");
	}
}

}; // namespace PopulationTools
