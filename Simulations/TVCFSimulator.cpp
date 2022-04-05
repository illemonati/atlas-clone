//
// Created by caduffm on 4/5/22.
//

#include "TVCFSimulator.h"

namespace Simulations {
using coretools::Probability;
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using coretools::str::toString;
using genometools::HighPrecisionPhredIntProbability;

//-------------------------------------------
// TVCFSimulationWriter
//-------------------------------------------

TVCFWriterSimulation::TVCFWriterSimulation(const std::string &Filename, const std::string &Source,
                                           const std::vector<std::string> &SampleNames, bool UsePhredScaledLikelihoods)
    : GLF::TGlfMultiReaderVcf(Filename, Source, SampleNames, UsePhredScaledLikelihoods) {}

void TVCFWriterSimulation::writeSite(const std::string &ChrName, uint32_t Position,
                                     genometools::PhredIntProbability VariantQuality, genometools::Base RefAllele,
                                     genometools::Base AltAllele, size_t TotalDepth, bool IsDiploid,
                                     const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods,
                                     const std::vector<size_t> &Depths) {
	_setMajorMinor(RefAllele, AltAllele);
	_writeSiteInformation(ChrName, Position, VariantQuality, RefAllele, AltAllele, TotalDepth);

	for (size_t i = 0; i < GenotypeLikelihoods.size(); ++i) {
		if (IsDiploid) {
			_writeDiploidIndividualToVCF(GenotypeLikelihoods[i], Depths[i]);
		} else {
			_writeHaploidIndividualToVCF(GenotypeLikelihoods[i], Depths[i]);
		}
	}

	// end of line
	_vcf << '\n';
}

void TVCFWriterSimulation::_writeHaploidIndividualToVCF(
    const GLF::TMultiGLFDataSampleOneAllelicCombination &GenotypeLikelihoods, size_t Depth) {
	if (Depth == 0) { return _writeMissing(); }
	const auto LL_haploidFirst  = GenotypeLikelihoods[genometools::BiallelicGenotype::haploidFirst];
	const auto LL_haploidSecond = GenotypeLikelihoods[genometools::BiallelicGenotype::haploidSecond];

	// write genotype and genotype quality
	const auto minQual = _writeGenotypeAndQualityHaploid({LL_haploidFirst, LL_haploidSecond});

	// write depth
	_vcf << Depth << ':';

	// write likelihoods
	_writeLikelihood(LL_haploidFirst - minQual);
	_vcf << ',';
	_writeLikelihood(LL_haploidSecond - minQual);
}

void TVCFWriterSimulation::_writeDiploidIndividualToVCF(
    const GLF::TMultiGLFDataSampleOneAllelicCombination &GenotypeLikelihoods, size_t Depth) {
	if (Depth == 0) { return _writeMissing(); }
	const auto LL_homoFirst  = GenotypeLikelihoods[genometools::BiallelicGenotype::homoFirst];
	const auto LL_het        = GenotypeLikelihoods[genometools::BiallelicGenotype::het];
	const auto LL_homoSecond = GenotypeLikelihoods[genometools::BiallelicGenotype::homoSecond];

	// write genotype and genotype quality
	const auto minQual = _writeGenotypeAndQualityDiploid({LL_homoFirst, LL_het, LL_homoSecond});

	// write depth
	_vcf << Depth << ':';

	// write likelihoods
	_writeLikelihood(LL_homoFirst - minQual);
	_vcf << ',';
	_writeLikelihood(LL_het - minQual);
	_vcf << ',';
	_writeLikelihood(LL_homoSecond - minQual);
}

//-------------------------------------------
// TVCFSimulator
//-------------------------------------------

TVCFSimulator::TVCFSimulator(const std::string &method) : TSimulator(method) {
	// read simulation parameters
	_error = parameters().getParameterWithDefault("error", 0.05);
	logfile().list("Will use a per allele genotyping error rate of " + coretools::str::toString(_error) + ".");

	const bool usePhredLikelihoods = parameters().parameterExists("phredLik");
	if (usePhredLikelihoods) {
		logfile().list("Will write phred-scaled likelihoods. (parameter phredLik)");
	} else {
		logfile().list("Will write raw likelihoods. (use phredLik to phred-scale)");
	}

	// generate sample names
	std::vector<std::string> sampleNames(_haploSimulator->sampleSize());
	for (size_t i = 0; i < _haploSimulator->sampleSize(); i++) { sampleNames[i] = "ind_" + coretools::str::toString(i); }

	// open VCF file and write header
	std::string filename = _outname + ".vcf.gz";
	logfile().listFlush("Writing vcf to file " + filename + "...");
	_vcf = std::make_unique<TVCFWriterSimulation>(filename, "ATLAS_simulations", sampleNames, usePhredLikelihoods);
}

GLF::TMultiGLFDataSampleOneAllelicCombination TVCFSimulator::_calculateGenotypeLikelihoods(size_t NumRef, size_t NumAlt,
                                                                                           bool IsDiploid) {
	const auto P0 =
	    HighPrecisionPhredIntProbability(Probability(pow(1. - _error, NumRef) * pow(_error, (double)NumAlt)));
	const auto P2 =
	    HighPrecisionPhredIntProbability(Probability(pow(_error, (double)NumRef) * pow(1. - _error, (double)NumAlt)));
	if (IsDiploid) {
		const auto P1 = HighPrecisionPhredIntProbability(Probability(pow(0.5, NumRef + NumAlt)));
		return {P0, P1, P2};
	} else {
		return {P0, P2};
	}
}

size_t TVCFSimulator::_simulateNumReadsWithReferenceAllele(genometools::Base a, genometools::Base b,
                                                           genometools::Base Ref, size_t Depth) {
	if (a == b && a == Ref) { // homozygous reference (haploid and diploid)
		return randomGenerator().getBinomialRand(1. - _error, Depth);
	} else if (a == b && a != Ref) { // homozygous alternative (haploid and diploid)
		return randomGenerator().getBinomialRand(_error, Depth);
	}
	return randomGenerator().getBinomialRand(0.5, Depth); // heterozygous
}

auto TVCFSimulator::_simulateDepthAndGTL(genometools::Base a, genometools::Base b, genometools::Base Ref,
                                         bool IsDiploid) {
	// simulate depth
	auto depth = randomGenerator().getPoissonRandom(_seqDepth);

	// simulate number of reference and alternative alleles
	auto numRef = _simulateNumReadsWithReferenceAllele(a, b, Ref, depth);
	auto numAlt = (int)depth - numRef;

	// get genotype likelihoods
	auto genotypeLikelihoods = _calculateGenotypeLikelihoods(numRef, numAlt, IsDiploid);

	return std::make_pair(depth, genotypeLikelihoods);
}

auto TVCFSimulator::_calculateLog10ProbThisAllelicCombination(
    const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods) {
	// estimate genotype frequencies
	genometools::TGenotypeFrequencies genotypeFrequencies;
	genotypeFrequencies.estimate(GenotypeLikelihoods, GenotypeLikelihoods.size(), 0.0000001);

	// calculate log10 likelihood
	return genotypeFrequencies.calculateLog10Likelihood(GenotypeLikelihoods, GenotypeLikelihoods.size());
}

auto TVCFSimulator::_calculateLog10ProbFixed(const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods,
                                             genometools::Base Major, genometools::Base Ref) {
	using coretools::Log10Probability;

	const bool majorIsRef = (Major == Ref);
	const auto ref =
	    majorIsRef ? genometools::BiallelicGenotype::haploidFirst : genometools::BiallelicGenotype::haploidSecond;
	const auto refHomo =
	    majorIsRef ? genometools::BiallelicGenotype::homoFirst : genometools::BiallelicGenotype::homoSecond;

	Log10Probability LL_fixed_glfPhred = 0.0;
	for (auto GenotypeLikelihood : GenotypeLikelihoods) {
		if (GenotypeLikelihood.isHaploid()) {
			LL_fixed_glfPhred += (Log10Probability)GenotypeLikelihood[ref];
		} else {
			LL_fixed_glfPhred += (Log10Probability)GenotypeLikelihood[refHomo];
		}
	}
	return LL_fixed_glfPhred;
}

auto TVCFSimulator::_calculateVariantQuality(const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods,
                                             genometools::Base Major, genometools::Base Ref) {
	const auto LL_allelicCombination = _calculateLog10ProbThisAllelicCombination(GenotypeLikelihoods);
	const auto LL_fixed              = _calculateLog10ProbFixed(GenotypeLikelihoods, Major, Ref);

	// calculate variant quality
	return genometools::PhredIntProbability(LL_fixed > LL_allelicCombination ? coretools::Log10Probability(0.0)
	                                                                         : LL_fixed - LL_allelicCombination);
}

auto TVCFSimulator::_findMajorMinorAllele(coretools::TStrongArray<size_t, genometools::Base, 4> AlleleCounts) {
	// first find allele with highest counts
	const auto max             = std::max_element(AlleleCounts.begin(), AlleleCounts.end());
	const auto maxAllele       = (genometools::Base)std::distance(AlleleCounts.begin(), max);
	// set it to zero to find allele with second-highest counts
	*max                       = 0;
	const auto secondMax       = std::max_element(AlleleCounts.begin(), AlleleCounts.end());
	const auto secondMaxAllele = (genometools::Base)std::distance(AlleleCounts.begin(), secondMax);

	return std::make_pair(maxAllele, secondMaxAllele);
}

void TVCFSimulator::_simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes) {
	for (size_t l = 0; l < Chromosome.length; ++l) {

		// prepare storage
		GLF::TMultiGLFDataOneAllelicCombination genotypeLikelihoods(_haploSimulator->sampleSize());
		std::vector<size_t> depths(_haploSimulator->sampleSize(), 0);
		coretools::TStrongArray<size_t, genometools::Base, 4> alleleCounts({0, 0, 0, 0});
		const bool isDiploid = Chromosome.ploidy == 2;

		for (size_t i = 0; i < _haploSimulator->sampleSize(); ++i) {
			// get haplotypes
			const auto hap1 = Haplotypes((int)i, 0, l);
			const auto hap2 = Haplotypes((int)i, 1, l);

			// simulate depth and genotype likelihoods
			auto [depth, GTL] = _simulateDepthAndGTL(hap1, hap2, _reference[l], isDiploid);

			// store
			genotypeLikelihoods[i] = GTL;
			depths[i]              = depth;
			alleleCounts[hap1]++;
			if (isDiploid) { alleleCounts[hap2]++; }
		}

		// find major and minor allele
		const auto [majorAllele, minorAllele] = _findMajorMinorAllele(alleleCounts);
		const auto refAllele                  = _reference[l];
		if (refAllele != majorAllele && refAllele != minorAllele) {
			throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": reference is not major nor minor allele!");
		}
		const auto altAllele = (refAllele == majorAllele) ? minorAllele : majorAllele;

		// calculate variant quality
		auto variantQuality = _calculateVariantQuality(genotypeLikelihoods, majorAllele, refAllele);

		// finally write site!
		_vcf->writeSite(Chromosome.name, l + 1, variantQuality, refAllele, altAllele, coretools::containerSum(depths),
		                isDiploid, genotypeLikelihoods, depths);
	}
	logfile().done();
}

} // namespace Simulations