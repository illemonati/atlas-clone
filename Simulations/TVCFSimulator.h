//
// Created by caduffm on 4/5/22.
//

#ifndef ATLAS_TVCFSIMULATOR_H
#define ATLAS_TVCFSIMULATOR_H

#include "TGenotypeFrequencies.h"
#include "TGlfMultiReader.h"
#include "TSimulator.h"
#include "probability.h"

namespace Simulations {

class TVCFWriterSimulation : public GLF::TGlfMultiReaderVcf {
private:
	void _writeHaploidIndividualToVCF(const GLF::TMultiGLFDataSampleOneAllelicCombination &GenotypeLikelihoods, size_t Depth);
	void _writeDiploidIndividualToVCF(const GLF::TMultiGLFDataSampleOneAllelicCombination &GenotypeLikelihoods, size_t Depth);

public:
	TVCFWriterSimulation(const std::string &Filename, const std::string &Source,
	                     const std::vector<std::string> &SampleNames, bool UsePhredScaledLikelihoods);

	void writeSite(const std::string &ChrName, uint32_t Position,
	          genometools::PhredIntProbability VariantQuality,
	          genometools::Base RefAllele, genometools::Base AltAllele, size_t TotalDepth,
	          bool IsDiploid,
	          const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods,
	          const std::vector<size_t> & Depths);
};

class TVCFSimulator : TSimulator {
private:
	coretools::Probability _error = 0.05;
	std::unique_ptr<TVCFWriterSimulation> _vcf;

	GLF::TMultiGLFDataSampleOneAllelicCombination _calculateGenotypeLikelihoods(size_t NumRef, size_t NumAlt,
	                                                                            bool IsDiploid);
	size_t _simulateNumReadsWithReferenceAllele(genometools::Base a, genometools::Base b, genometools::Base Ref,
	                                            size_t Depth);
	auto _simulateDepthAndGTL(genometools::Base a, genometools::Base b, genometools::Base Ref, bool IsDiploid);
	auto _calculateLog10ProbThisAllelicCombination(const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods);
	auto _calculateLog10ProbFixed(const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods,
	                              genometools::Base Major, genometools::Base Ref);
	auto _calculateVariantQuality(const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods,
	                              genometools::Base Major, genometools::Base Ref);
	auto _findMajorMinorAllele(coretools::TStrongArray<size_t, genometools::Base, 4> AlleleCounts);

public:
	TVCFSimulator(const std::string& method);

	// simulate reads and write bam files
	void _simulateAndWrite(const genometools::TChromosome & Chromosome, TSimulatorHaplotypes & Haplotypes) override;
};

} // namespace Simulations

#endif // ATLAS_TVCFSIMULATOR_H
