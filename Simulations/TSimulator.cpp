/*
 * TSimulator.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#include "TSimulator.h"

#include <math.h>
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <memory>
#include <numeric>
#include <ostream>
#include <stdexcept>
#include <type_traits>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include "coretools/Files/TFile.h"
#include "genometools/TGenotypeFrequencies.h"
#include "THaplotypeSimulator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/algorithms.h"
#include "coretools/Files/gzstream.h"
#include "coretools/Math/mathConstants.h"
#include "coretools/Types/probability.h"
#include "coretools/Main/progressTools.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/weakTypes.h"

namespace Simulations {
using coretools::Probability;
using coretools::LogProbability;
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using coretools::str::toString;
using genometools::Base;
using genometools::HighPrecisionPhredIntProbability;
using genometools::TChromosomes;
using genometools::TChromosome;
using genometools::TGenomePosition;
using BAM::RGInfo::TReadGroupInfo;
using BAM::RGInfo::TReadGroupInfoEntry;

//---------------------------------------------------------
// Helper functions
//---------------------------------------------------------

namespace /* anonymous */ {
std::unique_ptr<THaplotypeSimulator> makeHaploSimulator(const std::string &method, const TChromosomes &chs) {
	if (method == TSimulatorOne::name) {
		logfile().startIndent("Simulating a single individual (parameter 'type' = '", method, "')");
		return std::make_unique<TSimulatorOne>(chs.size());
	}
	if (method == TSimulatorPair::name) {
		logfile().startIndent("Simulating a pair of individual (parameter 'type' = '", method, "')");
		return std::make_unique<TSimulatorPair>();
	}
	if (method == TSimulatorSFS::name) {
		logfile().startIndent("Simulating individuals from an SFS (parameter 'type' = '", method, "')");
		return std::make_unique<TSimulatorSFS>(chs);
	}
	if (method == TSimulatorHW::name) {
		logfile().startIndent("Simulating individuals under Hardy-Weinberg (parameter 'type' = '", method, "')");
		return std::make_unique<TSimulatorHW>();
	}
	throw "Unknown simulation method '" + method + "'! Use '" + TSimulatorOne::name + "', '" + TSimulatorPair::name + "', '" + TSimulatorSFS::name + "' or '" + TSimulatorHW::name + "'";
	logfile().endIndent();
}

template <typename T>
std::vector<T> parse(const std::string & Argument, const std::string & Explanation, const T Default){
	std::vector<T> res;
	if(parameters().parameterExists(Argument)){
		std::vector<std::string> string_vec;
		parameters().fillParameterIntoContainer(Argument, string_vec, ',');

		coretools::str::repeatIndexes(string_vec, res);
		if(res.empty()) throw "Issue understanding " + Explanation + " '" + coretools::str::concatenateString(string_vec, ",") + "' (parameter '" + Argument + "')!";

		if(res.size() == 1){
			logfile().list("Will use ", Explanation, " ", res[0], ". (parameter '", Argument, "')");
		} else {
			//logfile().list("Will use ", Explanation, "s ", coretools::str::concatenateString(res, ","), ". (parameter '", Argument, "')");
			logfile().list("Will use ", Explanation, "s ", res, ". (parameter '", Argument, "')");
		}

	} else {
		res.push_back(Default); //default ploidy = 2
		logfile().list("Will use default chromosome length of ", res[0], ". (set with '", Argument, "')");
	}
	return res;
}

template <typename Vec>
void checkLength(Vec & vec, size_t numChr){
	if(vec.size() != numChr){
		if(vec.size() == 1){
			vec.resize(numChr, vec.front());
		} else {
			throw "Number of chromosomes implied by chrLength, ploidy and depth does not match!";
		}
	}
}

void makeChromosomes(TChromosomes & chs, std::vector<uint32_t> & depths){
	//TODO: make it possible to initialize chromosomes from a file with length, ploidy and depth
	logfile().startIndent("Parameters regarding chromosomes:");
	chs.clear();
	depths.clear();

	//parse chr details
	auto chrLengths = parse<uint32_t>("chrLength", "chromosome length", 1000000);
	depths = parse<uint32_t>("depth", "sequencing depth", 10);
	auto ploidies = parse<uint8_t>("ploidy", "ploidy", 2);

	//check if ploidy is supported
	for (auto &p : ploidies) {
		if (p != 1 && p != 2) { throw "Currently only ploidy 1 (haploid) or 2 (diploid) is supported!"; }
	}

	//check length
	size_t numChr = std::max( {chrLengths.size(), ploidies.size(), depths.size()} );
	checkLength(chrLengths, numChr);
	checkLength(ploidies, numChr);
	checkLength(depths, numChr);

	//report and create
	logfile().startIndent("List of ", chrLengths.size(), " chromosome(s) to simulate:");
	for(size_t i = 0; i < chrLengths.size(); ++i){
		//create chromosome
		const TChromosome& chr = chs.appendChromosome("chr" + coretools::str::toString(i + 1), chrLengths[i], ploidies[i]);

		std::string text = chr.name + " (";
		if (ploidies[i] == 1){
			text += "haploid) ";
		} else {
			text += "diploid) ";
		}
		text += "of length " + coretools::str::toString(chrLengths[i]) + " and average depth " + coretools::str::toString(depths[i]) + ".";
		logfile().list(text);
	}
	logfile().endIndent();
	logfile().endIndent();
}

} // namespace

//---------------------------------------------------------
// TSimulator
//---------------------------------------------------------

TSimulator::TSimulator(const std::string &method){
	// output settings
	logfile().startIndent("Output settings:");
	if(parameters().parameterExists("out")){
		_outname = parameters().getParameter<std::string>("out");
		logfile().list("Will write output files with tag '" + _outname + "'. (parameter 'out')");
	} else {
		_outname = "ATLAS_simulations";
		logfile().list("Will write output files with tag '" + _outname + "'. (set with 'out')");
	}

	_writeTrueGenotypes = parameters().parameterExists("writeTrueGenotypes");
	if(_writeTrueGenotypes){
		logfile().list("Will write true genotypes to file. (parameter 'writeTrueGenotypes')");
	} else {
		logfile().list("Will NOT write true genotypes to file. (request with 'writeTrueGenotypes')");
	}

	_writeVariantInvariantBedFiles = parameters().parameterExists("writeVariantBED");
	if(_writeVariantInvariantBedFiles){
		logfile().list("Will write BED files with variant and invariant positions. (parameter 'writeVariantBED')");
	} else {
		logfile().list("Will NOT write BED files with variant and invariant positions. (request with 'writeVariantBED')");
	}

	_reference.open(_outname + ".fasta");

	logfile().endIndent();

	//parse sequencing depth
	makeChromosomes(_chromosomes, _seqDepth);
	_haploSimulator = makeHaploSimulator(method, _chromosomes);
}

void TSimulator::runSimulations() {
	// prepare haplotypes and
	TSimulatorHaplotypes haplotypes(_haploSimulator->sampleSize());

	// open files to store extra info on sites
	if (_writeTrueGenotypes) {
		// open file for true genotypes
		const auto filename = _outname + "_trueGenotypes.vcf.gz";
		haplotypes.openTrueGenotypeVCF(filename);
	}

	TSimulatorVariantInvariantBedFiles bedFiles;
	if (_writeVariantInvariantBedFiles) bedFiles.open(_outname);

	// simulate sequences
	for(size_t i = 0; i < _chromosomes.size(); ++i){
		auto &chr = _chromosomes[i];
		logfile().startIndent("Simulating chromosome " + chr.name + ":");

		// update reference storage and update haplotype lengths
		_reference.setChr(chr.name, chr.length);
		haplotypes.setLength(chr.length);

		// simulate genotypes
		logfile().listFlush("Simulating genotypes ...");
		if (chr.ploidy == 1)
			_haploSimulator->simulateHaploid(haplotypes, _reference, chr);
		else
			_haploSimulator->simulateDiploid(haplotypes, _reference, chr);
		logfile().done();

		// write true genotypes
		if (_writeTrueGenotypes) {
			logfile().listFlush("Writing true genotypes ...");
			haplotypes.writeTrueGenotypes(chr.name, _reference);
			logfile().done();
		}

		// write BED files
		if (_writeVariantInvariantBedFiles) bedFiles.write(haplotypes, chr.name);

		// write bam / vcf files!
		_simulateAndWrite(chr, haplotypes, _seqDepth[i]);

		// end of chromosome
		logfile().endIndent();
	}
	logfile().endIndent();
}

//---------------------------------------------------
// TBamSimulator
//---------------------------------------------------

TBAMSimulator::TBAMSimulator(const std::string &method) : TSimulator(method) {
	using genometools::Base;
	_initializeReadSimulator();

	// open bam files
	_bamFiles =
	    std::make_unique<TSimulatorBamFiles>(_haploSimulator->sampleSize(), _outname, _readSimulators, _chromosomes);
}

void TBAMSimulator::_initializeReadSimulator(){
	logfile().startIndent("Parameters regarding sequencing:");
	//read RGInfo files from command line
	std::vector<std::string> filenames;
	if(parameters().parameterExists(BAM::RGInfo::TReadGroupInfo::_RGInfoArgument)){
		std::vector<std::string> tmp;
		parameters().fillParameterIntoContainer(BAM::RGInfo::TReadGroupInfo::_RGInfoArgument, tmp, ',');
		coretools::str::repeatIndexes(tmp, filenames);
	} else {
		filenames.push_back("");
	}

	//create read simulators
	_readSimulators.reserve(filenames.size());
	if(filenames.size() == 1){

		if(_haploSimulator->sampleSize() > 1){
			logfile().startIndent("Using one set of sequencing parameters for all ", _haploSimulator->sampleSize(), " individuals:");
		}
		_readSimulators.emplace_back(filenames.front());
		if(_haploSimulator->sampleSize() > 1){
			logfile().endIndent();
		}
	} else {
		logfile().startNumbering("Using individual-specific sequencing parameters:");
		//check sizes matches sample size
		if(_haploSimulator->sampleSize() != filenames.size()){
			UERROR("Number of read group info files does not match sample size!");
		}

		for(auto& s : filenames){
			logfile().numberWithIndent("Sequencing parameters for individual 1");
			_readSimulators.emplace_back(s);
		}
		logfile().endNumbering();
		logfile().endIndent();
	}

	//check if read length match chr length
	for(auto& chr : _chromosomes){
		for(auto& rs : _readSimulators){
			if(rs.maxFragmentLength() > chr.length){
				throw "Length of chromosome '" + chr.name + "' is less than the max fragment length of some read groups!";
			}
		}
	}
};

void TBAMSimulator::_simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, uint32_t avgDepth) {
	// now simulate and write reads
	logfile().startIndent("Simulating reads:");
	for (size_t i = 0; i < _haploSimulator->sampleSize(); ++i){
		if(_readSimulators.size() == 1){
			_simulateReadsFromHaplotypes(Chromosome,
										 Haplotypes.getHaplotypesOfIndividual(i),
										 _readSimulators.front(),
										 avgDepth,
										 (*_bamFiles)[i],
										 " for individual " + coretools::str::toString(i + 1));
		} else {
			_simulateReadsFromHaplotypes(Chromosome,
										 Haplotypes.getHaplotypesOfIndividual(i),
										 _readSimulators[i],
										 avgDepth,
										 (*_bamFiles)[i],
										 " for individual " + coretools::str::toString(i + 1));
		}
	}
	logfile().endIndent();
}

void TBAMSimulator::_simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr,
                                                 std::array<std::vector<Base>, 2> haplotypes,
												 TReadSimulators & readSimulators,
												 uint32_t avgDepth,
                                                 TSimulatorBamFile &bamFile,
												 const std::string &extraProgressText) {
	// Initialize probabilities to simulate reads
	const uint64_t numReads = thisChr.length * avgDepth / readSimulators.averageFragmentLength();
	const uint64_t chrLengthForStart = thisChr.length - readSimulators.maxFragmentLength() + 1;
	const double probReadPerSite     = 1.0 / chrLengthForStart;
	uint64_t numReadsSimulated       = 0;

	// initialize progress reporting
	coretools::TProgressReporter<uint64_t> reporter(numReads, "Simulating about " + coretools::str::toString(numReads) +
	                                                              " reads" + extraProgressText);

	// now simulate
	for(TGenomePosition pos(thisChr.refID(), 0); pos.position() < chrLengthForStart; ++pos){
		// write unwritten alignments
		readSimulators.writeUnwrittenAlignments(pos, bamFile);

		// draw random number to get number of reads starting at this position
		const auto numReadsHere = randomGenerator().getBinomialRand(probReadPerSite, numReads);

		// now simulate
		if (numReadsHere > 0) {
			numReadsSimulated += numReadsHere;
			for (uint32_t r = 0; r < numReadsHere; ++r) {
				readSimulators.simulate(pos, haplotypes[randomGenerator().sample(2)], bamFile);
			}
			// report progress
			reporter.next();
		}
	}
	// write unwritten alignments
	readSimulators.writeUnwrittenAlignments(TGenomePosition(thisChr.refID(), thisChr.length), bamFile);

	reporter.done();
	logfile().conclude("Simulated a total of ", numReadsSimulated, " reads.");
}

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
	_writeSiteInformation(ChrName, Position, VariantQuality, TotalDepth);

	for (size_t i = 0; i < GenotypeLikelihoods.size(); ++i) {
		const bool isMissing = (Depths[i] == 0);
		if (IsDiploid) {
			_writeCell<3>(isMissing, Depths[i],
			              {HighPrecisionPhredIntProbability(GenotypeLikelihoods[i][genometools::BiallelicGenotype::homoFirst]),
			               HighPrecisionPhredIntProbability(GenotypeLikelihoods[i][genometools::BiallelicGenotype::het]),
			               HighPrecisionPhredIntProbability(GenotypeLikelihoods[i][genometools::BiallelicGenotype::homoSecond])});
		} else {
			_writeCell<2>(isMissing, Depths[i],
			              {HighPrecisionPhredIntProbability(GenotypeLikelihoods[i][genometools::BiallelicGenotype::haploidFirst]),
			               HighPrecisionPhredIntProbability(GenotypeLikelihoods[i][genometools::BiallelicGenotype::haploidSecond])});
		}
	}

	// end of line
	_vcf.endln();
}

//-------------------------------------------
// Helper functions
//-------------------------------------------

namespace /* anonymous */ {

auto calculateLog10ProbThisAllelicCombination(const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods) {
	// estimate genotype frequencies
	genometools::TGenotypeFrequencies genotypeFrequencies;
	genotypeFrequencies.estimate<true>(GenotypeLikelihoods, GenotypeLikelihoods.size(), 0.0000001);

	// calculate log10 likelihood
	return genotypeFrequencies.calculateLog10Likelihood(GenotypeLikelihoods, GenotypeLikelihoods.size());
}

auto getRelevantBiallelicGenotype(genometools::Base Major, genometools::Base Ref, bool IsDiploid) {
	if (IsDiploid) {
		return (Major == Ref) ? genometools::BiallelicGenotype::homoFirst
							  : genometools::BiallelicGenotype::homoSecond;
	} else { // haploid
		return (Major == Ref) ? genometools::BiallelicGenotype::haploidFirst
		                      : genometools::BiallelicGenotype::haploidSecond;
	}
}

auto calculateLog10ProbFixed(const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods,
                             genometools::Base Major, genometools::Base Ref, bool IsDiploid) {
	using coretools::Log10Probability;

	auto BG = getRelevantBiallelicGenotype(Major, Ref, IsDiploid);

	Log10Probability LL_fixed_glfPhred = 0.0;
	for (auto GenotypeLikelihood : GenotypeLikelihoods) {
		LL_fixed_glfPhred += (Log10Probability)GenotypeLikelihood[BG];
	}
	return LL_fixed_glfPhred;
}

auto calculateVariantQuality(const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods,
                             genometools::Base Major, genometools::Base Ref, bool IsDiploid) {
	using coretools::Log10Probability;
	const auto LL_allelicCombination = calculateLog10ProbThisAllelicCombination(GenotypeLikelihoods);
	const auto LL_fixed              = calculateLog10ProbFixed(GenotypeLikelihoods, Major, Ref, IsDiploid);

	// calculate variant quality
	return genometools::PhredIntProbability(
		LL_fixed > LL_allelicCombination ? Log10Probability(0.0) : Log10Probability(LL_fixed - LL_allelicCombination));
}

} // end namespace

//-------------------------------------------
// TVCFSimulator
//-------------------------------------------

TVCFSimulator::TVCFSimulator(const std::string &method) : TSimulator(method) {
	// check if method is compatible with VCF: only allow bi-allelic haplotype simulators
	if (!_haploSimulator->simulatesBiallelic()) {
		throw "Can not simulate VCF files with method '" + method +
		    "': only bi-allelic haplotype simulators are allowed (parameter 'type'). Choose other method or simulate "
		    "BAM files instead.";
	}

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
	for (size_t i = 0; i < _haploSimulator->sampleSize(); i++) {
		sampleNames[i] = "ind_" + coretools::str::toString(i);
	}

	// open VCF file and write header
	std::string filename = _outname + ".vcf.gz";
	logfile().list("Writing vcf to file " + filename + "...");
	_vcf = std::make_unique<TVCFWriterSimulation>(filename, "ATLAS_simulations", sampleNames, usePhredLikelihoods);
}

GLF::TMultiGLFDataSampleOneAllelicCombination TVCFSimulator::_calculateGenotypeLikelihoods(size_t NumRef, size_t NumAlt,
                                                                                           bool IsDiploid) {
	const auto P0 =
	    HighPrecisionPhredIntProbability(LogProbability(NumRef * log(_error.complement()) + NumAlt * log(_error)));
	const auto P2 =
	    HighPrecisionPhredIntProbability(LogProbability(NumRef * log(_error) + NumAlt * log(_error.complement())));
	if (IsDiploid) {
		const auto P1 = HighPrecisionPhredIntProbability(LogProbability((NumRef + NumAlt) * lnOneHalf));
		return {P0, P1, P2};
	} else {
		return {P0, P2};
	}
}

size_t TVCFSimulator::_simulateNumReadsWithReferenceAllele(genometools::Base a, genometools::Base b,
                                                           genometools::Base Ref, size_t Depth, bool IsDiploid) {
	if ((IsDiploid && a == b && a == Ref) || (!IsDiploid && a == Ref)) {
		// homozygous reference (haploid and diploid)
		return randomGenerator().getBinomialRand(1. - _error, Depth);
	} else if ((IsDiploid && a == b && a != Ref) || (!IsDiploid && a != Ref)) {
		// homozygous alternative (haploid and diploid)
		return randomGenerator().getBinomialRand(_error, Depth);
	}
	// heterozygous
	return randomGenerator().getBinomialRand(0.5, Depth);
}

std::pair<size_t, GLF::TMultiGLFDataSampleOneAllelicCombination>
TVCFSimulator::_simulateDepthAndGTL(genometools::Base a, genometools::Base b, genometools::Base Ref, bool IsDiploid, uint32_t avgDepth) {
	// simulate depth
	auto depth = randomGenerator().getPoissonRandom(avgDepth);

	// simulate number of reference and alternative alleles
	auto numRef = _simulateNumReadsWithReferenceAllele(a, b, Ref, depth, IsDiploid);
	auto numAlt = (int)depth - numRef;

	// get genotype likelihoods
	auto genotypeLikelihoods = _calculateGenotypeLikelihoods(numRef, numAlt, IsDiploid);

	return std::make_pair(depth, genotypeLikelihoods);
}

std::pair<genometools::Base, genometools::Base>
TVCFSimulator::_findMajorMinorAllele(coretools::TStrongArray<size_t, genometools::Base, 4> AlleleCounts,
                                     genometools::Base RefAllele) {
	// major allele = allele with the highest counts
	const auto majorAllele = randomGenerator()
	                             .sampleIndexOfMaxima<coretools::TStrongArray<size_t, genometools::Base, 4>, Base,
	                                                  coretools::index(genometools::Base::max)>(AlleleCounts);

	if (majorAllele == RefAllele) {
		// choose minor allele: can be any allele (except major)
		// note: multiple alleles might have same counts
		std::vector<genometools::Base> minorAlleles;
		size_t secondMax = 0;
		for (auto b = genometools::Base::min; b < genometools::Base::max; ++b) {
			if (b == majorAllele) { continue; } // exclude major allele
			if (AlleleCounts[b] > secondMax) {  // found new max
				secondMax = AlleleCounts[b];
				minorAlleles.clear();
				minorAlleles.push_back(b);
			} else if (AlleleCounts[b] == secondMax) { // found allele that is equally good
				minorAlleles.push_back(b);
			}
		}
		// randomly sample minor allele, if there are multiple
		auto minorAllele = minorAlleles[randomGenerator().sample(minorAlleles.size())];
		return std::make_pair(majorAllele, minorAllele);
	} else {
		// major allele is not refAllele -> minor allele must be ref-allele!
		// quick check if this is valid (should always be true for bi-allelic simulators)
		for (auto b = genometools::Base::min; b < genometools::Base::max; ++b) {
			if (b != majorAllele && b != RefAllele && AlleleCounts[b] > AlleleCounts[RefAllele]) {
				throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": reference allele (" +
				                         toString(RefAllele) + ") is not major (" + toString(majorAllele) +
				                         ") nor minor (" + toString(AlleleCounts[b]) + ") allele!");
			}
		}
		return std::make_pair(majorAllele, RefAllele);
	}
}

void TVCFSimulator::_simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, uint32_t avgDepth) {
	logfile().startIndent("Simulating genotype likelihoods:");

	for (size_t l = 0; l < Chromosome.length; ++l) {

		// prepare storage
		GLF::TMultiGLFDataOneAllelicCombination genotypeLikelihoods(_haploSimulator->sampleSize());
		std::vector<size_t> depths(_haploSimulator->sampleSize(), 0);
		coretools::TStrongArray<size_t, genometools::Base, 4> alleleCounts({0, 0, 0, 0});
		const bool isDiploid = Chromosome.ploidy == 2;

		for (size_t i = 0; i < _haploSimulator->sampleSize(); ++i) {
			// get haplotypes
			const auto hap1 = Haplotypes(i, 0, l);
			const auto hap2 = Haplotypes(i, 1, l);

			// simulate depth and genotype likelihoods
			auto [depth, GTL] = _simulateDepthAndGTL(hap1, hap2, _reference[l], isDiploid, avgDepth);

			// store
			genotypeLikelihoods[i] = GTL;
			depths[i]              = depth;
			alleleCounts[hap1]++;
			if (isDiploid) { alleleCounts[hap2]++; }
		}

		// find major and minor allele
		const auto refAllele                  = _reference[l];
		const auto [majorAllele, minorAllele] = _findMajorMinorAllele(alleleCounts, refAllele);
		// quick check if ref allele is either major or minor allele. Should always be true if _findMajorMinorAllele is
		// working
		if (refAllele != majorAllele && refAllele != minorAllele) {
			throw std::runtime_error((std::string) __PRETTY_FUNCTION__ + ": Locus " + toString(l) +
			                         ": reference allele (" + toString(refAllele) + ") is not major (" +
			                         toString(majorAllele) + ") nor minor (" + toString(minorAllele) + ") allele!");
		}

		// calculate variant quality
		auto variantQuality = calculateVariantQuality(genotypeLikelihoods, majorAllele, refAllele, isDiploid);

		// finally write site!
		const auto altAllele = (refAllele == majorAllele) ? minorAllele : majorAllele;
		_vcf->writeSite(Chromosome.name, l, variantQuality, refAllele, altAllele, coretools::containerSum(depths),
		                isDiploid, genotypeLikelihoods, depths);
	}
	logfile().endIndent();
}

} // namespace Simulations
