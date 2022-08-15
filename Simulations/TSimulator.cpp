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

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TChromosomes.h"
#include "TFile.h"
#include "TGenotypeFrequencies.h"
#include "THaplotypeSimulator.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "algorithms.h"
#include "gzstream.h"
#include "mathConstants.h"
#include "progressTools.h"
#include "stringFunctions.h"
#include "weakTypes.h"

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

//---------------------------------------------------------
// Helper functions
//---------------------------------------------------------

namespace /* anonymous */ {
std::unique_ptr<THaplotypeSimulator> makeHaploSimulator(const std::string &method, const TChromosomes &chs) {
	if (method == "one") {
		logfile().list("Simulating a single individual (parameter --type one):");
		return std::make_unique<TSimulatorOne>(chs.size());
	}
	if (method == "pair") {
		logfile().list("Simulating a pair of individual (parameter --type pair):");
		return std::make_unique<TSimulatorPair>();
	}
	if (method == "SFS") {
		logfile().list("Simulating individuals from an SFS (parameter --type SFS):");
		return std::make_unique<TSimulatorSFS>(chs);
	}
	if (method == "HW") {
		logfile().list("Simulating a individuals under Hardy-Weinberg (parameter --type HW):");
		return std::make_unique<TSimulatorHW>();
	}
	throw "Unknown simulation method '" + method + "'!";
}



std::vector<uint8_t> parsePloidy(){
	//parse ploidy parameters
	std::vector<std::string> string_vec;
	std::vector<uint8_t> ploidy;
	if (parameters().parameterExists("ploidy")) {
		parameters().fillParameterIntoContainer("ploidy", string_vec, ',');
		coretools::str::repeatIndexes(string_vec, ploidy);
	} else {
		ploidy.push_back(2); //default ploidy = 2
		return ploidy;
	}

	//check if ploidy is supported
	for (auto &p : ploidy) {
		if (p != 1 && p != 2) { throw "Currently only ploidy 1 (haploid) or 2 (diploid) is supported!"; }
	}

	return ploidy;
}

std::vector<uint32_t> parseSeqDepth(){
	//parse ploidy parameters
	std::vector<std::string> string_vec;
	std::vector<uint32_t> depth;
	if (parameters().parameterExists("depth")) {
		parameters().fillParameterIntoContainer("depth", string_vec, ',');
		coretools::str::repeatIndexes(string_vec, depth);
	} else {
		depth.push_back(10); //default depth = 10
		return depth;
	}

	return depth;
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
	chs.clear();
	depths.clear();

	//parse chromosome lengths
	std::vector<std::string> string_vec;
	parameters().fillParameterIntoContainerWithDefault("chrLength", string_vec, ',', {"1000000"});

	std::vector<uint32_t> chrLengths;
	coretools::str::repeatIndexes(string_vec, chrLengths);
	if (chrLengths.empty()) throw "Issue understanding length of chromosomes!";

	//parse ploidies and depth
	std::vector<uint8_t> ploidies = parsePloidy();
	depths = parseSeqDepth();

	//check length
	size_t numChr = std::max( {chrLengths.size(), ploidies.size(), depths.size()} );
	checkLength(chrLengths, numChr);
	checkLength(ploidies, numChr);
	checkLength(depths, numChr);

	//report and create
	logfile().startIndent("Will simulate ", chrLengths.size(), " chromosome(s):");
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
}

} // namespace

//---------------------------------------------------------
// TSimulator
//---------------------------------------------------------

TSimulator::TSimulator(const std::string &method)
    : _outname(parameters().getParameterWithDefault<std::string>("out", "ATLAS_simulations")),
      _writeTrueGenotypes(parameters().parameterExists("writeTrueGenotypes")),
      _writeVariantInvariantBedFiles(parameters().parameterExists("writeVariantBED")), _reference(_outname + ".fasta")
{
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
	// depth
	logfile().list("Will write output files with tag '" + _outname + "'.");
	_initializeReadSimulator();

	// open bam files
	_bamFiles =
	    std::make_unique<TSimulatorBamFiles>(_haploSimulator->sampleSize(), _outname, _readGroups, _chromosomes);
}

void TBAMSimulator::_simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes, uint32_t avgDepth) {
	// now simulate and write reads
	logfile().startIndent("Simulating reads:");
	for (size_t i = 0; i < _haploSimulator->sampleSize(); ++i)
		_simulateReadsFromHaplotypes(Chromosome, Haplotypes.getHaplotypesOfIndividual(i), avgDepth, (*_bamFiles)[i],
		                             " for individual " + coretools::str::toString(i + 1));
	logfile().endIndent();
}

void TBAMSimulator::_initializeReadGroups(const Simulations::RGInfo::TSimulatorReadGroupInfo & RGinfo) {
	//fill TReadGroups
	RGinfo.createReadGroups(_readGroups);

	// create simulation read groups
	for(size_t i = 0; i < RGinfo.size(); ++i){
		std::string type = RGinfo[i].get(RGInfo::InfoType::seqType);
		if(type == "single"){
			_readSimulators.push_back(std::make_unique<TSimulatorSingleEndRead>(_readGroups[i], RGinfo[i]));
		} else if(type == "paired"){
			DEVERROR("Paired-end read groups not yet implemented!");
			//_readSimulators.push_back(std::make_unique<TSimulatorSingleEndRead>(RGinfo.readGroup(i), RGinfo[i]));
		} else {
			UERROR("Unable to understand read group type '" + type + "'! Use either 'single' or 'paired'.");
		}
	}
}

void TBAMSimulator::_initializePMD(){
	const std::string arg = "pmd";
	if (parameters().parameterExists(arg)) {
		const auto pmdString = parameters().getParameter<std::string>(arg);
		std::vector<uint16_t> ReadGroupsWithoutPMD;
		_PMD.initialize(pmdString, _readGroups, ReadGroupsWithoutPMD);

		// add PMD to simulators
		for (size_t r = 0; r < _readSimulators.size(); ++r) { _readSimulators[r]->setPMD(&_PMD[r]); }
	} else {
		logfile().list("Not simulating any PMD.");
	}
	logfile().endIndent();
}

void TBAMSimulator::_initializeQualityTransformations() {
	const std::string arg = "recal";
	if (parameters().parameterExists(arg)) {
		const std::string rhoString = parameters().getParameterWithDefault<std::string>("rho", "default");
		const auto recalString = parameters().getParameter<std::string>(ParameterName);
		_recal.initialize(recalString, rhoString, _readGroups);
		logfile().list("Will use '", recalString, "' for all read groups.");

		// add recal to simulators
		for (size_t r = 0; r < _readSimulators.size(); ++r) {
			_readSimulators[r]->setRecal(&_recal(r, false), &_recal(r, true));
		}
	} else {
		// add noRecal model. Is still needed for simulation of bases from base qualities!
		for (size_t r = 0; r < _readSimulators.size(); ++r) {
			_readSimulators[r]->setRecal(&_recal(r, false), &_recal(r, true));
		}
		logfile().list("Not simulating any quality transformation.");
	}
	logfile().endIndent();
}

void TBAMSimulator::_initializeReadSimulator() {
	// A) initialize read groups from RG Info / Command line
	Simulations::RGInfo::TSimulatorReadGroupInfo RGinfo;
	_initializeReadGroups(RGinfo);

	// B) initialize PMD
	//------------------
	_initializePMD();

	// E) initialize quality transformation
	//-------------------------------------
	_initializeQualityTransformations();

	// E) initialize contamination
	//----------------------------
	// TODO: Think about contamination object for both estimation and simulation

	// G) other things
	//----------------
	// initialize read group frequencies frequencies
	_initializeReadGroupFrequencies();

	logfile().endIndent();

	//report all read groups
	_printSimulationDetailsAllReadGroups();
}

void TBAMSimulator::_initializeReadGroupFrequencies() {
	_cumulSimGroupFrequenies.reserve(_readSimulators.size());
	_simGroupFrequencies.reserve(_readSimulators.size());
	if (parameters().parameterExists("readGroupFreq")) {
		// read frequencies
		std::vector<std::string> vec;
		parameters().fillParameterIntoContainer("readGroupFreq", vec, true);
		std::vector<double> freq;
		coretools::str::repeatIndexes(vec, freq);
		if (freq.size() != _readSimulators.size())
			throw "Provided read group frequencies do not match number of read groups!";

		// normalize and print
		const auto sum = std::accumulate(freq.cbegin(), freq.cend(), 0.);

		logfile().startIndent("Will simulate read groups with the following frequencies:");
		for (size_t i = 0; i < _readSimulators.size(); ++i) {
			_simGroupFrequencies[i] = freq[i] / sum;
			logfile().list(_simGroupFrequencies[i], " " + _readSimulators[i]->name());
		}
		logfile().endIndent();

		// fill cumulative
		_cumulSimGroupFrequenies[0] = _simGroupFrequencies[0];
		for (size_t i = 1; i < _readSimulators.size(); ++i)
			_cumulSimGroupFrequenies[i] = _cumulSimGroupFrequenies[i - 1] + _simGroupFrequencies[i];
		_cumulSimGroupFrequenies[_readSimulators.size() - 1] = 1.0; // ensure last entry is 1.0
	} else {
		// equal frequencies
		logfile().list("Will simulate reads equally distributed among read groups.");
		for (size_t i = 0; i < _readSimulators.size(); ++i) {
			_simGroupFrequencies[i]     = 1.0 / _readSimulators.size();
			_cumulSimGroupFrequenies[i] = (double)(i + 1) / _readSimulators.size();
		}
	}

	// precalculate some stuff
	_averageReadLength = 0;
	_maxReadLength     = 0;

	for (size_t i = 0; i < _readSimulators.size(); ++i) {
		_averageReadLength += _simGroupFrequencies[i] * _readSimulators[i]->meanReadLength();
		if (_readSimulators[i]->maxReadLength() > _maxReadLength) _maxReadLength = _readSimulators[i]->maxReadLength();
	}

	//check if read length match chr length
	for(auto& chr : _chromosomes){
		if(_maxReadLength > chr.length){
			throw "Length of chromosome '" + chr.name + "' is less than the max read length!";

		}
	}
}

void TBAMSimulator::_printSimulationDetailsAllReadGroups(){
	logfile().startIndent("Will simulate the following read groups:");

	for(size_t rs = 0; rs < _readSimulators.size(); ++rs){
		_readSimulators[rs]->printDetails(_simGroupFrequencies[rs]);
	}

	//TODO write file with settings
}

void TBAMSimulator::_simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr,
                                                 std::array<std::vector<Base>, 2> haplotypes,
												 uint32_t avgDepth,
                                                 TSimulatorBamFile &bamFile,
												 const std::string &extraProgressText) {
	// Initialize probabilities to simulate reads
	const uint64_t numReads = _averageReadLength == 0 ? 0 : thisChr.length * avgDepth / _averageReadLength;

	const uint64_t chrLengthForStart = thisChr.length - _maxReadLength + 1;
	const double probReadPerSite     = 1.0 / chrLengthForStart;
	uint64_t numReadsSimulated       = 0;

	// initialize progress reporting
	coretools::TProgressReporter<uint64_t> reporter(numReads, "Simulating about " + coretools::str::toString(numReads) +
	                                                              " reads" + extraProgressText);

	// now simulate
	for (uint32_t l = 0; l < chrLengthForStart; ++l) {
		// write unwritten alignments
		for (auto &rs : _readSimulators) rs->writeUnwrittenAlignments(l, bamFile);

		// draw random number to get number of reads starting at this position
		const auto numReadsHere = randomGenerator().getBinomialRand(probReadPerSite, numReads);

		// now simulate
		if (numReadsHere > 0) {
			numReadsSimulated += numReadsHere;
			for (uint32_t r = 0; r < numReadsHere; ++r) {
				const auto rg = randomGenerator().pickOne(_readSimulators.size(), _cumulSimGroupFrequenies.data());
				_readSimulators[rg]->simulate(haplotypes[randomGenerator().sample(2)], thisChr.refID(), l, bamFile);
			}
			// report progress
			reporter.next();
		}
	}
	// write unwritten alignments
	for (auto &rs : _readSimulators) rs->writeUnwrittenAlignments(thisChr.length, bamFile);

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
			              {GenotypeLikelihoods[i][genometools::BiallelicGenotype::homoFirst],
			               GenotypeLikelihoods[i][genometools::BiallelicGenotype::het],
			               GenotypeLikelihoods[i][genometools::BiallelicGenotype::homoSecond]});
		} else {
			_writeCell<2>(isMissing, Depths[i],
			              {GenotypeLikelihoods[i][genometools::BiallelicGenotype::haploidFirst],
			               GenotypeLikelihoods[i][genometools::BiallelicGenotype::haploidSecond]});
		}
	}

	// end of line
	_vcf << '\n';
}

//-------------------------------------------
// Helper functions
//-------------------------------------------

namespace /* anonymous */ {

auto calculateLog10ProbThisAllelicCombination(const GLF::TMultiGLFDataOneAllelicCombination &GenotypeLikelihoods) {
	// estimate genotype frequencies
	genometools::TGenotypeFrequencies genotypeFrequencies;
	genotypeFrequencies.estimate(GenotypeLikelihoods, GenotypeLikelihoods.size(), 0.0000001);

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
	const auto LL_allelicCombination = calculateLog10ProbThisAllelicCombination(GenotypeLikelihoods);
	const auto LL_fixed              = calculateLog10ProbFixed(GenotypeLikelihoods, Major, Ref, IsDiploid);

	// calculate variant quality
	return genometools::PhredIntProbability(LL_fixed > LL_allelicCombination ? coretools::Log10Probability(0.0)
	                                                                         : LL_fixed - LL_allelicCombination);
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
	                             .sampleIndexOfMaxima<coretools::TStrongArray<size_t, genometools::Base, 4>,
	                                                  genometools::Base, index(genometools::Base::max)>(AlleleCounts);

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
