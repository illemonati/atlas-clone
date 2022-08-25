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
		logfile().startIndent("Simulating a individuals under Hardy-Weinberg (parameter 'type' = '", method, "')");
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

void TBAMSimulator::_initializeReadGroups(const TReadGroupInfo & RGinfo) {
	// create simulation read groups
	using BAM::RGInfo::InfoType;
	for(size_t i = 0; i < RGinfo.size(); ++i){
		logfile().startIndent("Read group '", RGinfo[i][InfoType::RGName], "':");
		std::string type = RGinfo[i][InfoType::seqType];
		logfile().list("Sequencing type: ", type);
		logfile().list("Frequency: ", _simGroupFrequencies[i]);

		//initialize by type
		if(type == "single"){
			_readSimulators.push_back(std::make_unique<TSimulatorSingleEndRead>(_readGroups[i], RGinfo[i]));
		} else if(type == "paired"){
			_readSimulators.push_back(std::make_unique<TSimulatorPairedEndReads>(_readGroups[i], RGinfo[i]));
		} else {
			UERROR("Unable to understand read group type '" + type + "'! Use either 'single' or 'paired'.");
		}
		logfile().endIndent();
	}
}

void TBAMSimulator::_initializePMD(){
	const std::string arg = "pmd";
	if (parameters().parameterExists(arg)) {
		const auto pmdString = parameters().getParameter<std::string>(arg);
		_PMD.initialize(pmdString, _readGroups);

		// add PMD to simulators
		for (size_t r = 0; r < _readSimulators.size(); ++r) { _readSimulators[r]->setPMD(&_PMD[r]); }
	} else {
		logfile().list("Not simulating any PMD.");
	}
}

void TBAMSimulator::_initializeQualityTransformations() {
	const std::string arg = "recal";
	if (parameters().parameterExists(arg)) {
		const std::string rhoString = parameters().getParameterWithDefault<std::string>("rho", "default");
		const auto recalString = parameters().getParameter<std::string>(arg);
		_recal.initialize(recalString, rhoString, _readGroups);
		logfile().list("Will use '", recalString, "' for all read groups.");

		// add recal to simulators
		for (size_t r = 0; r < _readSimulators.size(); ++r) {
			_readSimulators[r]->setRecal(&_recal(r, false), &_recal(r, true));
		}
	} else {
		_recal.initializeNoRecal(_readGroups);
		// add noRecal model. Is still needed for simulation of bases from base qualities!
		_recal.initializeNoRecal(_readGroups);
		for (size_t r = 0; r < _readSimulators.size(); ++r) {
			_readSimulators[r]->setRecal(&_recal(r, false), &_recal(r, true));
		}
		logfile().list("Not simulating any quality transformation.");
	}
	logfile().endIndent();
}

void TBAMSimulator::_initializeReadSimulator() {
	logfile().startIndent("Parameters regarding sequencing:");
	// Read sequencing parameters from RG Info / Command line
	TReadGroupInfo RGinfo;
	_readGroups = RGinfo.readInfoAndCreateReadGroups();

	using BAM::RGInfo::InfoType;
	RGinfo.parse(InfoType::RGFrequency, InfoType::seqType, InfoType::cycles, InfoType::fragmentLength, InfoType::baseQuality, InfoType::mappingQuality, InfoType::softClipping);
	_initializeReadGroupFrequencies(RGinfo);
	logfile().endIndent();

	//Initialize read groups
	logfile().startIndent("Initializing ", _readGroups.size(), " read groups:");
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

	logfile().endIndent();

	//warn if read group info columns were not used
	RGinfo.warnAboutUnusedColumnsInFile();
	logfile().endIndent();

	//prepare simulations
	_prepareSimulations();
}

void TBAMSimulator::_initializeReadGroupFrequencies(const TReadGroupInfo & RGinfo) {
	_cumulSimGroupFrequenies.resize(RGinfo.size());
	_simGroupFrequencies.resize(RGinfo.size());

	using BAM::RGInfo::InfoType;
	if(RGinfo.hasInfo(InfoType::RGFrequency)){
		//fill frequencies and cumulative frequencies
		std::vector<double> tmp;
		RGinfo.fillContainerPerReadGroup(tmp, InfoType::RGFrequency);
		coretools::fillFromNormalized(_simGroupFrequencies, tmp);
	} else {
		Probability equal = 1.0 / (double) RGinfo.size();
		for (size_t i = 0; i < RGinfo.size(); ++i) {
			_simGroupFrequencies[i] = equal;
		}
	}
	coretools::fillCumulative(_simGroupFrequencies, _cumulSimGroupFrequenies);
}

void TBAMSimulator::_prepareSimulations(){
	// precalculate some stuff
	_averageReadLength = 0;
	_maxFragmentLength = 0;

	for (size_t i = 0; i < _readSimulators.size(); ++i) {
		_averageReadLength += _simGroupFrequencies[i] * _readSimulators[i]->meanReadLength();
		if (_readSimulators[i]->maxFragmentLength() > _maxFragmentLength){
			_maxFragmentLength = _readSimulators[i]->maxFragmentLength();
		}
	}

	//check if read length match chr length
	for(auto& chr : _chromosomes){
		if(_maxFragmentLength > chr.length){
			throw "Length of chromosome '" + chr.name + "' is less than the max fragment length!";

		}
	}
}

void TBAMSimulator::_simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr,
                                                 std::array<std::vector<Base>, 2> haplotypes,
												 uint32_t avgDepth,
                                                 TSimulatorBamFile &bamFile,
												 const std::string &extraProgressText) {
	// Initialize probabilities to simulate reads
	const uint64_t numReads = _averageReadLength == 0 ? 0 : thisChr.length * avgDepth / _averageReadLength;
	const uint64_t chrLengthForStart = thisChr.length - _maxFragmentLength + 1;
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
