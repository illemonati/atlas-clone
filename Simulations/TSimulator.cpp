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

//---------------------------------------------------------
// Helper functions
//---------------------------------------------------------

namespace /* anonymous */ {
std::unique_ptr<THaplotypeSimulator> makeHaploSimulator(const std::string &method, const TChromosomes &chs) {
	if (method == "one") {
		logfile().list("Simulating a single individual (parameter type=one):");
		return std::make_unique<TSimulatorOne>(chs.size());
	}
	if (method == "pair") {
		logfile().list("Simulating a pair of individual (parameter type=pair):");
		return std::make_unique<TSimulatorPair>();
	}
	if (method == "SFS") {
		logfile().list("Simulating individuals from an SFS (parameter type=SFS):");
		return std::make_unique<TSimulatorSFS>(chs);
	}
	if (method == "HW") {
		logfile().list("Simulating a individuals under Hardy-Weinberg (parameter type=HW):");
		return std::make_unique<TSimulatorHW>();
	}
	throw "Unknown simulation method '" + method + "'!";
}

TChromosomes makeChromosomes() {
	TChromosomes chs;
	chs.clear();
	std::vector<std::string> string_vec;
	parameters().fillParameterIntoContainerWithDefault("chrLength", string_vec, ',', {"1000000"});

	std::vector<uint32_t> chrLength;
	coretools::str::repeatIndexes(string_vec, chrLength);
	if (chrLength.empty()) throw "Issue understanding length of chromosomes!";

	std::vector<uint8_t> ploidy;
	if (parameters().parameterExists("ploidy")) {
		parameters().fillParameterIntoContainer("ploidy", string_vec, ',');
		coretools::str::repeatIndexes(string_vec, ploidy);
	} else {
		for (size_t i = 0; i < chrLength.size(); ++i) ploidy.push_back(2);
	}
	if (ploidy.size() != chrLength.size()) throw "List of chromosome lengths and ploidies differ in length!";
	for (auto &p : ploidy) {
		if (p != 1 && p != 2) { throw "Currently only ploidy 1 (haploid) or 2 (diploid) is supported!"; }
	}

	if (chrLength.size() == 1) {
		const auto numChr = parameters().getParameterWithDefault<int>("numChr", 1);
		std::string text  = "Will simulate " + coretools::str::toString(numChr);
		if (ploidy[0] == 1)
			text += " haploid";
		else
			text += " diploid";
		text += " chromosome(s) of length " + coretools::str::toString(chrLength[0]) + " each.";
		logfile().list(text);
		for (int i = 0; i < numChr; ++i) {
			chs.appendChromosome("chr" + coretools::str::toString(i + 1), chrLength[0], ploidy[0]);
		}
	} else {
		logfile().startIndent("Will simulate ", chrLength.size(), " chromosome(s) of the following length:");
		auto hIt = ploidy.begin();
		std::string text;
		for (auto it = chrLength.begin(); it != chrLength.end(); ++it, ++hIt) {
			text = coretools::str::toString(*it) + " (";
			if (*hIt == 1)
				text += "haploid)";
			else
				text += "diploid)";
			logfile().list(text);
		}

		for (size_t i = 0; i < chrLength.size(); ++i) {
			chs.appendChromosome("chr" + coretools::str::toString(i + 1), chrLength[i], ploidy[i]);
		}
		logfile().endIndent();
	}
	return chs;
}

} // namespace

//---------------------------------------------------------
// TSimulator
//---------------------------------------------------------

TSimulator::TSimulator(const std::string &method)
    : _outname(parameters().getParameterWithDefault<std::string>("out", "ATLAS_simulations")),
      _seqDepth(parameters().getParameterWithDefault("depth", 10.0)),
      _writeTrueGenotypes(parameters().parameterExists("writeTrueGenotypes")),
      _writeVariantInvariantBedFiles(parameters().parameterExists("writeVariantBED")), _reference(_outname + ".fasta"),
      _chromosomes(makeChromosomes()), _haploSimulator(makeHaploSimulator(method, _chromosomes)) {}

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
	for (auto &chr : _chromosomes) {
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
		_simulateAndWrite(chr, haplotypes);

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
	logfile().list("Will simulate to an average depth of ", _seqDepth, ".");
	logfile().list("Will write output files with tag '" + _outname + "'.");
	_initializeReadSimulator();

	// open bam files
	_bamFiles =
	    std::make_unique<TSimulatorBamFiles>(_haploSimulator->sampleSize(), _outname, _readGroups, _chromosomes);
}

void TBAMSimulator::_simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes) {
	// now simulate and write reads
	logfile().startIndent("Simulating reads:");
	for (size_t i = 0; i < _haploSimulator->sampleSize(); ++i)
		_simulateReadsFromHaplotypes(Chromosome, Haplotypes.getHaplotypesOfIndividual(i), (*_bamFiles)[i],
		                             " for individual " + coretools::str::toString(i + 1));
	logfile().endIndent();
}

std::vector<std::string> TBAMSimulator::_readSimInfoPerReadGroup(const std::string &Filename, const std::string &Column,
                                                                 const std::string &Name) {
	// function to initialize read groups
	logfile().listFlush("Reading " + Name + "(s) from file '" + Filename + "' ...");

	coretools::TInputFile in(Filename, {"ReadGroup", Column}, "\t", "//");
	std::vector<std::string> vec;
	std::vector<bool> found(_readGroups.size(), false);
	std::vector<std::string> ret(_readGroups.size());

	// now parse file
	while (in.read(vec)) {
		// find read group
		const auto rg = _readGroups.getId(vec[0]);
		found[rg]     = true;
		ret[rg]       = vec[1];
	}
	logfile().done();
	logfile().conclude("Read " + Name + "s for ", in.lineNumber(), " read groups.");

	// check if there was data for each read group
	for (size_t i = 0; i < found.size(); ++i) {
		if (!found[i])
			throw "No " + Name + " given for read group '" + _readGroups.getName(i) + "' in file '" + Filename + "'!";
	}
	return ret;
}

void TBAMSimulator::_initializeReadGroup(const std::string &readLengthString, const BAM::TReadGroup &ReadGroup) {
	// single or paired end? Is indicated at beginning of readLengthString!
	if (readLengthString.find("single:") == 0) {
		_readSimulators.push_back(std::make_unique<TSimulatorSingleEndRead>(ReadGroup));
	} else if (readLengthString.find("paired:") == 0) {
		_readSimulators.push_back(std::make_unique<TSimulatorPairedEndReads>(ReadGroup));
	} else
		throw "Unable to understand string '" + readLengthString + "'!";

	// add read Length distribution
	const auto readLengthDist = coretools::str::readAfterLast(readLengthString, ':');
	_readSimulators.back()->setReadLengthDistribution(readLengthDist);
}

void TBAMSimulator::_initializeReadGroupsFromReadLengthDistribution(const std::string &ParameterName,
                                                                    const std::string &DefaultValue,
                                                                    const std::string &Name) {
	logfile().startIndent("Parsing read length distribution (parameter '" + ParameterName + "'):");
	const auto s = parameters().getParameterWithDefault<std::string>(ParameterName, DefaultValue);
	_readSimulators.clear();

	// We allow for two options:
	//   1) initialized from the command line (one for all read groups)
	//   2) read-group specific as given in a file

	// check if it is a file (should not contain a ':')
	const auto pos = s.find(':');
	if (pos != std::string::npos) {
		// Option 1: a single read length distribution for all
		//---------------------------------------------------------------------
		logfile().list("Will use '" + s + " for all read groups.");

		// create read groups
		for (auto &r : _readGroups) { _initializeReadGroup(s, r); }
	} else {
		// Option 2: read group specific, given in a file
		//---------------------------------------------------------------------
		const std::vector<std::string> dist = _readSimInfoPerReadGroup(s, ParameterName, Name);

		for (size_t r = 0; r < dist.size(); ++r) { _initializeReadGroup(dist[r], _readGroups[r]); }
	}
	logfile().endIndent();
}

void TBAMSimulator::_initializeDistribution(const std::string &ParameterName, const std::string &DefaultValue,
                                            const std::string &Name,
                                            std::function<void(TSimulatorSingleEndRead &, std::string)> function) {
	logfile().startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");
	const auto s = parameters().getParameterWithDefault<std::string>(ParameterName, DefaultValue);

	// We allow for two options:
	//   1) initialized from the command line (one for all read groups)
	//   2) read-group specific as given in a file

	// check if it is a file (should not contain a '(')
	const auto pos = s.find('(');
	if (pos != std::string::npos) {
		// Option 1: a single read distribution for all
		//---------------------------------------------------------------------
		logfile().list("Will use '" + s + " for all read groups.");

		// create read groups
		for (auto &r : _readSimulators) { function(*r, s); }
	} else {
		// Option 2: read group specific, given in a file
		//---------------------------------------------------------------------
		const std::vector<std::string> dist = _readSimInfoPerReadGroup(s, ParameterName, Name);

		for (uint32_t r = 0; r < _readSimulators.size(); ++r) { function(*_readSimulators[r], dist[r]); }
	}
	logfile().endIndent();
}

void TBAMSimulator::_initializePMD(const std::string &ParameterName, const std::string &Name) {

	logfile().startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");

	if (parameters().parameterExists(ParameterName)) {
		const auto pmdString = parameters().getParameter<std::string>(ParameterName);
		std::vector<uint16_t> ReadGroupsWithoutPMD;
		_PMD.initialize(pmdString, _readGroups, ReadGroupsWithoutPMD);

		// add PMD to simulators
		for (size_t r = 0; r < _readSimulators.size(); ++r) { _readSimulators[r]->setPMD(&_PMD[r]); }
	} else {
		logfile().list("Not simulating any PMD.");
	}
	logfile().endIndent();
}

void TBAMSimulator::_initializeQualityTransformations(const std::string &ParameterName, const std::string &Name) {
	logfile().startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");

	if (parameters().parameterExists(ParameterName)) {
		const std::string rhoString = parameters().getParameterWithDefault<std::string>("rho", "default");
		const auto recalString = parameters().getParameter<std::string>(ParameterName);
		_recal.initialize(recalString, rhoString, _readGroups);

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

void addReadGroupsIfFile(const std::string &ParameterName, BAM::TReadGroups &ReadGroups) {
	// check if parameter is given
	if (parameters().parameterExists(ParameterName)) {
		const auto s = parameters().getParameter<std::string>(ParameterName);

		// check if string s provides a definition (contains a ':') or is a file (does not contain a ':')
		if (!coretools::str::stringContains(s, ":")) {
			// is probably a file -> try to open it
			if (std::filesystem::exists(s)) {
				coretools::TInputFile in(s, {"ReadGroup"}, "\t", "//");
				std::vector<std::string> tmp;
				while (in.read(tmp)) {
					// add all non-existing elemets
					ReadGroups.add(tmp[0]);
				}
			}
		}
	}
}

void TBAMSimulator::_initializeReadSimulator() {
	// For which read groups?
	// Check for each parameter if it is given per read group (a file) or common to all
	_readGroups.clear();
	addReadGroupsIfFile("readLength", _readGroups);
	addReadGroupsIfFile("qualityDist", _readGroups);
	addReadGroupsIfFile("pmd", _readGroups);
	addReadGroupsIfFile("recal", _readGroups);
	addReadGroupsIfFile("readGroupFreq", _readGroups);

	// any read groups specified?
	if (_readGroups.empty()) {
		const auto numRG = parameters().getParameterWithDefault<int>("numReadGroups", 1);
		for (int i = 0; i < numRG; ++i) { _readGroups.add("SimReadGroup" + coretools::str::toString(i + 1)); }
		// report
		if (numRG == 1) {
			logfile().startIndent("Initializing one read group (parameter 'numReadGroups'):");
		} else if (numRG > 1) {
			logfile().startIndent("Initializing ", numRG, " identical read groups (parameter 'numReadGroups'):");
		} else {
			throw "numReadGroups must be at least 1!";
		}

	} else {
		logfile().startIndent("Initializing ", _readGroups.size(), " individual read group(s):");
	}

	// A) read length: used to create simulator read groups
	//---------------
	_initializeReadGroupsFromReadLengthDistribution("readLength", "single:fixed(100)", "read length");

	// B) initialize quality distribution
	//-----------------------------------
	_initializeDistribution("baseQuality", "normal(30,10)[0,93]", "base quality distribution",
	                        &TSimulatorSingleEndRead::setQualityDistribution);

	// C) initialize mapping quality distribution
	//-----------------------------------
	_initializeDistribution("mappingQuality", "normal(60,10)[1,255]", "mapping quality distribution",
	                        &TSimulatorSingleEndRead::setMappingQualityDistribution);

	// D) initialize PMD
	//------------------
	_initializePMD("pmd", "post-mortem damage");

	// E) initialize quality transformation
	//-------------------------------------
	_initializeQualityTransformations("recal", "recalibration models");

	// E) initialize contamination
	//----------------------------
	// TODO: Think about contamination object for both estimation and simulation

	// F) other things
	//----------------

	// initialize read group frequencies frequencies
	_initializeReadGroupFrequencies();

	logfile().endIndent();
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
}

void TBAMSimulator::_simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr,
                                                 std::array<std::vector<Base>, 2> haplotypes,
                                                 TSimulatorBamFile &bamFile, const std::string &extraProgressText) {
	// Initialize probabilities to simulate reads
	const uint64_t numReads = _averageReadLength == 0 ? 0 : thisChr.length * _seqDepth / _averageReadLength;

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
TVCFSimulator::_simulateDepthAndGTL(genometools::Base a, genometools::Base b, genometools::Base Ref, bool IsDiploid) {
	// simulate depth
	auto depth = randomGenerator().getPoissonRandom(_seqDepth);

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

void TVCFSimulator::_simulateAndWrite(const genometools::TChromosome &Chromosome, TSimulatorHaplotypes &Haplotypes) {
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
			auto [depth, GTL] = _simulateDepthAndGTL(hap1, hap2, _reference[l], isDiploid);

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
