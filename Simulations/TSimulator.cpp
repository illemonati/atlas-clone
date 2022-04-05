/*
 * TSimulator.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#include "GenotypeTypes.h"

#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TSimulator.h"
#include <memory>
#include <numeric>

namespace Simulations {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using genometools::Base;

TSimulator::TSimulator()
    : _outname(parameters().getParameterWithDefault<std::string>("out", "ATLAS_simulations")),
      _seqDepth(parameters().getParameterWithDefault("depth", 10.0)),
      _writeTrueGenotypes(parameters().parameterExists("writeTrueGenotypes")),
      _writeVariantInvariantBedFiles(parameters().parameterExists("writeVariantBED")), _reference(_outname + ".fasta") {
}

//---------------------------------------------------
// TBamSimulator
//---------------------------------------------------

TBAMSimulator::TBAMSimulator() : TSimulator() {
	using genometools::Base;
	// depth
	logfile().list("Will simulate to an average depth of ", _seqDepth, ".");
	logfile().list("Will write output files with tag '" + _outname + "'.");

	_initializeReadSimulator();
	_initializeChromosomes();
}

//--------------------------------------------------------------
// Function to initialize read groups
//--------------------------------------------------------------
std::vector<std::string> TBamSimulator::_readSimInfoPerReadGroup(const std::string &Filename, const std::string &Column,
								 const std::string &Name) {
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

void TBamSimulator::_initializeReadGroup(const std::string &readLengthString, const BAM::TReadGroup &ReadGroup) {
	// single or paired end? Is indicated at beginning of readLengthString!
	if (readLengthString.find("single:") == 0) {
		_readSimulators.push_back(std::make_unique<TBamSimulatorSingleEndRead>(ReadGroup));
	} else if (readLengthString.find("paired:") == 0) {
		_readSimulators.push_back(std::make_unique<TBamSimulatorPairedEndReads>(ReadGroup));
	} else
		throw "Unable to understand string '" + readLengthString + "'!";

	// add read Length distribution
	const auto readLengthDist = coretools::str::readAfterLast(readLengthString, ':');
	_readSimulators.back()->setReadLengthDistribution(readLengthDist);
}

void TBamSimulator::_initializeReadGroupsFromReadLengthDistribution(const std::string &ParameterName,
								    const std::string &DefaultValue,
								    const std::string &Name) {
	logfile().startIndent("Parsing read length distribution (parameter '" + ParameterName + "'):");
	const auto s = parameters().getParameterWithDefault<std::string>(ParameterName, DefaultValue);
	_readSimulators.clear();

	// We allow for two options:
	//   1) initialized from the command line (one for all read groups)
	//   2) read-group specific as given in a file

	// check if it is a file (should not contain a ':')
	const auto pos = s.find(":");
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

void TBamSimulator::_initializeDistribution(const std::string &ParameterName, const std::string &DefaultValue,
					    const std::string &Name,
					    std::function<void(TBamSimulatorSingleEndRead &, std::string)> function) {
	logfile().startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");
	const auto s = parameters().getParameterWithDefault<std::string>(ParameterName, DefaultValue);

	// We allow for two options:
	//   1) initialized from the command line (one for all read groups)
	//   2) read-group specific as given in a file

	// check if it is a file (should not contain a '(')
	const auto pos = s.find("(");
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

void TBamSimulator::_initializePMD(const std::string &ParameterName, const std::string &Name) {

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
}

void TBamSimulator::_initializeQualityTransformations(const std::string &ParameterName, const std::string &Name) {
	logfile().startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");

	if (parameters().parameterExists(ParameterName)) {
		const auto recalString = parameters().getParameter<std::string>(ParameterName);
		_recal.initializeFromFile(recalString, _readGroups);

		// add recal to simulators
		for (size_t r = 0; r < _readSimulators.size(); ++r) {
			_readSimulators[r]->setRecal(&_recal(r, 0), &_recal(r, 1));
		}
	} else {
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

void TBamSimulator::_initializeReadSimulator() {
	// For which read groups?
	// Check for each parameter if it is given per read group (a file) or common to all
	_readGroups.clear();
	_addReadGroupsIfFile("readLength", _readGroups);
	_addReadGroupsIfFile("qualityDist", _readGroups);
	_addReadGroupsIfFile("pmd", _readGroups);
	_addReadGroupsIfFile("recal", _readGroups);
	_addReadGroupsIfFile("readGroupFreq", _readGroups);

	// any read groups specified?
	if (_readGroups.empty()) {
		const auto numRG = parameters().getParameterWithDefault<int>("numReadGroups", 1);
		for (int i = 0; i < numRG; ++i) { _readGroups.add("SimReadGroup" + coretools::str::toString(i + 1)); }
		// report
		if (numRG == 1) {
			logfile().startIndent("Initializing one read group (parameter 'numreadGroups'):");
		} else if (numRG > 1) {
			logfile().startIndent("Initializing ", numRG, " identical read groups (parameter 'numreadGroups'):");
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
				&TBamSimulatorSingleEndRead::setQualityDistribution);

	// C) initialize mapping quality distribution
	//-----------------------------------
	_initializeDistribution("mappingQuality", "normal(60,10)[1,255]", "mapping quality distribution",
				&TBamSimulatorSingleEndRead::setMappingQualityDistribution);

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
}

void TBamSimulator::_initializeReadGroupFrequencies() {
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

//--------------------------------------------------------------
// Initialize chromosomes, depth and base frequencies
//--------------------------------------------------------------
void TBamSimulator::_initializeChromosomes() {
	_chromosomes.clear();
	std::vector<std::string> string_vec;
	parameters().fillParameterIntoContainerWithDefault("chrLength", string_vec, ',', {"1000000"});

	std::vector<uint32_t> chrLength;
	coretools::str::repeatIndexes(string_vec, chrLength);
	if (chrLength.size() < 1) throw "Issue understanding length of chromosomes!";

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
			_chromosomes.appendChromosome("chr" + coretools::str::toString(i + 1), chrLength[0], ploidy[0]);
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
			_chromosomes.appendChromosome("chr" + coretools::str::toString(i + 1), chrLength[i], ploidy[i]);
		}
		logfile().endIndent();
	}
}

//--------------------------------------------------------------
// Run simulations
//--------------------------------------------------------------
Base sampleBase(const std::array<double, 4> &cumulProbs) {
	return genometools::Base(randomGenerator().pickOne(cumulProbs));
}

Base mutateBase(Base base, const std::array<double, 4> &cumulProbs) {
	using namespace genometools;
	return Base((index(base) + randomGenerator().pickOne(cumulProbs)) % index(Base::max));
}

void TBamSimulator::_simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr,
						 std::array<std::vector<Base>, 2> haplotypes,
						 TBamSimulatorBamFile &bamFile, std::string extraProgressText) {
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

void TBamSimulator::runSimulations() {
	// open bam files
	TBamSimulatorBamFiles bamFiles(_sampleSize, _outname, _readGroups, _chromosomes);

	// prepare haplotypes and
	TBamSimulatorHaplotypes haplotypes(_sampleSize);

	// open files to store extra info on sites
	if (_writeTrueGenotypes) {
		// open file for true genotypes
		const auto filename = _outname + "_trueGenotypes.vcf.gz";
		haplotypes.openTrueGenotypeVCF(filename);
	}

	TBamSimulatorVariantInvariantBedFiles bedFiles;
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
			_simulateHaplotypesHaploid(haplotypes, chr);
		else
			_simulateHaplotypesDiploid(haplotypes, chr);
		logfile().done();

		// write true genotypes
		if (_writeTrueGenotypes) {
			logfile().listFlush("Writing true genotypes ...");
			haplotypes.writeTrueGenotypes(chr.name, _reference);
			logfile().done();
		}

		// write BED files
		if (_writeVariantInvariantBedFiles) bedFiles.write(haplotypes, chr.name);

		// now simulate and write reads
		logfile().startIndent("Simulating reads:");
		for (int i = 0; i < _sampleSize; ++i)
			_simulateReadsFromHaplotypes(chr, haplotypes.getHaplotypesOfIndividual(i), bamFiles[i],
			                             " for individual " + coretools::str::toString(i + 1));
		logfile().endIndent();

		// end of chromosome
		logfile().endIndent();
	}
	logfile().endIndent();
}

} // namespace Simulations
