/*
 * TSimulator.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#include "GenotypeTypes.h"

#include "TChromosomes.h"
#include "THaplotypeSimulator.h"
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
using genometools::TChromosomes;

namespace /* anonymous */ {
std::unique_ptr<THaplotypeSimulator> makeHaploSimulator(const std::string &method, const TChromosomes chs) {
	if (method == "one") {
		logfile().startIndent("Simulating a single individual (parameter type=one):");
		return std::make_unique<TSimulatorOne>(chs.size());
	}
	if (method == "pair") {
		logfile().startIndent("Simulating a pair of individual (parameter type=pair):");
		return std::make_unique<TSimulatorPair>();
	}
	if (method == "SFS") {
		logfile().startIndent("Simulating individuals from an SFS (parameter type=SFS):");
		return std::make_unique<TSimulatorSFS>(chs);
	}
	if (method == "HW") {
		logfile().startIndent("Simulating a individuals under Hardy-Weinberg (parameter type=HW):");
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

TSimulator::TSimulator(const std::string &method)
    : _outname(parameters().getParameterWithDefault<std::string>("out", "ATLAS_simulations")),
      _seqDepth(parameters().getParameterWithDefault("depth", 10.0)),
      _writeTrueGenotypes(parameters().parameterExists("writeTrueGenotypes")),
      _writeVariantInvariantBedFiles(parameters().parameterExists("writeVariantBED")), _reference(_outname + ".fasta"),
      _chromosomes(makeChromosomes()),
      _haploSimulator(makeHaploSimulator(method, _chromosomes)) {}

//---------------------------------------------------
// TBamSimulator
//---------------------------------------------------

TBAMSimulator::TBAMSimulator(const std::string& method) : TSimulator(method) {
	using genometools::Base;
	// depth
	logfile().list("Will simulate to an average depth of ", _seqDepth, ".");
	logfile().list("Will write output files with tag '" + _outname + "'.");
	_initializeReadSimulator();
}

//--------------------------------------------------------------
// Function to initialize read groups
//--------------------------------------------------------------
std::vector<std::string> TBAMSimulator::_readSimInfoPerReadGroup(const std::string &Filename, const std::string &Column,
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

void TBAMSimulator::_initializeDistribution(const std::string &ParameterName, const std::string &DefaultValue,
					    const std::string &Name,
					    std::function<void(TSimulatorSingleEndRead &, std::string)> function) {
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
}

void TBAMSimulator::_initializeQualityTransformations(const std::string &ParameterName, const std::string &Name) {
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

//--------------------------------------------------------------
// Initialize chromosomes, depth and base frequencies
//--------------------------------------------------------------

//--------------------------------------------------------------
// Run simulations
//--------------------------------------------------------------

void TBAMSimulator::_simulateReadsFromHaplotypes(const genometools::TChromosome &thisChr,
						 std::array<std::vector<Base>, 2> haplotypes,
						 TSimulatorBamFile &bamFile, std::string extraProgressText) {
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

void TBAMSimulator::runSimulations() {
	// open bam files
	TSimulatorBamFiles bamFiles(_haploSimulator->sampleSize(), _outname, _readGroups, _chromosomes);

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

		// now simulate and write reads
		logfile().startIndent("Simulating reads:");
		for (size_t i = 0; i < _haploSimulator->sampleSize(); ++i)
			_simulateReadsFromHaplotypes(chr, haplotypes.getHaplotypesOfIndividual(i), bamFiles[i],
			                             " for individual " + coretools::str::toString(i + 1));
		logfile().endIndent();

		// end of chromosome
		logfile().endIndent();
	}
	logfile().endIndent();
}

} // namespace Simulations
