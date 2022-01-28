/*
 * TSimulator.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
 */

#include "TSimulator.h"
#include <memory>
#include <numeric>

namespace Simulations {

//---------------------------------------------------
// TSimulator
//---------------------------------------------------

TSimulator::TSimulator(TLog *Logfile, TParameters &params, TRandomGenerator *RandomGenerator)
	: _logfile(Logfile), _randomGenerator(RandomGenerator),
	  _outname(params.getParameterWithDefault<std::string>("out", "ATLAS_simulations")),
	  _referenceDivergence(params.getParameterWithDefault("refDiv", 0.01)),
	  _seqDepth(params.getParameterWithDefault("depth", 10.0)),
	  _writeTrueGenotypes(params.parameterExists("writeTrueGenotypes")),
	  _writeVariantInvariantBedFiles(params.parameterExists("writeVariantBED")),
	  _referenceObj(_outname + ".fasta", _logfile) {
	// depth
	_logfile->list("Will simulate to an average depth of ", _seqDepth, ".");
	_logfile->list("Will simulate data with reference divergence = ", _referenceDivergence, ".");
	_logfile->list("Will write output files with tag '" + _outname + "'.");

	// fill cumul table for reference divergence
	_cumulRef[0] = 1.0 - _referenceDivergence;
	_cumulRef[1] = _cumulRef[0] + _referenceDivergence / 3.0;
	_cumulRef[2] = _cumulRef[1] + _referenceDivergence / 3.0;
	_cumulRef[3] = 1.0;

	// base frequencies
	std::vector<double> freq;
	coretools::str::fillContainerFromString(
		params.getParameterWithDefault<std::string>("baseFreq", "0.25,0.25,0.25,0.25"), freq, ',');
	if (freq.size() != 4) throw "baseFreq vector must have size = 4!";

	const auto sum            = coretools::containerSum(freq);
	_baseFreq[genometools::A] = freq[0] / sum;
	_baseFreq[genometools::C] = freq[1] / sum;
	_baseFreq[genometools::G] = freq[2] / sum;
	_baseFreq[genometools::T] = freq[3] / sum;

	_cumulBaseFreq[0] = _baseFreq[genometools::A];
	_cumulBaseFreq[1] = _cumulBaseFreq[0] + _baseFreq[genometools::C];
	_cumulBaseFreq[2] = _cumulBaseFreq[1] + _baseFreq[genometools::G];
	_cumulBaseFreq[3] = 1.0;

	_logfile->list("Simulating with base frequencies " + (std::string)_baseFreq);

	_initializeReadSimulator(params);
	_initializeChromosomes(params);
}

//--------------------------------------------------------------
// Function to initialize read groups
//--------------------------------------------------------------
std::vector<std::string> TSimulator::_readSimInfoPerReadGroup(const std::string &Filename, const std::string &Column,
							      const std::string &Name) {
	_logfile->listFlush("Reading " + Name + "(s) from file '" + Filename + "' ...");

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
	_logfile->done();
	_logfile->conclude("Read " + Name + "s for ", in.lineNumber(), " read groups.");

	// check if there was data for each read group
	for (size_t i = 0; i < found.size(); ++i) {
		if (!found[i])
			throw "No " + Name + " given for read group '" + _readGroups.getName(i) + "' in file '" + Filename + "'!";
	}
	return ret;
}

void TSimulator::_initializeReadGroup(const std::string &readLengthString, const BAM::TReadGroup &ReadGroup) {
	// single or paired end? Is indicated at beginning of readLengthString!
	if (readLengthString.find("single:") == 0) {
		_readSimulators.push_back(std::make_unique<TSimulatorSingleEndRead>(ReadGroup, _randomGenerator));
	} else if (readLengthString.find("paired:") == 0) {
		_readSimulators.push_back(std::make_unique<TSimulatorPairedEndReads>(ReadGroup, _randomGenerator));
	} else
		throw "Unable to understand string '" + readLengthString + "'!";

	// add read Length distribution
	const auto readLengthDist = coretools::str::readAfterLast(readLengthString, ':');
	_readSimulators.back()->setReadLengthDistribution(readLengthDist, _logfile);
}

void TSimulator::_initializeReadGroupsFromReadLengthDistribution(TParameters &params, const std::string &ParameterName,
								 const std::string &DefaultValue,
								 const std::string &Name) {
	_logfile->startIndent("Parsing read length distribution (parameter '" + ParameterName + "'):");
	const auto s = params.getParameterWithDefault<std::string>(ParameterName, DefaultValue);
	_readSimulators.clear();

	// We allow for two options:
	//   1) initialized from the command line (one for all read groups)
	//   2) read-group specific as given in a file

	// check if it is a file (should not contain a ':')
	const auto pos = s.find(":");
	if (pos != std::string::npos) {
		// Option 1: a single read length distribution for all
		//---------------------------------------------------------------------
		_logfile->list("Will use '" + s + " for all read groups.");

		// create read groups
		for (auto &r : _readGroups) { _initializeReadGroup(s, r); }
	} else {
		// Option 2: read group specific, given in a file
		//---------------------------------------------------------------------
		const std::vector<std::string> dist = _readSimInfoPerReadGroup(s, ParameterName, Name);

		for (size_t r = 0; r < dist.size(); ++r) { _initializeReadGroup(dist[r], _readGroups[r]); }
	}
	_logfile->endIndent();
}

void TSimulator::_initializeDistribution(TParameters &params, const std::string &ParameterName,
					 const std::string &DefaultValue, const std::string &Name,
					 std::function<void(TSimulatorSingleEndRead &, std::string)> function) {
	_logfile->startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");
	const auto s = params.getParameterWithDefault<std::string>(ParameterName, DefaultValue);

	// We allow for two options:
	//   1) initialized from the command line (one for all read groups)
	//   2) read-group specific as given in a file

	// check if it is a file (should not contain a '(')
	const auto pos = s.find("(");
	if (pos != std::string::npos) {
		// Option 1: a single read distribution for all
		//---------------------------------------------------------------------
		_logfile->list("Will use '" + s + " for all read groups.");

		// create read groups
		for (auto &r : _readSimulators) { function(*r, s); }
	} else {
		// Option 2: read group specific, given in a file
		//---------------------------------------------------------------------
		const std::vector<std::string> dist = _readSimInfoPerReadGroup(s, ParameterName, Name);

		for (uint32_t r = 0; r < _readSimulators.size(); ++r) { function(*_readSimulators[r], dist[r]); }
	}
	_logfile->endIndent();
}

void TSimulator::_initializePMD(TParameters &params, const std::string &ParameterName, const std::string &Name) {

	_logfile->startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");

	if (params.parameterExists(ParameterName)) {
		const auto pmdString = params.getParameter<std::string>(ParameterName);
		std::vector<uint16_t> ReadGroupsWithoutPMD;
		_PMD.initialize(pmdString, _readGroups, _logfile, ReadGroupsWithoutPMD);

		// add PMD to simulators
		for (size_t r = 0; r < _readSimulators.size(); ++r) { _readSimulators[r]->setPMD(&_PMD[r]); }
	} else {
		_logfile->list("Not simulating any PMD.");
	}
}

void TSimulator::_initializeQualityTransformations(TParameters &params, const std::string &ParameterName,
						   const std::string &Name) {
	_logfile->startIndent("Parsing " + Name + " (parameter " + ParameterName + "):");

	if (params.parameterExists(ParameterName)) {
		const auto recalString = params.getParameter<std::string>(ParameterName);
		_recal.initializeFromFile(recalString, _readGroups, _logfile);

		// add recal to simulators
		for (size_t r = 0; r < _readSimulators.size(); ++r) { _readSimulators[r]->setPMD(&_PMD[r]); }
	} else {
		_logfile->list("Not simulating any quality transformation.");
	}
	_logfile->endIndent();
}

void TSimulator::_addReadGroupsIfFile(const std::string &ParameterName, TParameters &Parameters,
				      BAM::TReadGroups &ReadGroups) {
	// check if parameter is given
	if (Parameters.parameterExists(ParameterName)) {
		const auto s = Parameters.getParameter<std::string>(ParameterName);

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

void TSimulator::_initializeReadSimulator(TParameters &params) {
	// For which read groups?
	// Check for each parameter if it is given per read group (a file) or common to all
	_readGroups.clear();
	_addReadGroupsIfFile("readLength", params, _readGroups);
	_addReadGroupsIfFile("qualityDist", params, _readGroups);
	_addReadGroupsIfFile("pmd", params, _readGroups);
	_addReadGroupsIfFile("recal", params, _readGroups);
	_addReadGroupsIfFile("readGroupFreq", params, _readGroups);

	// any read groups specified?
	if (_readGroups.empty()) {
		const auto numRG = params.getParameterWithDefault<int>("numReadGroups", 1);
		for (int i = 0; i < numRG; ++i) { _readGroups.add("SimReadGroup" + coretools::str::toString(i + 1)); }
		// report
		if (numRG == 1) {
			_logfile->startIndent("Initializing one read group (parameter 'numreadGroups'):");
		} else if (numRG > 1) {
			_logfile->startIndent("Initializing ", numRG, " identical read groups (parameter 'numreadGroups'):");
		} else {
			throw "numReadGroups must be at least 1!";
		}

	} else {
		_logfile->startIndent("Initializing ", _readGroups.size(), " individual read group(s):");
	}

	// A) read length: used to create simulator read groups
	//---------------
	_initializeReadGroupsFromReadLengthDistribution(params, "readLength", "single:fixed(100)", "read length");

	// B) initialize quality distribution
	//-----------------------------------
	_initializeDistribution(params, "baseQuality", "normal(30,10)[0,93]", "base quality distribution",
				&TSimulatorSingleEndRead::setQualityDistribution);

	// C) initialize mapping quality distribution
	//-----------------------------------
	_initializeDistribution(params, "mappingQuality", "normal(60,10)[1,255]", "mapping quality distribution",
				&TSimulatorSingleEndRead::setMappingQualityDistribution);

	// D) initialize PMD
	//------------------
	_initializePMD(params, "pmd", "post-mortem damage");

	// E) initialize quality transformation
	//-------------------------------------
	_initializeQualityTransformations(params, "recal", "recalibration models");

	// E) initialize contamination
	//----------------------------
	// TODO: Think about contamination object for both estimation and simulation

	// F) other things
	//----------------

	// initialize read group frequencies frequencies
	_initializeReadGroupFrequencies(params);
}

void TSimulator::_initializeReadGroupFrequencies(TParameters &params) {
	_cumulSimGroupFrequenies.reserve(_readSimulators.size());
	_simGroupFrequencies.reserve(_readSimulators.size());
	if (params.parameterExists("readGroupFreq")) {
		// read frequencies
		std::vector<std::string> vec;
		params.fillParameterIntoContainer("readGroupFreq", vec, true);
		std::vector<double> freq;
		coretools::str::repeatIndexes(vec, freq);
		if (freq.size() != _readSimulators.size())
			throw "Provided read group frequencies do not match number of read groups!";

		// normalize and print
		const auto sum = std::accumulate(freq.cbegin(), freq.cend(), 0.);

		_logfile->startIndent("Will simulate read groups with the following frequencies:");
		for (size_t i = 0; i < _readSimulators.size(); ++i) {
			_simGroupFrequencies[i] = freq[i] / sum;
			_logfile->list(_simGroupFrequencies[i], " " + _readSimulators[i]->name());
		}
		_logfile->endIndent();

		// fill cumulative
		_cumulSimGroupFrequenies[0] = _simGroupFrequencies[0];
		for (size_t i = 1; i < _readSimulators.size(); ++i)
			_cumulSimGroupFrequenies[i] = _cumulSimGroupFrequenies[i - 1] + _simGroupFrequencies[i];
		_cumulSimGroupFrequenies[_readSimulators.size() - 1] = 1.0; // ensure last entry is 1.0
	} else {
		// equal frequencies
		_logfile->list("Will simulate reads equally distributed among read groups.");
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
void TSimulator::_initializeChromosomes(TParameters &params) {
	_chromosomes.clear();
	std::vector<std::string> string_vec;
	params.fillParameterIntoContainerWithDefault("chrLength", string_vec, ',', {"1000000"});

	std::vector<uint32_t> chrLength;
	coretools::str::repeatIndexes(string_vec, chrLength);
	if (chrLength.size() < 1) throw "Issue understanding length of chromosomes!";

	std::vector<uint8_t> ploidy;
	if (params.parameterExists("ploidy")) {
		params.fillParameterIntoContainer("ploidy", string_vec, ',');
		coretools::str::repeatIndexes(string_vec, ploidy);
	} else {
		for (size_t i = 0; i < chrLength.size(); ++i) ploidy.push_back(2);
	}
	if (ploidy.size() != chrLength.size()) throw "List of chromosome lengths and ploidies differ in length!";
	for (auto &p : ploidy) {
		if (p != 1 && p != 2) { throw "Currently only ploidy 1 (haploid) or 2 (diploid) is supported!"; }
	}

	if (chrLength.size() == 1) {
		const auto numChr = params.getParameterWithDefault<int>("numChr", 1);
		std::string text  = "Will simulate " + coretools::str::toString(numChr);
		if (ploidy[0] == 1)
			text += " haploid";
		else
			text += " diploid";
		text += " chromosome(s) of length " + coretools::str::toString(chrLength[0]) + " each.";
		_logfile->list(text);
		for (int i = 0; i < numChr; ++i) {
			_chromosomes.appendChromosome("chr" + coretools::str::toString(i + 1), chrLength[0], ploidy[0]);
		}
	} else {
		_logfile->startIndent("Will simulate ", chrLength.size(), " chromosome(s) of the following length:");
		auto hIt = ploidy.begin();
		std::string text;
		for (auto it = chrLength.begin(); it != chrLength.end(); ++it, ++hIt) {
			text = coretools::str::toString(*it) + " (";
			if (*hIt == 1)
				text += "haploid)";
			else
				text += "diploid)";
			_logfile->list(text);
		}

		for (size_t i = 0; i < chrLength.size(); ++i) {
			_chromosomes.appendChromosome("chr" + coretools::str::toString(i + 1), chrLength[i], ploidy[i]);
		}
		_logfile->endIndent();
	}
}

//--------------------------------------------------------------
// Run simulations
//--------------------------------------------------------------
Base TSimulator::_sampleBase(const std::array<double, 4> &cumulProbs) {
	return static_cast<Base>(_randomGenerator->pickOne(cumulProbs));
}

Base TSimulator::_mutateBase(const Base &base, const std::array<double, 4> &cumulProbs) {
	return static_cast<Base>(base.get() + _randomGenerator->pickOne(cumulProbs));
}

void TSimulator::_simulateReadsFromHaplotypes(const BAM::TChromosome &thisChr, std::array<std::vector<Base>,2> haplotypes,
					      TSimulatorBamFile &bamFile, std::string extraProgressText) {
	// Initialize probabilities to simulate reads
	const uint64_t numReads = _averageReadLength == 0 ? 0 : thisChr.length * _seqDepth / _averageReadLength;

	const uint64_t chrLengthForStart = thisChr.length - _maxReadLength + 1;
	const double probReadPerSite     = 1.0 / chrLengthForStart;
	uint64_t numReadsSimulated       = 0;

	// initialize progress reporting
	coretools::TProgressReporter<uint64_t> reporter(
		_logfile, numReads, "Simulating about " + coretools::str::toString(numReads) + " reads" + extraProgressText);

	// now simulate
	for (uint32_t l = 0; l < chrLengthForStart; ++l) {
		// write unwritten alignments
		for (auto &rs : _readSimulators) rs->writeUnwrittenAlignments(l, bamFile);

		// draw random number to get number of reads starting at this position
		const auto numReadsHere = _randomGenerator->getBinomialRand(probReadPerSite, numReads);
		// now simulate
		if (numReadsHere > 0) {
			numReadsSimulated += numReadsHere;
			for (uint32_t r = 0; r < numReadsHere; ++r) {
				const auto rg = _randomGenerator->pickOne(_readSimulators.size(), _cumulSimGroupFrequenies.data());
				_readSimulators[rg]->simulate(haplotypes[_randomGenerator->sample(2)], thisChr.refID(), l, bamFile);
			}
			// report progress
			reporter.next();
		}
	}
	// write unwritten alignments
	for (auto &rs : _readSimulators) rs->writeUnwrittenAlignments(thisChr.length, bamFile);

	reporter.done();
	_logfile->conclude("Simulated a total of ", numReadsSimulated, " reads.");
}

void TSimulator::runSimulations() {
	// open bam files
	TSimulatorBamFiles bamFiles(_sampleSize, _outname, _readGroups, _chromosomes, _logfile);

	// prepare haplotypes and
	TSimulatorHaplotypes haplotypes(_sampleSize);

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
		_logfile->startIndent("Simulating chromosome " + chr.name + ":");

		// update reference storage and update haplotype lengths
		_referenceObj.setChr(chr.name, chr.length);
		haplotypes.setLength(chr.length);

		// simulate genotypes
		_logfile->listFlush("Simulating genotypes ...");
		if (chr.ploidy == 1)
			_simulateHaplotypesHaploid(haplotypes, chr);
		else
			_simulateHaplotypesDiploid(haplotypes, chr);
		_logfile->done();

		// write true genotypes
		if (_writeTrueGenotypes) {
			_logfile->listFlush("Writing true genotypes ...");
			haplotypes.writeTrueGenotypes(chr.name, _referenceObj);
			_logfile->done();
		}

		// write BED files
		if (_writeVariantInvariantBedFiles) bedFiles.write(haplotypes, chr.name);

		// now simulate and write reads
		_logfile->startIndent("Simulating reads:");
		for (int i = 0; i < _sampleSize; ++i)
			_simulateReadsFromHaplotypes(chr, haplotypes.getHaplotypesOfIndividual(i), bamFiles[i],
						     " for individual " + coretools::str::toString(i + 1));
		_logfile->endIndent();

		// end of chromosome
		_logfile->endIndent();
	}
	_logfile->endIndent();
}

//---------------------------------------------------------
// TSimulatorOneIndividual
//---------------------------------------------------------
TSimulatorOne::TSimulatorOne(TLog *Logfile, TParameters &params, TRandomGenerator *RandomGenerator)
	: TSimulator(Logfile, params, RandomGenerator) {
	_logfile->startIndent("Reading parameters to simulate a single individual:");

	_sampleSize = 1;

	// now theta
	std::vector<std::string> tmp;
	params.fillParameterIntoContainerWithDefault("theta", tmp, ',', {"0.001"});
	coretools::str::repeatIndexes(tmp, _thetas);
	if (_thetas.size() == 1) {
		_logfile->list("Will simulate a single individual with theta = ", _thetas[0], ".");
		for (unsigned int i = 1; i < _chromosomes.size(); ++i) _thetas.push_back(_thetas[0]);
	} else {
		_logfile->list("Will simulate a single individual with chromosome specific thetas " +
			       coretools::str::concatenateString(_thetas, ", "));
	}

	// one theta per chromosome
	if (_thetas.size() != _chromosomes.size())
		throw "Number of theta values provided does not match number of chromosomes to simulate!";

	// done
	_logfile->endIndent();
}

void TSimulatorOne::_simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes,
							 const BAM::TChromosome &chromosome) {
	// fill mutation table
	TSimulatorMutationtable mutTable(_baseFreq, _thetas[chromosome.refID()]);

	for (uint64_t l = 0; l < chromosome.length; ++l) {
		haplotypes(0, 0, l) = _sampleBase(_cumulBaseFreq);
		haplotypes(0, 1, l) = _sampleBase(mutTable[haplotypes(0, 0, l)]);

		// decide on reference sequence
		if (haplotypes(0, 0, l) == haplotypes(0, 1, l)) {
			_referenceObj[l] = _mutateBase(haplotypes(0, 0, l), _cumulRef);
		} else {
			_referenceObj[l] = static_cast<Base>(haplotypes(0, _randomGenerator->sample(2), l));
		}
	}
}

void TSimulatorOne::_simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes,
							 const BAM::TChromosome &chromosome) {
	// fill mutation table
	TSimulatorMutationtable mutTable(_baseFreq, _thetas[chromosome.refID()]);

	// now simulate genotypes
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		haplotypes(0, 0, l) = _sampleBase(_cumulBaseFreq);
		haplotypes(0, 1, l) = haplotypes(0, 0, l);

		// decide on ref
		_referenceObj[l] = _mutateBase(haplotypes(0, 0, l), _cumulRef);
	}
}

//---------------------------------------------------------
// TSimulatorPairOfIndividuals
//---------------------------------------------------------
TSimulatorPair::TSimulatorPair(TLog *Logfile, TParameters &params,
							 TRandomGenerator *RandomGenerator)
	: TSimulator(Logfile, params, RandomGenerator) {
	_logfile->startIndent("Reading parameters to simulate two individuals with a specific genetic distance:");

	_sampleSize = 2;

	// Initialize phis
	std::vector<std::string> tmp;
	params.fillParameterIntoContainer("phi", tmp, ',');
	coretools::str::repeatIndexes(tmp, _phis);

	if (_phis.size() != 9)
		throw "Wrong number of phi! Required are nine values for genotype combinations 00/00, 00/01, 01/00, 00/11, "
		      "01/01, 01/02, 00/12, 01/22, 01/23";

	// normalize phis
	const double sum = std::accumulate(_phis.cbegin(), _phis.cend(), 0.);
	if (sum != 1.0) {
		_logfile->list("Normalizing phi to sum to one (currently summing to ", sum, ").");
		for (auto & ph: _phis) ph /= sum;
	}
	_logfile->list("Used phi are: " + coretools::str::concatenateString(_phis, ", "));

	// initializes tables
	_fillTables();

	// done
	_logfile->endIndent();
}

void TSimulatorPair::_fillTables() {
	// file cumulative frequencies of cases (phis)
	double sum   = 0.0;
	for (size_t i = 0; i < _phis.size(); ++i) {
		sum += _phis[i];
		_cumulGenoCaseFrequencies[i] = sum;
	}
	_cumulGenoCaseFrequencies.back() = 1.0;
	if (fabs(sum - 1.0) > 0.0000000001)
		throw "Phis do not sum to 1.0! They sum to " + coretools::str::toString(sum) + ".";


	// case 0: aa/aa
	//-----------------------------------------
	sum = 0.0;
	for (size_t a = 0; a < 4; ++a) {
		sum += _baseFreq[static_cast<Base>(a)];
		_cumulGenoCombinationFreq[0].push_back(sum);
		_genoTrans[0].push_back({static_cast<Base>(a), static_cast<Base>(a), static_cast<Base>(a), static_cast<Base>(a)});
	}

	// cases 1 to 3: aa/ab, ab/aa, aa/bb
	//-----------------------------------------
	sum           = 0.0;
	for (size_t a = 0; a < 4; ++a) {
		for (size_t b = 0; b < 4; ++b) {
			if (a == b) continue;

			sum += _baseFreq[static_cast<Base>(a)] * _baseFreq[static_cast<Base>(b)];

			_cumulGenoCombinationFreq[1].push_back(sum);
			_cumulGenoCombinationFreq[2].push_back(sum);
			_cumulGenoCombinationFreq[3].push_back(sum);
			_genoTrans[1].push_back({static_cast<Base>(a), static_cast<Base>(a), static_cast<Base>(a), static_cast<Base>(b)});
			_genoTrans[2].push_back({static_cast<Base>(a), static_cast<Base>(b), static_cast<Base>(a), static_cast<Base>(a)});
			_genoTrans[3].push_back({static_cast<Base>(a), static_cast<Base>(a), static_cast<Base>(b), static_cast<Base>(b)});
		}
	}
	for (auto& c: _cumulGenoCombinationFreq[1]) c /= sum;
	for (auto& c: _cumulGenoCombinationFreq[2]) c /= sum;
	for (auto& c: _cumulGenoCombinationFreq[3]) c /= sum;

	// cases 4: ab/ab
	//-----------------------------------------
	sum          = 0.0;
	for (size_t a = 0; a < 3; ++a) {
		for (size_t b = a + 1; b < 4; ++b) {
			sum += _baseFreq[static_cast<Base>(a)] * _baseFreq[static_cast<Base>(b)];

			_cumulGenoCombinationFreq[4].push_back(sum);
			_genoTrans[4].push_back({static_cast<Base>(a), static_cast<Base>(b), static_cast<Base>(a), static_cast<Base>(b)});
		}
	}
	for (auto& c: _cumulGenoCombinationFreq[4]) c /= sum;

	// case 5: ab/ac
	//-----------------------------------------
	sum = 0.0;
	for (size_t a = 0; a < 4; ++a) {
		for (size_t b = 0; b < 4; ++b) {
			if (a == b) continue;
			for (size_t c = 0; c < 4; ++c) {
				if (c == a || c == b) continue;
				sum += _baseFreq[static_cast<Base>(a)] * _baseFreq[static_cast<Base>(b)] * _baseFreq[static_cast<Base>(c)];

				_cumulGenoCombinationFreq[5].push_back(sum);
				_genoTrans[5].push_back(
					{static_cast<Base>(a), static_cast<Base>(b), static_cast<Base>(a), static_cast<Base>(c)});
			}
		}
	}
	for (auto &c : _cumulGenoCombinationFreq[5]) c /= sum;

	// cases 6 and 7: aa/bc, ab/cc
	//-----------------------------------------
	sum   = 0.0;
	for (size_t a = 0; a < 4; ++a) {
		for (size_t b = 0; b < 4; ++b) {
			if (a == b) continue;
			for (size_t c = 0; c < 4; ++c) {
				if (c == a || c == b) continue;
				sum += _baseFreq[static_cast<Base>(a)] * _baseFreq[static_cast<Base>(b)] * _baseFreq[static_cast<Base>(c)];
				_cumulGenoCombinationFreq[6].push_back(sum);
				_cumulGenoCombinationFreq[7].push_back(sum);
				_genoTrans[6].push_back(
					{static_cast<Base>(a), static_cast<Base>(a), static_cast<Base>(b), static_cast<Base>(c)});
				_genoTrans[7].push_back(
					{static_cast<Base>(a), static_cast<Base>(b), static_cast<Base>(c), static_cast<Base>(c)});
			}
		}
	}
	for (auto &c : _cumulGenoCombinationFreq[6]) c /= sum;
	for (auto &c : _cumulGenoCombinationFreq[7]) c /= sum;

	// case 8: ab/cd
	//-----------------------------------------
	sum = 0.;
	for (size_t a = 0; a < 4; ++a) {
		for (size_t b = 0; b < 4; ++b) {
			if (a == b) continue;
			for (size_t c = 0; c < 4; ++c) {
				if (c == a || c == b) continue;
				for (uint8_t d = 0; d < 4; ++d) {
					if (d == a || d == b || d == c) continue;

					sum += 1./24;
					_cumulGenoCombinationFreq[8].push_back(sum);
				_genoTrans[8].push_back(
					{static_cast<Base>(a), static_cast<Base>(b), static_cast<Base>(c), static_cast<Base>(d)});
				}
			}
		}
	}
}

void TSimulatorPair::_simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes,
							     const BAM::TChromosome &chromosome) {
	// first run diploid
	_simulateHaplotypesDiploid(haplotypes, chromosome);

	// now set homozygous
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		// assign to haplotypes
		haplotypes(0, 1, l) = haplotypes(0, 0, l);
		haplotypes(1, 1, l) = haplotypes(1, 0, l);
	}
}

void TSimulatorPair::_simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes,
							     const BAM::TChromosome &chromosome) {
	// run across loci
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		// pick a case
		const int c = _randomGenerator->pickOne(_cumulGenoCaseFrequencies);

		// pick genotypes
		const int g = _randomGenerator->pickOne(_cumulGenoCombinationFreq[c]);

		// pick order
		const int o = _randomGenerator->sample(4);

		// assign to haplotypes
		haplotypes(0, 0, l) = _genoTrans[c][g][_orderLookup[o][0]];
		haplotypes(0, 1, l) = _genoTrans[c][g][_orderLookup[o][1]];
		haplotypes(1, 0, l) = _genoTrans[c][g][_orderLookup[o][2]];
		haplotypes(1, 1, l) = _genoTrans[c][g][_orderLookup[o][3]];

		// simulate reference
		if (c == 0) {
			_referenceObj[l].mutateWithStep(_randomGenerator->pickOne(_cumulRef));
		} else {
			const int r      = _randomGenerator->sample(4);
			_referenceObj[l] = _genoTrans[c][g][r];
		}
	}
}

//---------------------------------------------------------
// TSimulatorSFS
//---------------------------------------------------------
TSimulatorSFS::TSimulatorSFS(TLog *Logfile, TParameters &params, TRandomGenerator *RandomGenerator)
	: TSimulator(Logfile, params, RandomGenerator), _mutTable(_baseFreq) {
	_logfile->startIndent("Reading parameters to simulate a population sample given an SFS:");

	// sample size
	_sampleSize = params.getParameterWithDefault<int>("sampleSize", 10);

	// read SFS
	_logfile->startIndent("Initializing SFS:");
	if (params.parameterExists("sfs")) {
		_logfile->startIndent("Reading SFS from files:");

		std::vector<std::string> tmp;
		std::vector<std::string> sfsFileNames;
		params.fillParameterIntoContainer("sfs", tmp, ',');
		coretools::str::repeatIndexes(tmp, sfsFileNames);

		// if a single SFS is given: use it for all chromosomes
		if (sfsFileNames.size() == 1) {
			for (size_t _ = 1; _ < _chromosomes.size(); ++_) sfsFileNames.push_back(sfsFileNames.front());
		}

		// check if numbe rof chromosomes given matches number of chromosomes
		if (sfsFileNames.size() != _chromosomes.size())
			throw "Number of SFS files does not match number of chromosomes!";

		// initialize SFS from files
		const bool folded = params.parameterExists("folded");
		_initializeSFS(sfsFileNames, folded);
	} else if (params.parameterExists("theta")) {
		// parse theta from command line
		std::vector<std::string> tmp;
		params.fillParameterIntoContainer("theta", tmp, ',');
		std::vector<double> thetas;
		coretools::str::repeatIndexes(tmp, thetas);
		if (thetas.size() == 1) {
			_logfile->list("Will simulate from SFS with theta = ", thetas.front(), ".");
			for (unsigned int _ = 1; _ < _chromosomes.size(); ++_) thetas.push_back(thetas.front());
		} else {
			_logfile->list("Will simulate data from chromosome specific SFS with thetas " +
				       coretools::str::concatenateString(thetas, ", "));
		}
		_initializeSFS(thetas);
	} else
		throw "Either argument sfs or theta must be provided to simulate population samples!";

	// done
	_logfile->endIndent();
}

void TSimulatorSFS::_initializeSFS(const std::vector<double> &thetas) {
	if (thetas.size() != _chromosomes.size()) throw "Number of theta values does not match number of chromosomes!";

	// generate SFS for each chromosome
	_logfile->listFlush("Initializing SFS ...");
	for (size_t i = 0; i < _chromosomes.size(); ++i) {
		_sfs.push_back(std::make_unique<SFS>(_chromosomes[i].ploidy * _sampleSize, (float)thetas[i]));

		// save true SFS
		const auto filename = _outname + "_trueSFS_chr" + coretools::str::toString(_chromosomes[i].refID()) + ".txt";
		_sfs.back()->writeToFile(filename);
	}
	_logfile->done();
	_logfile->conclude("True SFS written to '" + _outname + "_trueSFS_chr*.txt'.");
}

void TSimulatorSFS::_initializeSFS(const std::vector<std::string> &sfsFileNames, bool folded) {
	if (sfsFileNames.size() != _chromosomes.size()) throw "Number of SFS files does not match number of chromosomes!";

	// read the SFS of each chromosome from the corresponding file
	for (size_t i = 0; i < _chromosomes.size(); ++i) {
		_logfile->listFlush("Reading the sfs of chromosome '" + _chromosomes[i].name + "' from file '" + sfsFileNames[i] + "' ...");
		if (folded)
			_sfs.push_back(std::make_unique<SFSfolded>(sfsFileNames[i]));
		else
			_sfs.push_back(std::make_unique<SFS>(sfsFileNames[i]));
		_logfile->done();

		const uint32_t nChr = _chromosomes[i].ploidy * _sampleSize;
		if (_sfs.back()->numChromosomes() != nChr) {
			throw coretools::str::toString("SFS does not match sample size! It contains data for ",
						       (*_sfs.rbegin())->numChromosomes(), " instead of ", nChr, " chromosomes.");
		}
	}
}

int is_odd(int x) { return x & 1; }

void TSimulatorSFS::_simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) {
	// now simulate haplotypes
	for (uint32_t l = 0; l < chromosome.length; ++l) {
		// pick alleles
		const Base ancestral = _sampleBase(_cumulBaseFreq);
		const Base derived   = _sampleBase(_mutTable[ancestral]);

		// pick derived allele frequency
		const int alleleCount = _sfs[chromosome.refID()]->getRandomAlleleCount(_randomGenerator);

		// pick haplotypes that are derived
		int numNeeded = alleleCount;
		for (int i = 0; i < _sampleSize; ++i) {
			if (_randomGenerator->getRand() < (double)numNeeded/(_sampleSize - i)) {
				haplotypes(i, 0, l) = derived;
				--numNeeded;
				if (numNeeded == 0) {
					for (int j = i + 1; j < _sampleSize; ++j) {
						haplotypes(i, 0, l) = ancestral;
						haplotypes(i, 1, l) = ancestral;
					}
					break;
				}
			} else
				haplotypes(i, 0, l) = ancestral;
			// make homozygous
			haplotypes(i, 1, l) = haplotypes(i, 0, l);
		}

		// decide on reference sequence
		if (alleleCount > 0) {
			if (_randomGenerator->getRand() < (double)alleleCount/_sampleSize)
				_referenceObj[l] = derived;
			else
				_referenceObj[l] = ancestral;
		} else {
			_referenceObj[l] = _mutateBase(ancestral, _cumulRef);
		}
	}
}

void TSimulatorSFS::_simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes, const BAM::TChromosome &chromosome) {
	const int numHaplotypes = 2 * _sampleSize;
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		// pick alleles
		const Base ancestral = _sampleBase(_cumulBaseFreq);
		const Base derived   = _sampleBase(_mutTable[ancestral]);

		// pick derived allele frequency
		const int alleleCount = _sfs[chromosome.refID()]->getRandomAlleleCount(_randomGenerator);
		// oo << alleleCount << "\n";

		// pick haplotypes that are derived
		if (alleleCount == 0) {
			for (int i = 0; i < _sampleSize; ++i) {
				haplotypes(i, 0, l) = ancestral;
				haplotypes(i, 1, l) = ancestral;
			}
			// decide on reference sequence
			_referenceObj[l] = _mutateBase(ancestral, _cumulRef);
		} else {
			int numNeeded = alleleCount;
			for (int i = 0; i < numHaplotypes; ++i) {
				double prob = (double)numNeeded/(numHaplotypes - i);
				if (_randomGenerator->getRand() < prob) {
					haplotypes(i / 2, is_odd(i), l) = derived;
					--numNeeded;
					if (numNeeded == 0) {
						for (int j = i + 1; j < numHaplotypes; ++j) haplotypes(j / 2, is_odd(j), l) = ancestral;
						break;
					}
				} else
					haplotypes(i / 2, is_odd(i), l) = ancestral;
			}

			// decide on reference sequence
			if (_randomGenerator->getRand() < (double)alleleCount / (double)numHaplotypes)
				_referenceObj[l] = derived;
			else
				_referenceObj[l] = ancestral;
		}
	}
}

//---------------------------------------------------------
// TSimulatorHardyWeinberg
//---------------------------------------------------------
TSimulatorHW::TSimulatorHW(coretools::TLog *Logfile, TParameters &params,
						 coretools::TRandomGenerator *RandomGenerator)
	: TSimulator(Logfile, params, RandomGenerator), _fracPoly(params.getParameterWithDefault("fracPoly", 0.1)),
	  _alpha(params.getParameterWithDefault("alpha", 0.5)), _beta(params.getParameterWithDefault("beta", 0.5)),
	  _F(params.getParameterWithDefault("F", 0.0)), _mutTable(_baseFreq) {
	_logfile->startIndent("Reading parameters to simulate a population sample under Hardy-Weinberg equilibrium:");

	// sample size
	_sampleSize = params.getParameterWithDefault<int>("sampleSize", 10);

	// parameters of beta distribution
	_logfile->list("Will simulate ", _fracPoly, " of all sites as polymorphic. (parameter fracPoly)");
	if (_alpha <= 0.0) throw "Alpha must be > 0!";
	if (_beta <= 0.0) throw "Beta must be > 0!";
	_logfile->list("Polymoprhic sites will have derived allele frequencies f~Beta(", _alpha, ", ", _beta,
		       "). (parameters alpha, beta)");
	if (_F == 0.0) {
		_logfile->list("Will assume no inbreeding. (parameter F=0)");
	} else {
		_logfile->list("Will use an inbreeding coefficient of ", _F, ". (parameter F)");
		if (_F < 0.0 || _F > 1.0) throw "Inbreeding coefficient F must be within [0,1]!";
	}

	// write true allele freq?
	if (params.parameterExists("writeTrueAlleleFreq")) {
		const auto alleleFreqFile = _outname + "_trueAlleleFreq.txt.gz";
		_logfile->list("Will write true allele frequencies to file '" + alleleFreqFile + "'.");
		_trueFreqFile.open(alleleFreqFile);
		_trueFreqFile.writeHeader({"Chr", "Pos", "Ancestral", "Derived", "derivedFreq", "MAF"});
		_writeTrueAlleleFreq = true;
	}

	// done
	_logfile->endIndent();
}

void TSimulatorHW::_fillCumulGenoProb(double f) {
	const double oneMinus_f = 1.0 - f;
	_cumulGenoProb[0]        = _F * oneMinus_f + (1.0 - _F) * oneMinus_f * oneMinus_f;
	_cumulGenoProb[1]        = _cumulGenoProb[0] + (1.0 - _F) * 2.0 * f * oneMinus_f;
	_cumulGenoProb[2]        = 1.0;
}

void TSimulatorHW::_simulateSite(TSimulatorHWSite &site, const std::string &chr, uint64_t pos) {
	// simulate bases
	site.reference   = static_cast<Base>(_randomGenerator->pickOne(_cumulBaseFreq));
	site.alternative = static_cast<Base>(_randomGenerator->pickOne(_mutTable[site.reference.get()]));

	// is the site polymorphic?
	if (_randomGenerator->getRand() < _fracPoly) {
		site.isPolymorphic = true;
		site.f             = _randomGenerator->getBetaRandom(_alpha, _beta);

		// reference is a random sample from pop with frequency f: flip if ref is alt!
		if (_randomGenerator->getRand() < site.f) {
			std::swap(site.reference, site.alternative);
			site.f = 1 - site.f;
		}
	} else {
		site.isPolymorphic = false;
		// is reference diverged?
		site.f = _randomGenerator->getRand() < _referenceDivergence ? 1. : 0.;
	}
	// store reference
	_referenceObj[pos] = site.reference;

	// write true frequency: pos is 1 based!
	if (_writeTrueAlleleFreq) {
		_trueFreqFile << chr << pos + 1 << site.reference << site.alternative << site.f;
		_trueFreqFile << (site.f < 0.5 ? site.f : 1. - site.f) << std::endl;
	}
}

void TSimulatorHW::_fillhaplotypesMonomoprhic(TSimulatorHaplotypes &haplotypes, uint64_t locus,
							 const TSimulatorHWSite &site) {
	if (site.f == 0.0) {
		for (int i = 0; i < _sampleSize; ++i) {
			haplotypes(i, 0, locus) = site.reference;
			haplotypes(i, 1, locus) = site.reference;
		}
	} else {
		for (int i = 0; i < _sampleSize; ++i) {
			haplotypes(i, 0, locus) = site.alternative;
			haplotypes(i, 1, locus) = site.alternative;
		}
	}
}

void TSimulatorHW::_simulateHaplotypesHaploid(TSimulatorHaplotypes &haplotypes,
							 const BAM::TChromosome &chromosome) {
	// storage

	// now simulate haplotypes
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		TSimulatorHWSite site;
		// simulate site
		_simulateSite(site, chromosome.name, l);

		// polymoprhic or not?
		if (site.isPolymorphic) {
			// simulate genotypes
			for (int i = 0; i < _sampleSize; ++i) {
				if (_randomGenerator->getRand() < site.f) {
					haplotypes(i, 0, l) = site.alternative;
					haplotypes(i, 1, l) = site.alternative;
				} else {
					haplotypes(i, 0, l) = site.reference;
					haplotypes(i, 1, l) = site.reference;
				}
			}
		} else {
			_fillhaplotypesMonomoprhic(haplotypes, l, site);
		}
	}
}

void TSimulatorHW::_simulateHaplotypesDiploid(TSimulatorHaplotypes &haplotypes,
							 const BAM::TChromosome &chromosome) {
	// storage

	// now simulate haplotypes
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		TSimulatorHWSite site;
		// simulate site
		_simulateSite(site, chromosome.name, l);

		// polymoprhic or not?
		if (site.isPolymorphic) {
			_fillCumulGenoProb(site.f);

			// simulate genotypes
			for (int i = 0; i < _sampleSize; ++i) {
				int geno = _randomGenerator->pickOne(3, _cumulGenoProb);
				if (geno == 0) {
					haplotypes(i, 0, l) = site.reference;
					haplotypes(i, 1, l) = site.reference;
				} else if (geno == 1) {
					if (_randomGenerator->getRand() < 0.5) {
						haplotypes(i, 0, l) = site.reference;
						haplotypes(i, 1, l) = site.alternative;
					} else {
						haplotypes(i, 0, l) = site.alternative;
						haplotypes(i, 1, l) = site.reference;
					}
				} else {
					haplotypes(i, 0, l) = site.alternative;
					haplotypes(i, 1, l) = site.alternative;
				}
			}
		} else {
			_fillhaplotypesMonomoprhic(haplotypes, l, site);
		}
	}
}

//--------------------------------------------------------------------
// Functions to simulate pooled data
//--------------------------------------------------------------------
// TODO: Need to switch to haplotype model

/*
void TSimulator::simulatePooledData(int sampleSize, SFS & sfs, std::string outname){
    //open BAM file
    openBamFile(outname + ".bam");

    //open FASTA file for reference sequences
    std::string filename = outname + ".fasta";
    openFastaFile(filename);

    //prepare variables
    float* altFreq = NULL;
    long numReads;
    long chrLengthForStart;
    double probReadPerSite;
    int numReadsHere;
    long numReadsSimulated;
    initializeQualToErrorTable();

    //open frequency file
    filename = outname + "_frequencies.txt";
    std::ofstream freqFile(filename.c_str());

    //simulate sequences
    int refId = 0;
    for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++refId){
        logfile->startIndent("Simulating chromosome " + chrIt->name + ":");

        //simulate reference and alternative sequence
        simulateReferenceAndAlternativeSequenceCurChromosome();

        //simulate alternative frequencies (and write to file)
        logfile->listFlush("Simulating alternative allele frequencies ...");
        delete[] altFreq;
        altFreq = new float[chrIt->length];
        for(int l=0; l<chrIt->length; ++l){
            altFreq[l] = sfs.getRandomFrequency(randomGenerator);
            freqFile << chrIt->name << "\t" << l+1 << altFreq[l] << "\n";
        }
        logfile->done();

        //simulating reads
        numReads = chrIt->length * seqDepth / readLength;
        chrLengthForStart = chrIt->length - readLength;
        probReadPerSite = 1.0 / (double) chrLengthForStart;
        numReadsSimulated = 0;
        bamAlignment.RefID = refId;
        int prog;
        int oldProg = 0;
        std::string progressString = "Simulating about " + toString(numReads) + " reads ...";
        logfile->listFlush(progressString);
        for(long l=0; l<chrLengthForStart; ++l){
            //draw random number to get number of reads starting at this position
            numReadsHere = randomGenerator->getBiomialRand(probReadPerSite, numReads);

            //now simulate
            if(numReadsHere > 0){
                simulateReads(numReadsHere, l, altFreq);
                numReadsSimulated += numReadsHere;

                //report progress
                prog = 100.0 * (float) numReadsSimulated / (float) numReads;
                if(prog > oldProg){
                    oldProg = prog;
                    logfile->listOverFlush(progressString + "(" + toString(prog) + "%)");
                }
            }
        }
        logfile->overList(progressString + " done!  ");
        logfile->conclude("Simulated a total of " + toString(numReadsSimulated) + " reads.");
        logfile->endIndent();
    }

    //close stuff
    closeBamFile();
    closeFastaFile();
    freqFile.close();

    //clear memory
    delete[] altFreq;
}
*/

} // namespace Simulations
