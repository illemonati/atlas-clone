/*
 * THaplotypeSimulator.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: phaentu
*/


#include "THaplotypeSimulator.h"

#include <math.h>
#include <algorithm>
#include <exception>
#include <numeric>
#include <ostream>
#include <utility>

#include "TChromosomes.h"
#include "TGenotypeData.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TSimulatorAuxiliaryTools.h"
#include "probability.h"
#include "stringFunctions.h"
#include "weakTypes.h"

namespace Simulations {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::instances::randomGenerator;
using genometools::Base;
using GenotypeLikelihoods::TBaseProbabilities;

namespace impl {
std::string toString(const TBaseProbabilities &probs) {
	std::string s = "";
	for (auto b = Base::min; b < Base::max; ++b) {
		s += genometools::toString(b) + ": " + coretools::str::toString(probs[b]) + ", ";
	}
	return s.substr(0, s.size() - 1);
}
} // namespace impl

Base sampleBase(const coretools::TStrongArray<double, Base> &cumulProbs) {
	return genometools::Base(randomGenerator().pickOne(cumulProbs));
}

Base mutateBase(Base base, const coretools::TStrongArray<double, Base> &cumulProbs) {
	using namespace genometools;
	return Base((coretools::index(base) + randomGenerator().pickOne(cumulProbs)) % coretools::index(Base::max));
}

THaplotypeSimulator::THaplotypeSimulator()
    : _referenceDivergence(parameters().getParameterWithDefault("refDiv", 0.01)) {
	logfile().list("Will simulate data with reference divergence = ", _referenceDivergence, ".");
	_cumulRef[Base::A] = 1.0 - _referenceDivergence;
	_cumulRef[Base::C] = _cumulRef[Base::A] + _referenceDivergence / 3.0;
	_cumulRef[Base::G] = _cumulRef[Base::C] + _referenceDivergence / 3.0;
	_cumulRef[Base::T] = 1.0;

	// base frequencies
	std::vector<double> freq;
	coretools::str::fillContainerFromString(
	    parameters().getParameterWithDefault<std::string>("baseFreq", "0.25,0.25,0.25,0.25"), freq, ',');
	if (freq.size() != 4) throw "baseFreq vector must have size = 4!";

	_baseFreq = GenotypeLikelihoods::TBaseProbabilities{freq};

	logfile().list("Simulating with base frequencies " + impl::toString(_baseFreq));

	_cumulBaseFreq[Base::A] = _baseFreq[Base::A];
	_cumulBaseFreq[Base::C] = _cumulBaseFreq[Base::A] + _baseFreq[Base::C];
	_cumulBaseFreq[Base::G] = _cumulBaseFreq[Base::C] + _baseFreq[Base::G];
	_cumulBaseFreq[Base::T] = 1.0;
}

//---------------------------------------------------------
// TSimulatorOneIndividual
//---------------------------------------------------------
TSimulatorOne::TSimulatorOne(size_t nChoromosomes) : THaplotypeSimulator() {
	logfile().startIndent("Reading parameters to simulate a single individual:");

	// now theta
	std::vector<std::string> tmp;
	parameters().fillParameterIntoContainerWithDefault("theta", tmp, ',', {"0.001"});
	coretools::str::repeatIndexes(tmp, _thetas);
	if (_thetas.size() == 1) {
		logfile().list("Will simulate a single individual with theta = ", _thetas[0], ".");
		for (unsigned int i = 1; i < nChoromosomes; ++i) _thetas.push_back(_thetas[0]);
	} else {
		logfile().list("Will simulate a single individual with chromosome specific thetas " +
			       coretools::str::concatenateString(_thetas, ", "));
	}

	// one theta per chromosome
	if (_thetas.size() != nChoromosomes)
		throw "Number of theta values provided does not match number of chromosomes to simulate!";

	// done
	logfile().endIndent();
}

	void TSimulatorOne::simulateDiploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					       const genometools::TChromosome &chromosome) {
	// fill mutation table
	TSimulatorMutationtable mutTable(_baseFreq, _thetas[chromosome.refID()]);

	for (uint64_t l = 0; l < chromosome.length; ++l) {
		haplotypes(0, 0, l) = sampleBase(_cumulBaseFreq);
		haplotypes(0, 1, l) = sampleBase(mutTable[haplotypes(0, 0, l)]);

		// decide on reference sequence
		if (haplotypes(0, 0, l) == haplotypes(0, 1, l)) {
			reference[l] = mutateBase(haplotypes(0, 0, l), _cumulRef);
		} else {
			reference[l] = haplotypes(0, randomGenerator().sample(2), l);
		}
	}
}

void TSimulatorOne::simulateHaploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					       const genometools::TChromosome &chromosome) {
	// fill mutation table
	TSimulatorMutationtable mutTable(_baseFreq, _thetas[chromosome.refID()]);

	// now simulate genotypes
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		haplotypes(0, 0, l) = sampleBase(_cumulBaseFreq);
		haplotypes(0, 1, l) = haplotypes(0, 0, l);

		// decide on ref
		reference[l] = mutateBase(haplotypes(0, 0, l), _cumulRef);
	}
}

//---------------------------------------------------------
// TSimulatorPairOfIndividuals
//---------------------------------------------------------
TSimulatorPair::TSimulatorPair() : THaplotypeSimulator() {
	logfile().startIndent("Reading parameters to simulate two individuals with a specific genetic distance:");

	// Initialize phis
	std::vector<std::string> tmp;
	parameters().fillParameterIntoContainer("phi", tmp, ',');
	coretools::str::repeatIndexes(tmp, _phis);

	if (_phis.size() != 9)
		throw "Wrong number of phi! Required are nine values for genotype combinations 00/00, 00/01, 01/00, 00/11, "
		      "01/01, 01/02, 00/12, 01/22, 01/23";

	// normalize phis
	const double sum = std::accumulate(_phis.cbegin(), _phis.cend(), 0.);
	if (sum != 1.0) {
		logfile().list("Normalizing phi to sum to one (currently summing to ", sum, ").");
		for (auto &ph : _phis) ph /= sum;
	}
	logfile().list("Used phi are: " + coretools::str::concatenateString(_phis, ", "));

	// initializes tables
	_fillTables();

	// done
	logfile().endIndent();
}

void TSimulatorPair::_fillTables() {
	// file cumulative frequencies of cases (phis)
	double sum = 0.0;
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
	for (Base a = Base::min; a < Base::max; ++a) {
		sum += _baseFreq[a];
		_cumulGenoCombinationFreq[0].push_back(sum);
		_genoTrans[0].push_back({a, a, a, a});
	}

	// cases 1 to 3: aa/ab, ab/aa, aa/bb
	//-----------------------------------------
	sum = 0.0;
	for (Base a = Base::min; a < Base::max; ++a) {
		for (Base b = Base::min; b < Base::max; ++b) {
			if (a == b) continue;

			sum += _baseFreq[a] * _baseFreq[b];

			_cumulGenoCombinationFreq[1].push_back(sum);
			_cumulGenoCombinationFreq[2].push_back(sum);
			_cumulGenoCombinationFreq[3].push_back(sum);
			_genoTrans[1].push_back({a, a, a, b});
			_genoTrans[2].push_back({a, b, a, a});
			_genoTrans[3].push_back({a, a, b, b});
		}
	}
	for (auto &c : _cumulGenoCombinationFreq[1]) c /= sum;
	for (auto &c : _cumulGenoCombinationFreq[2]) c /= sum;
	for (auto &c : _cumulGenoCombinationFreq[3]) c /= sum;

	// cases 4: ab/ab
	//-----------------------------------------
	sum = 0.0;
	for (Base a = Base::min; a < genometools::Base::T; ++a) {
		for (Base b = coretools::next(a); b < Base::max; ++b) {
			sum += _baseFreq[a] * _baseFreq[b];

			_cumulGenoCombinationFreq[4].push_back(sum);
			_genoTrans[4].push_back({a, b, a, b});
		}
	}
	for (auto &c : _cumulGenoCombinationFreq[4]) c /= sum;

	// case 5: ab/ac
	//-----------------------------------------
	sum = 0.0;
	for (Base a = Base::min; a < Base::max; ++a) {
		for (Base b = Base::min; b < Base::max; ++b) {
			if (a == b) continue;
			for (Base c = Base::min; c < Base::max; ++c) {
				if (c == a || c == b) continue;
				sum += _baseFreq[a] * _baseFreq[b] * _baseFreq[c];

				_cumulGenoCombinationFreq[5].push_back(sum);
				_genoTrans[5].push_back({a, b, a, c});
			}
		}
	}
	for (auto &c : _cumulGenoCombinationFreq[5]) c /= sum;

	// cases 6 and 7: aa/bc, ab/cc
	//-----------------------------------------
	sum = 0.0;
	for (Base a = Base::min; a < Base::max; ++a) {
		for (Base b = Base::min; b < Base::max; ++b) {
			if (a == b) continue;
			for (Base c = Base::min; c < Base::max; ++c) {
				if (c == a || c == b) continue;
				sum += _baseFreq[a] * _baseFreq[b] * _baseFreq[c];
				_cumulGenoCombinationFreq[6].push_back(sum);
				_cumulGenoCombinationFreq[7].push_back(sum);
				_genoTrans[6].push_back({a, a, b, c});
				_genoTrans[7].push_back({a, b, c, c});
			}
		}
	}
	for (auto &c : _cumulGenoCombinationFreq[6]) c /= sum;
	for (auto &c : _cumulGenoCombinationFreq[7]) c /= sum;

	// case 8: ab/cd
	//-----------------------------------------
	sum = 0.;
	for (Base a = Base::min; a < Base::max; ++a) {
		for (Base b = Base::min; b < Base::max; ++b) {
			if (a == b) continue;
			for (Base c = Base::min; c < Base::max; ++c) {
				if (c == a || c == b) continue;
				for (Base d = Base::min; d < Base::max; ++d) {
					if (d == a || d == b || d == c) continue;

					sum += 1. / 24;
					_cumulGenoCombinationFreq[8].push_back(sum);
					_genoTrans[8].push_back({a, b, c, d});
				}
			}
		}
	}
}

void TSimulatorPair::simulateHaploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
						const genometools::TChromosome &chromosome) {
	// first run diploid
	simulateDiploid(haplotypes, reference, chromosome);

	// now set homozygous
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		// assign to haplotypes
		haplotypes(0, 1, l) = haplotypes(0, 0, l);
		haplotypes(1, 1, l) = haplotypes(1, 0, l);
	}
}

void TSimulatorPair::simulateDiploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
						const genometools::TChromosome &chromosome) {
	// run across loci
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		// pick a case
		const int c = randomGenerator().pickOne(_cumulGenoCaseFrequencies);

		// pick genotypes
		const int g = randomGenerator().pickOne(_cumulGenoCombinationFreq[c]);

		// pick order
		const int o = randomGenerator().sample(4);

		// assign to haplotypes
		haplotypes(0, 0, l) = _genoTrans[c][g][_orderLookup[o][0]];
		haplotypes(0, 1, l) = _genoTrans[c][g][_orderLookup[o][1]];
		haplotypes(1, 0, l) = _genoTrans[c][g][_orderLookup[o][2]];
		haplotypes(1, 1, l) = _genoTrans[c][g][_orderLookup[o][3]];

		// simulate reference
		if (c == 0) {
			reference[l] = mutateBase(reference[l], _cumulRef);
		} else {
			const int r   = randomGenerator().sample(4);
			reference[l] = _genoTrans[c][g][r];
		}
	}
}

//---------------------------------------------------------
// TSimulatorSFS
//---------------------------------------------------------
	TSimulatorSFS::TSimulatorSFS(const genometools::TChromosomes& chromosomes) : THaplotypeSimulator(), _mutTable(_baseFreq) {
	logfile().startIndent("Reading parameters to simulate a population sample given an SFS:");
	const auto nChromosomes = chromosomes.size();

	// sample size
	_sampleSize = parameters().getParameterWithDefault<int>("sampleSize", 10);

	// read SFS
	logfile().startIndent("Initializing SFS:");
	if (parameters().parameterExists("sfs")) {
		logfile().startIndent("Reading SFS from files:");

		std::vector<std::string> tmp;
		std::vector<std::string> sfsFileNames;
		parameters().fillParameterIntoContainer("sfs", tmp, ',');
		coretools::str::repeatIndexes(tmp, sfsFileNames);

		// if a single SFS is given: use it for all chromosomes
		if (sfsFileNames.size() == 1) {
			for (size_t _ = 1; _ < nChromosomes; ++_) sfsFileNames.push_back(sfsFileNames.front());
		}

		// check if numbe rof chromosomes given matches number of chromosomes
		if (sfsFileNames.size() != nChromosomes)
			throw "Number of SFS files does not match number of chromosomes!";

		// initialize SFS from files
		const bool folded = parameters().parameterExists("folded");
		_initializeSFS(chromosomes, sfsFileNames, folded);
	} else if (parameters().parameterExists("theta")) {
		// parse theta from command line
		std::vector<std::string> tmp;
		parameters().fillParameterIntoContainer("theta", tmp, ',');
		std::vector<double> thetas;
		coretools::str::repeatIndexes(tmp, thetas);
		if (thetas.size() == 1) {
			logfile().list("Will simulate from SFS with theta = ", thetas.front(), ".");
			for (unsigned int _ = 1; _ < nChromosomes; ++_) thetas.push_back(thetas.front());
		} else {
			logfile().list("Will simulate data from chromosome specific SFS with thetas " +
				       coretools::str::concatenateString(thetas, ", "));
		}
		const bool folded = parameters().parameterExists("folded");
		_initializeSFS(chromosomes, thetas, folded);
	} else
		throw "Either argument sfs or theta must be provided to simulate population samples!";

	// done
	logfile().endIndent();
}

void TSimulatorSFS::_initializeSFS(const genometools::TChromosomes& chromosomes, const std::vector<double> &thetas, bool folded) {
	if (thetas.size() != chromosomes.size()) throw "Number of theta values does not match number of chromosomes!";
	const auto outname = parameters().getParameterWithDefault<std::string>("out", "ATLAS_simulations");

	// generate SFS for each chromosome
	logfile().listFlush("Initializing SFS ...");
	for (size_t i = 0; i < chromosomes.size(); ++i) {
		if (folded){
			_sfs.push_back(std::make_unique<SFSfolded>(_sampleSize, (float)thetas[i]));
			/*// save true SFS
			const auto filename = outname + "_trueSFS_chr" + coretools::str::toString(chromosomes[i].refID() + 1) + ".txt";
			_sfs.back()->writeToFile(filename);*/
		} else
			_sfs.push_back(std::make_unique<SFS>(chromosomes[i].ploidy * _sampleSize, (float)thetas[i]));

		// save true SFS}
		const auto filename = outname + "_trueSFS_chr" + coretools::str::toString(chromosomes[i].refID() + 1) + ".txt";
		_sfs.back()->writeToFile(filename);

	}
	logfile().done();
	logfile().conclude("True SFS written to '" + outname + "_trueSFS_chr*.txt'.");
}

void TSimulatorSFS::_initializeSFS(const genometools::TChromosomes& chromosomes, const std::vector<std::string> &sfsFileNames, bool folded) {
	if (sfsFileNames.size() != chromosomes.size()) throw "Number of SFS files does not match number of chromosomes!";

	// read the SFS of each chromosome from the corresponding file
	for (size_t i = 0; i < chromosomes.size(); ++i) {
		logfile().listFlush("Reading the sfs of chromosome '" + chromosomes[i].name + "' from file '" +
				    sfsFileNames[i] + "' ...");
		if (folded)
			_sfs.push_back(std::make_unique<SFSfolded>(sfsFileNames[i]));
		else
			_sfs.push_back(std::make_unique<SFS>(sfsFileNames[i]));
		logfile().done();

		const uint32_t nChr = chromosomes[i].ploidy * _sampleSize;
		if (_sfs.back()->numChromosomes() != nChr) {
			throw coretools::str::toString("SFS does not match sample size! It contains data for ",
						       (*_sfs.rbegin())->numChromosomes(), " instead of ", nChr, " chromosomes.");
		}
	}
}

constexpr int is_odd(int x) noexcept { return x & 1; }

void TSimulatorSFS::simulateHaploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					       const genometools::TChromosome &chromosome) {
	// now simulate haplotypes
	for (uint32_t l = 0; l < chromosome.length; ++l) {
		// pick alleles
		const Base ancestral = sampleBase(_cumulBaseFreq);
		const Base derived   = sampleBase(_mutTable[ancestral]);

		// pick derived allele frequency
		const int alleleCount = _sfs[chromosome.refID()]->getRandomAlleleCount();

		// pick haplotypes that are derived
		int numNeeded = alleleCount;
		for (int i = 0; i < _sampleSize; ++i) {
			if (randomGenerator().getRand() < (double)numNeeded / (_sampleSize - i)) {
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
			if (randomGenerator().getRand() < (double)alleleCount / _sampleSize)
				reference[l] = derived;
			else
				reference[l] = ancestral;
		} else {
			reference[l] = mutateBase(ancestral, _cumulRef);
		}
	}
}

void TSimulatorSFS::simulateDiploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					       const genometools::TChromosome &chromosome) {
	const int numHaplotypes = 2 * _sampleSize;
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		// pick alleles
		const Base ancestral = sampleBase(_cumulBaseFreq);
		const Base derived   = sampleBase(_mutTable[ancestral]);

		// pick derived allele frequency
		const int alleleCount = _sfs[chromosome.refID()]->getRandomAlleleCount();
		// oo << alleleCount << "\n";

		// pick haplotypes that are derived
		if (alleleCount == 0) {
			for (int i = 0; i < _sampleSize; ++i) {
				haplotypes(i, 0, l) = ancestral;
				haplotypes(i, 1, l) = ancestral;
			}
			// decide on reference sequence
			reference[l] = mutateBase(ancestral, _cumulRef);
		} else {
			int numNeeded = alleleCount;
			for (int i = 0; i < numHaplotypes; ++i) {
				double prob = (double)numNeeded / (numHaplotypes - i);
				if (randomGenerator().getRand() < prob) {
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
			if (randomGenerator().getRand() < (double)alleleCount / (double)numHaplotypes)
				reference[l] = derived;
			else
				reference[l] = ancestral;
		}
	}
}

//---------------------------------------------------------
// TSimulatorHardyWeinberg
//---------------------------------------------------------
TSimulatorHW::TSimulatorHW()
    : THaplotypeSimulator(), _fracPoly(parameters().getParameterWithDefault("fracPoly", 0.1)),
      _alpha(parameters().getParameterWithDefault("alpha", 0.5)),
      _beta(parameters().getParameterWithDefault("beta", 0.5)), _F(parameters().getParameterWithDefault("F", 0.0)),
      _mutTable(_baseFreq) {
	logfile().startIndent("Reading parameters to simulate a population sample under Hardy-Weinberg equilibrium:");

	// sample size
	_sampleSize = parameters().getParameterWithDefault<int>("sampleSize", 10);

	// parameters of beta distribution
	logfile().list("Will simulate ", _fracPoly, " of all sites as polymorphic. (parameter fracPoly)");
	if (_alpha <= 0.0) throw "Alpha must be > 0!";
	if (_beta <= 0.0) throw "Beta must be > 0!";
	logfile().list("Polymoprhic sites will have derived allele frequencies f~Beta(", _alpha, ", ", _beta,
		       "). (parameters alpha, beta)");
	if (_F == 0.0) {
		logfile().list("Will assume no inbreeding. (parameter F=0)");
	} else {
		logfile().list("Will use an inbreeding coefficient of ", _F, ". (parameter F)");
		if (_F < 0.0 || _F > 1.0) throw "Inbreeding coefficient F must be within [0,1]!";
	}

	// write true allele freq?
	if (parameters().parameterExists("writeTrueAlleleFreq")) {
		const auto outname = parameters().getParameterWithDefault<std::string>("out", "ATLAS_simulations");
		const auto alleleFreqFile = outname + "_trueAlleleFreq.txt.gz";
		logfile().list("Will write true allele frequencies to file '" + alleleFreqFile + "'.");
		_trueFreqFile.open(alleleFreqFile);
		_trueFreqFile.writeHeader({"Chr", "Pos", "Ancestral", "Derived", "derivedFreq", "MAF"});
		_writeTrueAlleleFreq = true;
	}

	// done
	logfile().endIndent();
}

void TSimulatorHW::_fillCumulGenoProb(double f) {
	const double oneMinus_f = 1.0 - f;
	_cumulGenoProb[0]       = _F * oneMinus_f + (1.0 - _F) * oneMinus_f * oneMinus_f;
	_cumulGenoProb[1]       = _cumulGenoProb[0] + (1.0 - _F) * 2.0 * f * oneMinus_f;
	_cumulGenoProb[2]       = 1.0;
}

void TSimulatorHW::_simulateSite(TSimulatorHWSite &site, TSimulatorReference &reference, const std::string &chr, uint64_t pos) {
	// simulate bases
	site.reference   = genometools::Base(randomGenerator().pickOne(_cumulBaseFreq));
	site.alternative = genometools::Base(randomGenerator().pickOne(_mutTable[site.reference]));

	// is the site polymorphic?
	if (randomGenerator().getRand() < _fracPoly) {
		site.isPolymorphic = true;
		site.f             = randomGenerator().getBetaRandom(_alpha, _beta);

		// reference is a random sample from pop with frequency f: flip if ref is alt!
		if (randomGenerator().getRand() < site.f) {
			std::swap(site.reference, site.alternative);
			site.f = 1 - site.f;
		}
	} else {
		site.isPolymorphic = false;
		// is reference diverged?
		site.f             = randomGenerator().getRand() < _referenceDivergence ? 1. : 0.;
	}
	// store reference
	reference[pos] = site.reference;

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

void TSimulatorHW::simulateHaploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					      const genometools::TChromosome &chromosome) {
	// storage

	// now simulate haplotypes
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		TSimulatorHWSite site;
		// simulate site
		_simulateSite(site, reference, chromosome.name, l);

		// polymoprhic or not?
		if (site.isPolymorphic) {
			// simulate genotypes
			for (int i = 0; i < _sampleSize; ++i) {
				if (randomGenerator().getRand() < site.f) {
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

void TSimulatorHW::simulateDiploid(TSimulatorHaplotypes &haplotypes, TSimulatorReference &reference,
					      const genometools::TChromosome &chromosome) {
	// storage

	// now simulate haplotypes
	for (uint64_t l = 0; l < chromosome.length; ++l) {
		TSimulatorHWSite site;
		// simulate site
		_simulateSite(site, reference, chromosome.name, l);

		// polymoprhic or not?
		if (site.isPolymorphic) {
			_fillCumulGenoProb(site.f);

			// simulate genotypes
			for (int i = 0; i < _sampleSize; ++i) {
				int geno = randomGenerator().pickOne(3, _cumulGenoProb);
				if (geno == 0) {
					haplotypes(i, 0, l) = site.reference;
					haplotypes(i, 1, l) = site.reference;
				} else if (geno == 1) {
					if (randomGenerator().getRand() < 0.5) {
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
