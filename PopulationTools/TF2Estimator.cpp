//
// Created by reynac on 2/22/22.
//

#include "TF2Estimator.h"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <ostream>

#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/TTimer.h"
#include "genometools/VCF/TVcfParser.h"
#include "coretools/Strings/stringFunctions.h"

namespace PopulationTools {
using coretools::TOutputFile;
using namespace coretools::instances;

//------------------------------------------------
// TF2Estimator
//------------------------------------------------

void TF2Estimator::run() {
	_openVCF();
	_matchSamples();

	_readVCF();
}

void TF2Estimator::_openVCF() {
	// open VCF
	std::string vcfFilename = parameters().get("vcf");
	_reader.initialize(false);
	_reader.openVCF(vcfFilename);

	// read output name
	auto tmp = coretools::str::readBeforeLast(vcfFilename, ".vcf");
	_outname = parameters().get("out", tmp);
	logfile().list("Writing output files with prefix '", _outname, "'. (specify with 'out')");
}

void TF2Estimator::_matchSamples() {
	if (parameters().exists("samples")) {
		_samples.readSamples(parameters().get<std::string>("samples"));
	}

	if (_samples.hasSamples()) {
		_samples.fillVCFOrder(_reader.getSampleVCFNames());
	} else {
		_samples.readSamplesFromVCFNames(_reader.getSampleVCFNames());
	}
}

void TF2Estimator::_readVCF() {

	// caclulate total comparisons and diff sites
	std::vector<size_t> countsDiff(_samples.numSamples() * _samples.numSamples());
	// traverse VCF
	logfile().startIndent("Parsing VCF file:");
	while (_reader.readDataFromVCF(_samples)) {
		// compare genotypes one ind vs all others at current line
		for (size_t s1 = 0; s1 < _samples.numSamples() - 1; ++s1) {
			if (!_reader.sampleIsMissing(_samples, s1)) {
				auto genotype_s1 = _reader.biallelicGenotype(_samples, s1);

				for (size_t s2 = s1 + 1; s2 < _samples.numSamples(); ++s2) {
					if (!_reader.sampleIsMissing(_samples, s2)) {
						auto genotype_s2 = _reader.biallelicGenotype(_samples, s2);
						// store total # comparison per combination lower triangle
						++countsDiff[(s2 * _samples.numSamples()) + s1];
						if (genotype_s1 != genotype_s2) {
							// store # diff Sites per combination upper triangle
							++countsDiff[(s1 * _samples.numSamples()) + s2];
						}
					}
				}
			}
		}
	}
	logfile().endIndent();
	_reader.concludeFilters();

	calculateF2(countsDiff);
}

void TF2Estimator::calculateF2(const std::vector<size_t> &countsDiff) {
	// calculate sample F2
	std::vector<double> sampleF2(_samples.numSamples() * _samples.numSamples());
	for (size_t s1 = 0; s1 < _samples.numSamples() - 1; ++s1) {
		for (size_t s2 = s1 + 1; s2 < _samples.numSamples(); ++s2) {
			// diff Sites / total compared sites
			if (countsDiff[(s2 * _samples.numSamples()) + s1] != 0) {
				sampleF2[(s1 * _samples.numSamples()) + s2] =
				    static_cast<float>(countsDiff[(s1 * _samples.numSamples()) + s2]) /
				    countsDiff[(s2 * _samples.numSamples()) + s1];
				sampleF2[(s2 * _samples.numSamples()) + s1] = sampleF2[(s1 * _samples.numSamples()) + s2];
			}
		}
	}

	// calculate within and between population average F2
	std::vector<double> popF2(_samples.numPopulations() * _samples.numPopulations());
	for (size_t p1 = 0; p1 < _samples.numPopulations(); ++p1) {
		for (size_t p2 = p1; p2 < _samples.numPopulations(); ++p2) {

			if (p1 == p2) {
				size_t counts = 0;
				for (size_t s1 = _samples.startIndex(p1);
				     s1 < _samples.startIndex(p1) + _samples.numSamplesInPop(p1) - 1; ++s1) {
					for (size_t s2 = s1 + 1; s2 < _samples.startIndex(p2) + _samples.numSamplesInPop(p2); ++s2) {
						popF2[p1 * _samples.numPopulations() + p2] += sampleF2[s1 * _samples.numSamples() + s2];
						++counts;
					}
				}
				if (counts > 0) { popF2[p1 * _samples.numPopulations() + p2] /= counts; }
			} else {
				size_t counts = 0;
				for (size_t s1 = _samples.startIndex(p1); s1 < _samples.startIndex(p1) + _samples.numSamplesInPop(p1);
				     ++s1) {
					for (size_t s2 = _samples.startIndex(p2);
					     s2 < _samples.startIndex(p2) + _samples.numSamplesInPop(p2); ++s2) {
						popF2[p1 * _samples.numPopulations() + p2] += sampleF2[s1 * _samples.numSamples() + s2];
						++counts;
					}
				}
				popF2[p1 * _samples.numPopulations() + p2] /= counts;
				popF2[p2 * _samples.numPopulations() + p1] = popF2[p1 * _samples.numPopulations() + p2];
			}
		}
	}
	writeF2(countsDiff, sampleF2, popF2);
};

void TF2Estimator::writeF2(const std::vector<size_t> &countsDiff, const std::vector<double> &sampleF2,
                           const std::vector<double> &popF2) {
	// open output file for counts
	std::string filename = _outname + "_counts.txt";
	logfile().listFlush("Writing counts of different/compared sites in the upper/lower triangle to file '" + filename +
	                    "' ...");
	std::vector<std::string> header = {"Sample"};
	for (size_t s = 0; s < _samples.numSamples(); ++s) { header.emplace_back(_samples.sampleName(s)); }
	TOutputFile out(filename, header);

	for (size_t s = 0; s < _samples.numSamples(); ++s) {
		size_t tmp                    = s * _samples.numSamples();
		std::vector<size_t> subvector = {countsDiff.begin() + tmp, countsDiff.begin() + tmp + (_samples.numSamples())};
		out.writeln(_samples.sampleName(s), subvector);
	}
	logfile().done();

	// Sample F2
	filename = _outname + "_sampleF2.txt";
	logfile().listFlush("Writing sample F2 results to file '" + filename + "' ...");
	// open output file for sample F2
	TOutputFile outF2(filename, header);
	for (size_t s = 0; s < _samples.numSamples(); ++s) {
		const size_t tmp              = s * _samples.numSamples();
		std::vector<double> subvector = {sampleF2.begin() + tmp, sampleF2.begin() + tmp + (_samples.numSamples())};
		outF2.writeln(_samples.sampleName(s), subvector);
	}
	logfile().done();

	// Pop F2
	filename = _outname + "_popF2.txt";
	logfile().listFlush("Writing population F2 results to file '" + filename + "' ...");
	// open output file for population F2
	std::vector<std::string> Pops = {"Population"};
	for (size_t p = 0; p < _samples.numPopulations(); ++p) { Pops.emplace_back(_samples.getPopulationName(p)); }
	TOutputFile outPopF2(filename, Pops);
	for (size_t p = 0; p < _samples.numPopulations(); ++p) {
		const size_t tmp              = p * _samples.numPopulations();
		std::vector<double> subvector = {popF2.begin() + tmp, popF2.begin() + tmp + _samples.numPopulations()};
		outPopF2.writeln(_samples.getPopulationName(p), subvector);
	}
	logfile().done();
	logfile().endIndent();
}

} // namespace PopulationTools
