/*
 * TPMDEstimator.cpp
 *
 *  Created on: Jun 3, 2020
 *      Author: phaentu
 */

#include "TPMDEstimator.h"

#include <exception>
#include <memory>
#include <vector>

#include "coretools/Files/TOutputFile.h"
#include "genometools/GenotypeTypes.h"
#include "TAlignment.h"
#include "TBamFile.h"
#include "TGenome.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TReadGroups.h"

namespace GenomeTasks{
using namespace std::literals;
using coretools::instances::logfile;
using coretools::instances::parameters;

//----------------------------------------
// TPMDEstimator.h
//----------------------------------------
TPMDEstimator::TPMDEstimator(): TGenome_parsed(), _readGroupMap(_bamFile.readGroups(), parameters().getParameter<std::string>("poolReadGroups", false)) {
	//make sure there is pmd
	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.postMortemDamageModels();
	if (_genotypeLikelihoodCalculator.hasPMD() && !parameters().parameterExists("reestimate")) {
		logfile().list("PMD model already exists, will reestimate it.");
	}
	if (!_genotypeLikelihoodCalculator.hasPMD()) {
		pmd.initialize(parameters().getParameterWithDefault("pmd", "doubleStrand:Empiric:Empiric"s), _bamFile.readGroups());
	}

	//make sure it has a reference
	_openReference(true);

	//parse estimation parameters
	logfile().startIndent("Parameters for PMD Estimation:");
	_tableSize = parameters().getParameterWithDefault<int>("length", 50) + 1;
	logfile().list("Estimating PMD from the first ", _tableSize - 1, " positions.");

	//create PMD tables
	//_pmdTables.initialize(&_bamFile.readGroups(), _maxLengthForInference, &_readGroupMap);
	const auto &readGroups = _bamFile.readGroups();
	_pmdTables.resize(readGroups.size());
	for (auto rg = readGroups.cbegin(); rg != readGroups.cend(); ++rg) {
		const auto rg_i = _readGroupMap.pooledIndex(rg->id());
		_pmdTables[rg_i].resize(_tableSize);
	}
};

void TPMDEstimator::_handleAlignment() {
	using GenotypeLikelihoods::ReadEnd;
	using genometools::flipped;
	for (size_t d = 0; d < _alignment.size(); ++d) {
		if (_alignment.isAlignedAtInternalPos(d)) {
			const auto data     = _alignment[d];
			const auto from     = _alignment.referenceAtInternalPos(d);
			const auto to       = data.base;
			const auto rg_i     = _readGroupMap.pooledIndex(data.readGroupID);
			const auto from3    = data.distFrom3Prime < data.distFrom5Prime;
			const auto pos      = std::min<size_t>(_tableSize - 1, from3 ? data.distFrom3Prime : data.distFrom5Prime);
			if (data.isReverseStrand()) {
				const auto readEnd = from3 ? ReadEnd::reverse3 : ReadEnd::reverse5;
				_pmdTables[rg_i][pos][readEnd][flipped(from)][flipped(to)]++;
			} else {
				const auto readEnd = from3 ? ReadEnd::forward3 : ReadEnd::forward5;
				_pmdTables[rg_i][pos][readEnd][from][to]++;
			}
		}
	}
};

void TPMDEstimator::_write(bool normalized) {
	using GenotypeLikelihoods::ReadEnd;
	using namespace genometools;
	const auto fn = _outputName + "_PMD_counts" + (normalized ? "Normalized" : "") + ".txt";
	coretools::TOutputFile out(fn, _tableSize + 5);
	// Header
	out.write("readGroup", "strand", "fromEnd", "ref", "data");
	for (size_t i = 1; i < _tableSize; ++i) {
		out.writeNoDelim("pos_");
		out.write(i);
	}
	out.writeNoDelim("pos>");
	out.writeln(_tableSize - 1);

	constexpr auto directions = []() {
		coretools::TStrongArray<std::array<std::string_view, 2>, ReadEnd> ar{};
		ar[ReadEnd::forward3] = {"forward", "3'"};
		ar[ReadEnd::forward5] = {"forward", "5'"};
		ar[ReadEnd::reverse3] = {"reverse", "3'"};
		ar[ReadEnd::reverse5] = {"reverse", "5'"};
		return ar;
	}();

	const auto &rg = _bamFile.readGroups();
	for (size_t rg_i = 0; rg_i < rg.size(); ++rg_i) {
		if (rg.readGroupInUse(rg_i)) {
			for (auto j = ReadEnd::min; j < ReadEnd::max; ++j) {
				for (Base f = Base::min; f <= Base::max; ++f) {
					std::vector<size_t> sums(_tableSize, 0.);
					for (size_t i = 0; i < _tableSize; ++i) {
						for (Base t = Base::min; t < Base::max; ++t) { sums[i] += _pmdTables[rg_i][i][j][f][t]; }
					}

					for (Base t = Base::min; t <= Base::max; ++t) {
						out.write(rg.getName(rg_i), directions[j], f, t);
						for (size_t i = 0; i < _tableSize; ++i) {
							if (normalized)
								out.write(static_cast<double>(_pmdTables[rg_i][i][j][f][t]) / sums[i]);
							else
								out.write(_pmdTables[rg_i][i][j][f][t]);
						}
						out.endln();
					}
					if (!normalized) {
						out.writeln(rg.getName(rg_i), directions[j], f, "sum", sums);
					}
				}
			}
		}
	}
}

void TPMDEstimator::run(){
	// 1) traverse BAM to assemble PMD Tables
	_traverseBAMPassedQC();

	// 2) print PMD tables
	logfile().startIndent("Writing PMD tables:");

	//counts
	std::string filename = _outputName + "_PMD_counts.txt";
	logfile().listFlush("Writing PMD counts to '" + filename + "' ...");
	_write(false);
	//_pmdTables.write(filename, false);
	logfile().done();

	//normalized counts
	filename = _outputName + "_PMD_countsNormalized.txt";
	logfile().listFlush("Writing PMD normalized counts to '" + filename + "' ...");
	_write(true);
	//_pmdTables.write(filename, true);
	logfile().done();
	logfile().endIndent();

	// 3) estimate models
	using coretools::instances::parameters;
	GenotypeLikelihoods::TPostMortemDamage& pmd = _genotypeLikelihoodCalculator.postMortemDamageModels();

	//estimate all models with data, i.e. only one model per pool
	pmd.estimate(_readGroupMap, _pmdTables);

	//writing PMD file
	filename = _outputName + "_PMD.txt";
	pmd.writeToFile(_bamFile.readGroups(), _readGroupMap, filename);
};

}; // end namespace
