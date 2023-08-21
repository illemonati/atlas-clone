/*
 * TPostMortemDamage.cpp
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */

/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#include "TModels.h"

#include <filesystem>
#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>

#include "TFunction.h"
#include "TModel.h"
#include "genometools/GenotypeTypes.h"
#include "coretools/Main/TError.h"
#include "coretools/Files/TOutputFile.h"
#include "TGenotypeData.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TReadGroupInfo.h"
#include "TSequencedBase.h"
#include "coretools/Types/probability.h"


namespace GenotypeLikelihoods::oldPMD {

using coretools::instances::logfile;
using namespace coretools::str;

void TModels::writeToFile(const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
				    std::string_view outputName) const {
	using coretools::TOutputFile;
	TOutputFile out_PMD(std::string(outputName) + "_PMD.txt", {"readGroup", "type", "pmd"});
	std::array<TOutputFile,2> out_counts; 
	out_counts.front().open(std::string(outputName) + "_PMD_counts.txt", _tableSize + 5);
	out_counts.back().open(std::string(outputName) + "_PMD_countsNormalized.txt", _tableSize + 5);

	for (auto& out_count : out_counts) {
		out_count.write("readGroup", "strand", "fromEnd", "ref", "data");
		for (size_t i = 1; i < _tableSize; ++i) {
			out_count.writeNoDelim("pos_");
			out_count.write(i);
		}
		out_count.writeNoDelim("pos>");
		out_count.writeln(_tableSize - 1);
	}
	// write for each read group
	for (auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r) {
		out_PMD.writeln(r->name_ID, _models[r->id]->typeString(),
					_models[ReadGroupMap.pooledIndex(r->id)]->functionString());
		if (ReadGroupMap.inUse(r->id)) {
			_models[ReadGroupMap.pooledIndex(r->id)]->writeTable(ReadGroups.getName(r->id),out_counts);
		}
	}
}

void TModels::_initializeFromString(const std::string &pmdString) {
	// not a file: initialize all read groups have the same pmd
	logfile().startIndent("PMD function used for all read groups:");
	for (auto &p : _models) { p.reset(makeType(pmdString)); }

	logfile().list(_models[0]->functionString());
	logfile().endIndent();
}

std::vector<size_t> TModels::_initializeFromFile(const BAM::TReadGroups &ReadGroups, const std::string &filename) {
	// create an array of TPMD objects for each read group
	// also works if no parameters are provided (e.g. for estimation)
	// read from file for each read group
	logfile().listFlush("Initializing PMD from file '" + filename + "' ...");
	coretools::TInputFile in(filename, {"readGroup", "pmd"}, "\t", "//");

	// parse file that has structure: ReadGroup, Type, pmd
	std::vector<std::string> vec;
	while (in.read(vec)) {
		if (ReadGroups.readGroupExists(vec[0])) { // ignore if it does not exist
			// get read group
			const size_t readGroupId = ReadGroups.getId(vec[0]);

			// create type
			_models[readGroupId].reset(makeType(vec[1]));
		}
	}
	logfile().done();

	// check if we have a function for all read groups
	// create no-PMD types for all remaining ones and return their indexes
	std::vector<size_t> readGroupsWithoutPMD;
	for (size_t i = 0; i < ReadGroups.size(); ++i) {
		if (!_models[i]) {
			_models[i].reset(makeType("non"));
			readGroupsWithoutPMD.push_back(i);
		}
	}
	return readGroupsWithoutPMD;
}

void TModels::_setHasDamage() {
	// check if there is PMD for at least one read group
	_hasPMD = false;
	for (auto &p : _models) {
		if (p->hasDamage()) {
			_hasPMD = true;
			break;
		}
	}
}

std::vector<size_t> TModels::initialize(const std::string &pmdString, const BAM::TReadGroups &ReadGroups) {
	if (_hasPMD) {
		DEVERROR("Models already initialized!");
	}

	// prepare objects
	std::vector<size_t> readGroupsWithoutPMD;
	_models.resize(ReadGroups.size());

	if (!std::filesystem::exists(pmdString)) {
		_initializeFromString(pmdString);
	} else {
		readGroupsWithoutPMD = _initializeFromFile(ReadGroups, pmdString);
	}

	// check if there is PMD for at least one read group
	_setHasDamage();
	return readGroupsWithoutPMD;
}

void TModels::initialize(BAM::RGInfo::TReadGroupInfo &RgInfo) {
	using BAM::RGInfo::InfoType;
	_models.resize(RgInfo.size());

	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		_models[rg].reset(makeType(RgInfo.getString(rg, InfoType::pmd, TNo::name)));
	}
}

void TModels::resize(const BAM::TReadGroupMap& ReadGroupMap) {
	_tableSize = coretools::instances::parameters().getParameterWithDefault<int>("length", 50) + 1;
	logfile().list("Estimating PMD from the first ", _tableSize - 1, " positions.");
		for (auto &r : ReadGroupMap.readGroupsInUse()) { _models[r]->resize(_tableSize); }
	}

TBaseLikelihoods TModels::P_dij(const BAM::TSequencedBase &data,
                                            const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return _hasPMD
		? _models[data.readGroupID]->baseLikelihoods(data, baseLikelihoodsNoPMD)
		: baseLikelihoodsNoPMD;
}

TBaseProbabilities TModels::P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
													  const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return _hasPMD ? _models[data.readGroupID]->massFunction(b, data, baseLikelihoodsNoPMD)
	               : TNoPMD::massFunctions[b];
}

TBaseProbabilities TModels::P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
													  const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return _hasPMD ? _models[data.readGroupID]->massFunction(g, data, baseLikelihoodsNoPMD)
		: TNoPMD::massFunction(g, baseLikelihoodsNoPMD);
}

}; // namespace GenotypeLikelihoods
