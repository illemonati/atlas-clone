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

#include "TPostMortemDamage.h"

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

#include "TPMDFunction.h"
#include "TPMDType.h"
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


namespace GenotypeLikelihoods {

using coretools::instances::logfile;
using namespace coretools::str;


void TPostMortemDamage::writeToFile(const BAM::TReadGroups &ReadGroups, const std::string filename) const {
	std::vector<std::string> header = {"readGroup", "type", "pmd"};
	coretools::TOutputFile out(filename, header);

	// write for each read group
	for (auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r)
		out.writeln(r->name_ID, _pmdObjects[r->id()]->typeString(), _pmdObjects[r->id()]->functionString());
}

void TPostMortemDamage::writeToFile(const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
				    const std::string filename) const {
	std::vector<std::string> header = {"readGroup", "type", "pmd"};
	coretools::TOutputFile out(filename, header);

	// write for each read group
	for (auto r = ReadGroups.cbegin(); r != ReadGroups.cend(); ++r)
		out.writeln(r->name_ID, _pmdObjects[r->id()]->typeString(),
					_pmdObjects[ReadGroupMap.pooledIndex(r->id())]->functionString());
}

void TPostMortemDamage::_initializeFromString(const std::string &pmdString) {
	// not a file: initialize all read groups have the same pmd
	logfile().startIndent("PMD function used for all read groups:");

	for (auto &p : _pmdObjects) { p.reset(makeType(pmdString)); }

	// report
	logfile().list(_pmdObjects[0]->functionString());
	logfile().endIndent();
}

std::vector<size_t> TPostMortemDamage::_initializeFromFile(const BAM::TReadGroups &ReadGroups, const std::string &filename) {
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
			uint16_t readGroupId = ReadGroups.getId(vec[0]);

			// create type
			_pmdObjects[readGroupId].reset(makeType(vec[1]));
		}
	}
	logfile().done();

	// check if we have a function for all read groups
	// create no-PMD types for all remaining ones and return their indexes
	std::vector<size_t> readGroupsWithoutPMD;
	for (size_t i = 0; i < ReadGroups.size(); ++i) {
		if (!_pmdObjects[i]) {
			_pmdObjects[i].reset(makeType("non"));
			readGroupsWithoutPMD.push_back(i);
		}
	}
	return readGroupsWithoutPMD;
}


void TPostMortemDamage::_setHasDamage() {
	// check if there is PMD for at least one read group
	_hasPMD = false;
	for (auto &p : _pmdObjects) {
		if (p->hasDamage()) {
			_hasPMD = true;
			break;
		}
	}
}

std::vector<size_t> TPostMortemDamage::initialize(const std::string &pmdString, const BAM::TReadGroups &ReadGroups) {
	if (_hasPMD) {
		DEVERROR("Models already initialized!");
	}

	// prepare objects
	std::vector<size_t> readGroupsWithoutPMD;
	_pmdObjects.resize(ReadGroups.size());

	if (!std::filesystem::exists(pmdString)) {
		_initializeFromString(pmdString);
	} else {
		readGroupsWithoutPMD = _initializeFromFile(ReadGroups, pmdString);
	}

	// check if there is PMD for at least one read group
	_setHasDamage();
	return readGroupsWithoutPMD;
}

void TPostMortemDamage::initialize(BAM::RGInfo::TReadGroupInfo &RgInfo) {
	using BAM::RGInfo::InfoType;

	_pmdObjects.resize(RgInfo.size());

	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		_pmdObjects[rg].reset(makeType(RgInfo.getString(rg, InfoType::pmd, TPMDFunctionNoPMD::name)));
	}
}

TBaseLikelihoods TPostMortemDamage::baseLikelihoods(const BAM::TSequencedBase &data,
                                            const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return _hasPMD
		? _pmdObjects[data.readGroupID]->baseLikelihoods(data, baseLikelihoodsNoPMD)
		: baseLikelihoodsNoPMD;
}

TBaseProbabilities TPostMortemDamage::massFunction(genometools::Base b, const BAM::TSequencedBase &data,
													  const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return _hasPMD ? _pmdObjects[data.readGroupID]->massFunction(b, data, baseLikelihoodsNoPMD)
	               : TPMDTypeNone::massFunctions[b];
}

TBaseProbabilities TPostMortemDamage::massFunction(genometools::Genotype g, const BAM::TSequencedBase &data,
													  const TBaseLikelihoods &baseLikelihoodsNoPMD) const {
	return _hasPMD ? _pmdObjects[data.readGroupID]->massFunction(g, data, baseLikelihoodsNoPMD)
		: TPMDTypeNone::massFunction(g, baseLikelihoodsNoPMD);
}

}; // namespace GenotypeLikelihoods
