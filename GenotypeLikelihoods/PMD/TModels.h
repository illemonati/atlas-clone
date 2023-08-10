/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef PMD_TMODELS_H_
#define PMD_TMODELS_H_

#include <array>
#include <map>
#include <memory>
#include <stdint.h>
#include <string>
#include <vector>

#include <armadillo>

#include "TReadGroupInfo.h"
#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TFunction.h"
#include "TModel.h"
#include "TReadGroups.h"
#include "coretools/Strings/stringFunctions.h"

namespace BAM {
class TSequencedBase;
}

namespace GenotypeLikelihoods::PMD {
class TModels {
private:
	std::vector<std::unique_ptr<TModel>> _models;
	bool _hasPMD      = false;
	size_t _tableSize = 0;

	void _initializeFromString(const std::string &pmdString);
	std::vector<size_t> _initializeFromFile(const BAM::TReadGroups &ReadGroups, const std::string &filename);
	void _setHasDamage();

public:
	TModels() = default;
	TModels(const std::string &pmdString, const BAM::TReadGroups &ReadGroups,
					  std::vector<size_t> &ReadGroupsWithoutPMD) {
		ReadGroupsWithoutPMD = initialize(pmdString, ReadGroups);
	}
	constexpr bool hasPMD() const noexcept { return _hasPMD; };
	const TModel &operator[](size_t ReadGroupIndex) const noexcept { return *_models[ReadGroupIndex]; }
	TModel &operator[](size_t ReadGroupIndex) noexcept { return *_models[ReadGroupIndex]; }

	std::vector<size_t> initialize(const std::string &pmdString, const BAM::TReadGroups &ReadGroups);
	void initialize(BAM::RGInfo::TReadGroupInfo & RgInfo);
	void writeToFile(const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
	                 std::string_view outputName) const;
	TBaseLikelihoods P_dij(const BAM::TSequencedBase &data,
	                                    const TBaseLikelihoods &baseLikelihoodsNoPMD) const;
	TBaseProbabilities P_bbar(genometools::Base b, const BAM::TSequencedBase &data,
									   const TBaseLikelihoods &baseLikelihoodsNoPMD) const;
	TBaseProbabilities P_bbar(genometools::Genotype g, const BAM::TSequencedBase &data,
									   const TBaseLikelihoods &baseLikelihoodsNoPMD) const;

	void resize(const BAM::TReadGroupMap& ReadGroupMap);
	void add(const BAM::TReadGroupMap& ReadGroupMap, genometools::Base from, BAM::TSequencedBase data) {
		_models[ReadGroupMap.pooledIndex(data.readGroupID)]->add(from, data);
	}

	void estimate(const BAM::TReadGroupMap& ReadGroupMap) {
		for (auto &r : ReadGroupMap.readGroupsInUse()) { _models[r]->estimate(); }
	}

	std::string functionString() const noexcept {
		std::string r;
		for (auto &p: _models) {
			r.append(p->functionString()).append(1, ';');
		}
		return r;
	}
};

}; // namespace GenotypeLikelihoods

#endif /* TPOSTMORTEMDAMAGE_H_ */
