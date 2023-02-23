/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TPOSTMORTEMDAMAGE_H_
#define TPOSTMORTEMDAMAGE_H_

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
#include "TPMDTables.h"
#include "TPMDFunction.h"
#include "TPMDType.h"
#include "TReadGroups.h"
#include "coretools/Strings/stringFunctions.h"

namespace BAM {
class TSequencedBase;
}

namespace GenotypeLikelihoods {

class TPostMortemDamage {
private:
	std::vector<std::unique_ptr<TPMDType>> _pmdObjects;
	bool _hasPMD = false;

	void _initializeFromString(const std::string &pmdString);
	std::vector<uint16_t> _initializeFromFile(const BAM::TReadGroups &ReadGroups, const std::string &filename);
	void _setHasDamage();

public:
	TPostMortemDamage() = default;
	TPostMortemDamage(const std::string &pmdString, const BAM::TReadGroups &ReadGroups,
					  std::vector<uint16_t> &ReadGroupsWithoutPMD) {
		ReadGroupsWithoutPMD = initialize(pmdString, ReadGroups);
	}
	constexpr bool hasPMD() const noexcept { return _hasPMD; };
	const TPMDType &operator[](uint16_t ReadGroupIndex) const noexcept { return *_pmdObjects[ReadGroupIndex]; }
	TPMDType &operator[](uint16_t ReadGroupIndex) noexcept { return *_pmdObjects[ReadGroupIndex]; }

	std::vector<uint16_t> initialize(const std::string &pmdString, const BAM::TReadGroups &ReadGroups);
	void initialize(BAM::RGInfo::TReadGroupInfo & RgInfo);
	void writeToFile(const BAM::TReadGroups &ReadGroups, const std::string filename) const;
	void writeToFile(const BAM::TReadGroups &ReadGroups, const BAM::TReadGroupMap &ReadGroupMap,
	                 const std::string filename) const;
	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &data,
	                                    const TBaseLikelihoods &baseLikelihoodsNoPMD) const;
	TBaseProbabilities massFunction(genometools::Base b, const BAM::TSequencedBase &data,
									   const TBaseLikelihoods &baseLikelihoodsNoPMD) const;

	std::string functionString() const noexcept {
		std::string r;
		for (auto &p: _pmdObjects) {
			r.append(p->functionString()).append(1, ';');
		}
		return r;
	}

};

}; // namespace GenotypeLikelihoods

#endif /* TPOSTMORTEMDAMAGE_H_ */
