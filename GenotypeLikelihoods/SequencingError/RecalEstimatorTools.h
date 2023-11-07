/*
 * RecalEstimatorTools.h
 *
 *  Created on: May 4, 2021
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_
#define GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_

#include <stddef.h>
#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "SequencingError/TCovariate.h"
#include "TSequencedBase.h"
#include "coretools/Containers/TBitSet.h"
#include "TReadGroups.h"
#include "coretools/Containers/TStrongArray.h"

namespace BAM { class TSequencedBase; }
namespace GenotypeLikelihoods { class TSite; }

namespace GenotypeLikelihoods::RecalEstimatorTools {

class TRecalDataTable {
private:
	uint64_t _counts = 0;
	//all vectors are uint16_t, which is used by seq error models for all covariates
// Object to store for which qualities and positions data is available.
	coretools::TStrongArray<std::vector<size_t>, SequencingError::Covariates> _tables;

public:
	void add(const BAM::TSequencedBase & base);

	constexpr size_t size() const noexcept { return _counts; }
	const std::vector<size_t>& operator[](SequencingError::Covariates cov) const {return _tables[cov];}
};

using TRecalDataTableOneReadGroup = coretools::TStrongArray<TRecalDataTable, BAM::Mate>;

class TRecalDataTables{
private:
	const BAM::TReadGroupMap* _readGroupMap = nullptr;
	std::vector<TRecalDataTableOneReadGroup> _tables; //tables[readGroup][first/second mate]
	size_t _size = 0;
	size_t _N_g1 = 0;

public:
	TRecalDataTables(const BAM::TReadGroupMap &ReadGroupMapObject)
	    : _readGroupMap(&ReadGroupMapObject), _tables(_readGroupMap->numReadGroupsInUse()){}

	void add(const TSite & site);

	constexpr size_t size() const noexcept {return _size;}
	constexpr size_t nSites_g1() const noexcept {return _N_g1;}
	const TRecalDataTableOneReadGroup& operator[](size_t readGroupId) const;

	void write(std::string_view Name) const;
};

}; // namespace GenotypeLikelihoods::RecalEstimatorTools

#endif /* GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_ */
