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
	std::vector<bool> _positions;
	std::vector<bool> _fragmentLengths;
	std::vector<bool> _qualities;
	std::vector<bool> _mappingQualities;

public:
	void add(const BAM::TSequencedBase & base);

	constexpr size_t size() const noexcept { return _counts; };
	const std::vector<bool> &positions() const noexcept { return _positions; };
	const std::vector<bool> &fragmentLengths() const noexcept { return _fragmentLengths; };
	const std::vector<bool> &qualities() const noexcept { return _qualities; };
	const std::vector<bool> &mappingQualities() const noexcept { return _mappingQualities; };
};

using TRecalDataTableOneReadGroup = coretools::TStrongArray<TRecalDataTable, BAM::Mate>;

class TRecalDataTables{
private:
	const BAM::TReadGroupMap* _readGroupMap = nullptr;
	std::vector<TRecalDataTableOneReadGroup> _tables; //tables[readGroup][first/second mate]
	size_t _size = 0;
	size_t _N_g1 = 0;

public:
	TRecalDataTables() = default;
	TRecalDataTables(const BAM::TReadGroupMap &ReadGroupMapObject)
	    : _readGroupMap(&ReadGroupMapObject), _tables(_readGroupMap->numReadGroupsInUse()){};
	TRecalDataTables(const BAM::TReadGroupMap &ReadGroupMapObject, const std::vector<std::vector<TSite>> & sites)
	    : _readGroupMap(&ReadGroupMapObject), _tables(_readGroupMap->numReadGroupsInUse()){
		for (const auto& s: sites) add(s);
	};

	void initialize(const BAM::TReadGroupMap* ReadGroupMapObject);
	void add(const std::vector<TSite> & sites);

	constexpr size_t size() const noexcept {return _size;};
	constexpr size_t nSites_g1() const noexcept {return _N_g1;};
	const TRecalDataTableOneReadGroup& operator[](uint16_t readGroupId) const;
};

}; // namespace GenotypeLikelihoods::RecalEstimatorTools

#endif /* GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_ */
