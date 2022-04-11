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

#include "TBitSet.h"
#include "TReadGroups.h"
#include "TStrongArray.h"

namespace BAM { class TSequencedBase; }
namespace GenotypeLikelihoods { class TSite; }

namespace GenotypeLikelihoods::RecalEstimatorTools {

std::vector<uint16_t> vectorOfUsed(const std::vector<uint32_t>& counts);
inline uint16_t max(const std::vector<uint32_t> &counts) { return counts.size() - 1; };

class TRecalDataTable {
private:
	uint64_t _counts = 0;
	//all vectors are uint16_t, which is used by seq error models for all covariates
// Object to store for which qualities and positions data is available.
	std::vector<uint32_t> _positions;
	std::vector<uint32_t> _fragmentLengths;
	std::vector<uint32_t> _qualities;
	std::vector<uint32_t> _mappingQualities;

public:
	void clear() noexcept;
	void add(const BAM::TSequencedBase & base);

	constexpr size_t size() const noexcept { return _counts; };
	const std::vector<uint32_t> &positions() const noexcept { return _positions; };
	const std::vector<uint32_t> &fragmentLengths() const noexcept { return _fragmentLengths; };
	const std::vector<uint32_t> &qualities() const noexcept { return _qualities; };
	const std::vector<uint32_t> &mappingQualities() const noexcept { return _mappingQualities; };
};

using TRecalDataTableOneReadGroup = std::array<TRecalDataTable, 2>;

class TRecalDataTables{
private:
	const BAM::TReadGroups* _readGroups = nullptr;
	const BAM::TReadGroupMap* _readGroupMap = nullptr;
	std::vector<TRecalDataTableOneReadGroup> _tables; //tables[readGroup][first/second mate]
	uint64_t _totalCounts = 0;
public:
	TRecalDataTables() = default;
	TRecalDataTables(const BAM::TReadGroups *ReadGroups, const BAM::TReadGroupMap *ReadGroupMapObject)
	    : _readGroups(ReadGroups), _readGroupMap(ReadGroupMapObject), _tables(_readGroupMap->numReadGroupsInUse()){};

	void clear();
	void initialize(const BAM::TReadGroups* ReadGroups, const BAM::TReadGroupMap* ReadGroupMapObject);
	void reset();
	void add(const BAM::TSequencedBase & base);
	void add(const TSite & site);
	void add(const std::vector<TSite> & sites);

	constexpr uint64_t size() const noexcept {return _totalCounts;};
	const TRecalDataTableOneReadGroup& operator[](uint16_t readGroupId) const;
};

//------------------------------------------------
// Classes to keep track of models to estimate
//------------------------------------------------
class TModelStatusEntry{
private:
	coretools::TBitSet<2> _bs{};
public:
	constexpr uint16_t size() const noexcept { return _bs.get<0>() + _bs.get<1>(); };
	constexpr void set(const bool &IsSecondMate) noexcept {
		_bs[IsSecondMate] = true;
	};
	std::string getString() const;
};

enum class ModelStatusTypes : uint8_t {copied=0, noData, littleData, dataButNoRecal};
using TModelStatus = coretools::TStrongArray<TModelStatusEntry, ModelStatusTypes, 4>;

class TModelStati{
private:
	std::unordered_map<uint16_t, TModelStatus> modelStatus;
public:
	void add(uint16_t ReadGroupId);
	TModelStatus& operator[](uint16_t ReadGroupId);
	uint16_t num(ModelStatusTypes Type) const;
	void report(ModelStatusTypes Type, const std::string & Title, const BAM::TReadGroups & ReadGroups) const;
};

}; // namespace GenotypeLikelihoods::RecalEstimatorTools

#endif /* GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_ */
