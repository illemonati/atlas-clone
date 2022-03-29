/*
 * RecalEstimatorTools.h
 *
 *  Created on: May 4, 2021
 *      Author: phaentu
 */

#ifndef GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_
#define GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_

#include "TSite.h"
#include "TReadGroups.h"
#include "TStrongArray.h"

namespace GenotypeLikelihoods{

namespace RecalEstimatorTools {

//--------------------------------------------------------------------
// TRecalibrationEMDataTable
// Object to store for which qualities and positions data is available.
//--------------------------------------------------------------------
template <typename T>
class TRecalDataVector{
private:
	std::vector<uint32_t> _counts;
public:
	void add(const T & value){
		if(_counts.size() <= value){
			_counts.resize(value+1, 0);
		}
		++_counts[value];
	};

	const T& operator[](const T & value) const{
		return _counts[value];
	};

	void clear(){
		_counts.clear();
	};

	std::vector<T> vectorOfUsed() const{
		//insert all with counts
		std::vector<T> vec;
		for(size_t i = 0; i < _counts.size(); ++i){
			if(_counts[i] > 0){
				vec.push_back(i);
			}
		}
		return vec;
	};

	uint16_t max() const {
		return _counts.size() - 1;
	};
};

std::vector<uint16_t> vectorOfUsed(const std::vector<uint32_t>& counts);

constexpr uint16_t max(const std::vector<uint32_t> &counts) { return counts.size() - 1; };

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

	constexpr size_t size() const noexcept {return _counts;};
	const std::vector<uint32_t>& positions() const;
	const std::vector<uint32_t>& fragmentLengths() const;
	const std::vector<uint32_t>& qualities() const;
	const std::vector<uint32_t>& mappingQualities() const;
};

using TRecalDataTableOneReadGroup = std::array<TRecalDataTable, 2>;

class TRecalDataTables{
private:
	const BAM::TReadGroupMap* _readGroupMap;
	const BAM::TReadGroups* _readGroups;
	std::vector<TRecalDataTableOneReadGroup> _tables; //tables[readGroup][first/second mate]
	uint64_t _totalCounts;

public:
	TRecalDataTables();
	TRecalDataTables(const BAM::TReadGroups* ReadGroups, const BAM::TReadGroupMap* ReadGroupMapObject);
	~TRecalDataTables() = default;

	void clear();
	void initialize(const BAM::TReadGroups* ReadGroups, const BAM::TReadGroupMap* ReadGroupMapObject);
	void reset();
	void add(const BAM::TSequencedBase & base);
	void add(const TSite & site);
	void add(const std::vector<TSite> & sites);

	uint64_t size() const;
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
using TModelStatus = coretools::StrongTypes::TStrongArray<TModelStatusEntry, ModelStatusTypes, 4>;


class TModelStati{
private:
	std::unordered_map<uint16_t, TModelStatus> modelStatus;
public:
	void add(uint16_t ReadGroupId);
	TModelStatus& operator[](uint16_t ReadGroupId);
	uint16_t num(ModelStatusTypes Type);
	void report(ModelStatusTypes Type, const std::string & Title, const BAM::TReadGroups & ReadGroups);
};

}; //end namespace RecalEstimatorTools

}; //end namespace GenotypeLikelihoods


#endif /* GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_ */
