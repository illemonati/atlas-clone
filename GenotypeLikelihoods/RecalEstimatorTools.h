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

class TRecalDataTable{
private:
	uint64_t _counts;
	//all vectors are uint16_t, which is used by seq error models for all covariates
	TRecalDataVector<uint16_t> _positions;
	TRecalDataVector<uint16_t> _fragmentLengths;
	TRecalDataVector<uint16_t> _qualities;
	TRecalDataVector<uint16_t> _mappingQualities;

public:

	TRecalDataTable();
	~TRecalDataTable() = default;

	void clear();
	void add(const BAM::TSequencedBase & base);

	size_t size() const;
	const TRecalDataVector<uint16_t>& positions() const;
	const TRecalDataVector<uint16_t>& fragmentLengths() const;
	const TRecalDataVector<uint16_t>& qualities() const;
	const TRecalDataVector<uint16_t>& mappingQualities() const;
};

class TRecalDataTableOneReadGroup{
private:
	std::array<TRecalDataTable, 2> _tables;

public:
	const TRecalDataTable& operator[](const bool & IsSecondMate) const;
	void add(const BAM::TSequencedBase & base);
	void clear();
};

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
	bool _first;
	bool _second;

public:
	TModelStatusEntry();

	uint16_t size();
	void set(const bool & IsSecondMate);
	std::string getString() const;
};

enum class ModelStatusTypes : uint8_t {copied=0, noData, littleData, dataButNoRecal};

class TModelStatus{
private:
	std::array<TModelStatusEntry, 4> _status;

public:
	TModelStatus() = default;
	TModelStatusEntry& operator[](const ModelStatusTypes & Type);
};

class TModelStati{
private:
	std::map<uint16_t, TModelStatus> modelStatus;

public:

	void add(uint16_t ReadGroupId);
	TModelStatus& operator[](uint16_t ReadGroupId);
	uint16_t num(const ModelStatusTypes & Type);
	void report(const ModelStatusTypes & Type, const std::string & Title, const BAM::TReadGroups & ReadGroups, coretools::TLog* Logfile);
};

}; //end namespace RecalEstimatorTools

}; //end namespace GenotypeLikelihoods


#endif /* GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_ */
