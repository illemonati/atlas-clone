/*
 * RecalEstimatorTools.h
 *
 *  Created on: May 4, 2021
 *      Author: phaentu
 */

#ifndef TRECALDATATABLES_h_
#define TRECALDATATABLES_h_

#include <vector>

#include "TReadGroupMap.h"
#include "coretools/Containers/TStrongArray.h"
#include "SequencingError/TRecalDataTable.h"

namespace GenotypeLikelihoods { class TSite; }

namespace GenotypeLikelihoods::RecalEstimatorTools {


using TRecalDataTableOneReadGroup = coretools::TStrongArray<TRecalDataTable, BAM::Mate>;

class TRecalDataTables{
private:
	const BAM::TReadGroupMap *_readGroupMap = nullptr;
	std::vector<TRecalDataTableOneReadGroup> _tables; //tables[readGroup][first/second mate]
	size_t _size         = 0;
	size_t _N_g1         = 0;

public:
	TRecalDataTables(const BAM::TReadGroupMap &ReadGroupMapObject)
	    : _readGroupMap(&ReadGroupMapObject), _tables(_readGroupMap->size()){}

	void add(const TSite & site);

	constexpr size_t size() const noexcept {return _size;}
	constexpr size_t nSites_g1() const noexcept {return _N_g1;}
	const TRecalDataTableOneReadGroup& operator[](size_t readGroupId) const;

	void write(std::string_view Name) const;
};

}; // namespace GenotypeLikelihoods::RecalEstimatorTools

#endif /* GENOTYPELIKELIHOODS_RECALESTIMATORTOOLS_H_ */
