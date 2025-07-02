/*
 * RecalEstimatorTools.cpp
 *
 *  Created on: May 4, 2021
 *      Author: phaentu
 */

#include "TRecalDataTables.h"

#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Types/TPseudoInt.h"
#include "genometools/Genotypes/Base.h"

#include "TSite.h"
#include "SequencingError/TCovariate.h"

namespace GenotypeLikelihoods::RecalEstimatorTools {
using SequencingError::Covariates;

namespace impl {
void writeTransformed(Covariates C, uint8_t Value, coretools::TOutputFile & OFile) {
	switch(C) {
	case Covariates::Context:
		OFile.write(genometools::base2char(genometools::Base(Value)));
		break;
	case Covariates::FragmentLength:
		OFile.write(coretools::TLogInt::fromLog(Value).linear());
		break;
	case Covariates::MappingQuality:
	case Covariates::Quality:
		OFile.write(Value);
		break;
	case Covariates::Position:
		OFile.write(coretools::TPseudoInt::fromPseudo(Value).linear());
		break;
	default: throw coretools::TDevError("This Covariate does not exist");
	}
}
}

void TRecalDataTables::add(const TSite &site) {
	static size_t Id = 0;
	_size += site.depth();
	if (site.depth() > 1) ++_N_g1;

	for (const auto &b : site) {
		_tables[_readGroupMap->pooledIndex(b.readGroupID)][b.mate()].add(b, Id);
	}
	++Id;
}

const TRecalDataTableOneReadGroup& TRecalDataTables::operator[](size_t readGroupId) const{
	return _tables[ _readGroupMap->pooledIndex(readGroupId) ];
}

void TRecalDataTables::write(std::string_view Name) const {
	using coretools::TStrongArray;
	TStrongArray<coretools::TOutputFile, Covariates> files;

	files[Covariates::Context].open(std::string(Name).append("_contexts.txt.gz"));
	files[Covariates::Context].write("context");

	files[Covariates::FragmentLength].open(std::string(Name).append("_fragmentLengths.txt.gz"));
	files[Covariates::FragmentLength].write("fragmentLength");

	files[Covariates::MappingQuality].open(std::string(Name).append("_mappingQualities.txt.gz"));
	files[Covariates::MappingQuality].write("mappingQuality");

	files[Covariates::Position].open(std::string(Name).append("_positions.txt.gz"));
	files[Covariates::Position].write("position");

	files[Covariates::Quality].open(std::string(Name).append("_qualities.txt.gz"));
	files[Covariates::Quality].write("quality");

	std::vector<const TRecalDataTable*> usedTables;
	for (auto rg: _readGroupMap->readGroupsInUse()) {
		const auto & table = _tables[rg];
		if (table.back().size() > 0) {
			for (auto& f: files) f.writeNoDelim("RG_", rg, "_mate1").writeDelim();
			for (auto& f: files) f.writeNoDelim("RG_", rg, "_mate2").writeDelim();
			usedTables.push_back(&table.front());
			usedTables.push_back(&table.back());
		}
		else if (table.front().size() > 0) {
			for (auto& f: files) f.writeNoDelim("RG_", rg).writeDelim();
			usedTables.push_back(&table.front());
		}
	}
	for (auto& f: files) f.endln();

	for (auto c = Covariates::min; c < Covariates::max; ++c) {
		size_t N = 0;
		for (auto pt : usedTables) {
			const auto & table = *pt;
			N = std::max(N, table[c].size());
		}
		N = std::min<size_t>(N, std::numeric_limits<uint8_t>::max());

		for (uint8_t i = 0; i < N; ++i) {
			impl::writeTransformed(c, i, files[c]);
			for (auto pt : usedTables) {
				const auto & table = *pt;
				if (i < table[c].size()) {
					files[c].write(table[c][i]);
				}
				else files[c].write(0);
			}
			files[c].endln();
		}
	}
}

} // namespace GenotypeLikelihoods::RecalEstimatorTools
