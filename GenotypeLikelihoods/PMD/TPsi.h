/*
 * PMD/TModels.h
 */

#ifndef GENOTYPELIKELIHOODS_PMD_TPSI_H_
#define GENOTYPELIKELIHOODS_PMD_TPSI_H_

#include <cstdint>
#include <vector>

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"
#include "TReadGroupInfo.h"

#include "TSequencedData.h"

namespace GenotypeLikelihoods::PMD {

enum class Type : size_t {min, CT=min, GA, max};

using TBaseBaseProbabilities =
	coretools::TStrongArray<coretools::TStrongArray<coretools::Probability, genometools::Base>, genometools::Base>;

class TPsi {
	static constexpr coretools::TStrongArray<Type, Type> _flip{{Type::GA, Type::CT}};
	static constexpr coretools::TStrongArray<genometools::Base, Type> _from{{genometools::Base::C, genometools::Base::G}};
	static constexpr coretools::TStrongArray<genometools::Base, Type> _to{{genometools::Base::T, genometools::Base::A}};

	struct SumType {
		struct NumDenom {
			double num;
			double denom;
		};
		struct FromTo {
			int32_t fromTo;
			int32_t fromSum;
			int32_t toFrom;
			int32_t toSum;
		};
		union {
			FromTo fromTo;
			NumDenom numDenom;
		};
	};

	std::vector<coretools::TStrongArray<std::vector<coretools::Probability>, Type>> _tables;
	std::vector<coretools::TStrongArray<std::vector<SumType>, Type>> _tableSums;

	size_t _nSingle = 0;
	size_t _nPaired = 0;
	size_t _S       = 5;
	size_t _Cmax    = 0;

	void _printTable(std::string_view OutputName);
	void _initEnd(size_t Idx, int32_t MinData);
	void _joinTables(size_t From, size_t To) noexcept;
	void _fromString(std::string_view Psi);
	void _parse(const BAM::RGInfo::TInfo & info);
	void _add(Type From_To, const BAM::TSequencedData &data, genometools::Base ref);
	size_t _tIndex(const BAM::TSequencedData &data) const noexcept;

public:

	TPsi();
	TPsi(const BAM::RGInfo::TInfo & info);

	BAM::RGInfo::TInfo info() const;

	coretools::Probability prob(Type From_To, const BAM::TSequencedData &data) const noexcept;

	void add(Type From_To, const BAM::TSequencedData &data, coretools::Probability P_g_I_di,
			 const TBaseBaseProbabilities &P_b_bbar_I_gdij);

	void add(const BAM::TSequencedData &data, genometools::Base ref) noexcept {
		if (data.get<BAM::Flags::Paired>()) ++_nPaired;
		else ++_nSingle;
		_add(Type::CT, data, ref);
		_add(Type::GA, data, ref);
	}

	void estimate() noexcept;
	void estimateInit(std::string_view OutputName, size_t MinData) noexcept;

	void log() const noexcept;
	bool paired() const noexcept {return _nPaired > _nSingle;}
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
