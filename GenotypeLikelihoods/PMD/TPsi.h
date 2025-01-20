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

	coretools::TStrongArray<coretools::TStrongArray<std::vector<coretools::Probability>, Type>, BAM::End> _tables;
	coretools::TStrongArray<coretools::TStrongArray<std::vector<SumType>, Type>, BAM::End> _tableSums;

	size_t _nSingle = 0;
	size_t _nPaired = 0;

	void _printTable(std::string_view FName);
	void _initEnd(BAM::End e, int32_t MinData);
	void _joinTables() noexcept;
	void _fromString(std::string_view Psi);
	void _parse(const BAM::RGInfo::TInfo & info);

	template<Type From_To>
	void _add(const BAM::TSequencedData &data, coretools::Probability &P_g_I_di,
			  const TBaseBaseProbabilities &P_b_bbar_I_gdij) {
		using genometools::Base;
		using BAM::End;

		constexpr auto From = _from[From_To];
		constexpr auto To   = _to[From_To];

		const auto end      = paired() ? End::from5 : data.end();
		const auto realType = data.get<BAM::Flags::ReversedStrand>() ? _flip[From_To] : From_To;
		auto &tSum          = _tableSums[end][realType];
		if (tSum.empty()) return; // wrong pattern

		const auto pos  = std::min<size_t>(tSum.size() - 1, data.dist(end).pseudo());

		tSum[pos].numDenom.num   += P_b_bbar_I_gdij[From][To] * P_g_I_di;
		tSum[pos].numDenom.denom += (P_b_bbar_I_gdij[From][To] + P_b_bbar_I_gdij[From][From]) * P_g_I_di;
	}

	template<Type From_To>
	void _add(const BAM::TSequencedData &data, genometools::Base ref) {
		using genometools::Base;
		using BAM::End;

		constexpr auto From = _from[From_To];
		constexpr auto To   = _to[From_To];

		const auto end  = data.end();
		const auto pos  = data.dist(end).pseudo();
		const auto base = data.base;

		const auto realType = data.get<BAM::Flags::ReversedStrand>() ? _flip[From_To] : From_To;
		auto &tSum          = _tableSums[end][realType];
		if (tSum.size() <= pos) tSum.resize(pos + 1, {{{0,0,0,0}}});

		if (ref == From) {
			if (base == To) tSum[pos].fromTo.fromTo += 1.;
			tSum[pos].fromTo.fromSum += 1.;
		} else if (ref == To) {
			if (base == From) tSum[pos].fromTo.toFrom += 1.;
			tSum[pos].fromTo.toSum += 1.;
		}
	}

public:

	TPsi();
	TPsi(const BAM::RGInfo::TInfo & info);

	BAM::RGInfo::TInfo info() const;

	void reset(const BAM::RGInfo::TInfo & info);

	template<Type From_To>
	coretools::Probability prob(const BAM::TSequencedData &data) const noexcept {
		using BAM::End;
		const auto realType = data.get<BAM::Flags::ReversedStrand>() ? _flip[From_To] : From_To;
		const auto end      = data.get<BAM::Flags::Paired>() ? End::from5 : data.end();
		const auto pos      = data.dist(end).pseudo();

		const auto &table = _tables[end][realType];
		if (pos >= table.size()) return table.back();
		return table[pos];
	}

	void addGA(const BAM::TSequencedData &data, coretools::Probability P_g_I_di,
			 const TBaseBaseProbabilities& bbProbs) noexcept {
		_add<Type::GA>(data, P_g_I_di, bbProbs);
	}

	void addCT(const BAM::TSequencedData &data, coretools::Probability P_g_I_di,
			 const TBaseBaseProbabilities& bbProbs) noexcept {
		_add<Type::CT>(data, P_g_I_di, bbProbs);
	}

	void add(const BAM::TSequencedData &data, genometools::Base ref) noexcept {
		if (data.get<BAM::Flags::Paired>()) ++_nPaired;
		else ++_nSingle;
		_add<Type::CT>(data, ref);
		_add<Type::GA>(data, ref);
	}

	void estimate() noexcept;
	void estimateInit(std::string_view OutputName, size_t MinData) noexcept;

	void log() const noexcept;

	bool paired() const noexcept {return _nPaired > _nSingle;}

	const std::vector<coretools::Probability>& vals(BAM::End E, Type T) const noexcept {
		return _tables[E][T];
	}
};
} // namespace GenotypeLikelihoods::PMD

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORMODELS_H_ */
