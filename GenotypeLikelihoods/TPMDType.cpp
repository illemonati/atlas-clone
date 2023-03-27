#include "TPMDType.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"
#include <utility>

namespace GenotypeLikelihoods {

using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using genometools::Base;
using namespace coretools::str;

namespace impl {
template<size_t End>
constexpr ReadEnd makeForward() {
	if constexpr (End == 3) return ReadEnd::forward3;
	if constexpr (End == 5) return ReadEnd::forward5;
	else static_assert(End == 3);
}

template<size_t End>
constexpr ReadEnd makeReverse() {
	if constexpr (End == 3) return ReadEnd::reverse3;
	if constexpr (End == 5) return ReadEnd::reverse5;
	else static_assert(End == 3);
}

template<size_t End, Base From, Base To>
auto makeFromTo(const TPMDType::PMDTable &table) {
	// Assumption: From->To pattern is the same for forward and reverse reads from their respective Ends
	const auto N = table.size();
	std::vector<double> from_to;
	from_to.reserve(N);
	std::vector<double> to_from;
	to_from.reserve(N);

	constexpr auto forward = makeForward<End>();
	constexpr auto reverse = makeReverse<End>();

	for (size_t i = 0; i < N; ++i) {
		// CT
		from_to.push_back(table[i][forward][From][To]);
		from_to.back() += table[i][reverse][From][To];
		if (from_to.back() < 100) {
			if (i < 10) {
				UERROR("Not sufficient ", From, "-", To, " data to estimate PMD model at position ", i, ": ", from_to.back(),
					   ", the first 10 positions must have > 100 data points!\nConsider pooling read groups (parameter "
					   "poolReadGroups).");
			}
			from_to.pop_back();
			break;
		}
		double s_from = 0;
		double s_to   = 0;
		for (auto b = Base::min; b < Base::max; ++b) {
			s_from += table[i][forward][From][b];
			s_from += table[i][reverse][From][b];
			s_to += table[i][forward][To][b];
			s_to += table[i][reverse][To][b];
		}
		from_to.back() /= s_from;

		to_from.push_back(table[i][forward][To][From]);
		to_from.back() += table[i][reverse][To][From];
		to_from.back() /= s_to;
	}
	return std::make_pair(from_to, to_from);
}

} // namespace impl

TBaseProbabilities TPMDTypeNone::massFunction(genometools::Genotype g, const TBaseLikelihoods &baseLikelihoodsNoPMD) {
		using namespace genometools;
		const Base a = first(g);
		const Base b = second(g);
		TBaseLikelihoods mf{0};
		mf[a] = baseLikelihoodsNoPMD[a];
		mf[b] = baseLikelihoodsNoPMD[b];
		return TBaseProbabilities::normalize(mf);
	}

enum class Strand : size_t {min, Single=min, Double, max};

template<Strand strand> class TPMDTypeStrand final : public TPMDType {
private:
	std::unique_ptr<TPMDFunction> _pmd5;
	std::unique_ptr<TPMDFunction> _pmd3;

	coretools::Probability _probCT(const BAM::TSequencedBase &data) const noexcept {
		using coretools::Probability;
		if constexpr (strand == Strand::Single) {
			if (data.isReverseStrand()) return coretools::Probability{};
			return data.distFrom3Prime < data.distFrom5Prime ? _pmd3->prob(data.distFrom3Prime)
															 : _pmd5->prob(data.distFrom5Prime);

		} else {
			if (data.distFrom3Prime < data.distFrom5Prime) { // from 3
				return !data.isReverseStrand() ? Probability{} : _pmd3->prob(data.distFrom3Prime);
			} else { // from 5
				return !data.isReverseStrand() ? _pmd5->prob(data.distFrom5Prime) : Probability{};
			}
		}
	}

	coretools::Probability _probGA(const BAM::TSequencedBase &data) const noexcept {
		using coretools::Probability;
		if constexpr (strand == Strand::Single) {
			if (!data.isReverseStrand()) return coretools::Probability{};
			return data.distFrom3Prime < data.distFrom5Prime ? _pmd3->prob(data.distFrom3Prime)
															 : _pmd5->prob(data.distFrom5Prime);

		} else {
			if (data.distFrom3Prime < data.distFrom5Prime) { // from 3
				return !data.isReverseStrand() ? _pmd3->prob(data.distFrom3Prime) : Probability{};
			} else { // from 5
				return !data.isReverseStrand() ? Probability{} : _pmd5->prob(data.distFrom5Prime);
			}
		}
	}

	static constexpr coretools::TStrongArray<std::string_view, Strand> _names{{"singleStrand", "doubleStrand"}};
public:
	static constexpr std::string_view name = _names[strand];
	TPMDTypeStrand(const std::vector<std::string> &Details) {
		constexpr size_t nDetails = 3;
		if (Details.size() != nDetails) {
			UERROR("Cannot initialize PMD type ", name, ": expect ", nDetails, " entries but found ", Details.size(),
				   "!", "\nProvided string: '", concatenateString(Details, ':'), "'.", "\nExpect string of the form '",
				   name, "':functionCT:functionGA'.");
		}
		_pmd5.reset(makeFunction(Details[1]));
		_pmd3.reset(makeFunction(Details[2]));
	}

	bool hasDamage() const noexcept override { return _pmd5->hasDamage() || _pmd3->hasDamage(); };
	std::string functionString() const noexcept override {
		return std::string(name).append(1, ':').append(_pmd5->string()).append(1, ':').append(_pmd3->string());
	}
	std::string_view typeString() const noexcept override { return name; }

	void estimate(const PMDTable &table) override {
		logfile().startIndent("Learning 5' C-T pattern:");
		const auto [C_T5, T_C5] = impl::makeFromTo<5, Base::C, Base::T>(table);
		_pmd5->learn(C_T5, T_C5);
		logfile().endIndent();
		if constexpr (strand == Strand::Single) {
			logfile().startIndent("Learning 3' C-T pattern:");
			const auto [C_T3, T_C3] = impl::makeFromTo<3, Base::C, Base::T>(table);
			_pmd3->learn(C_T3, T_C3);
			logfile().endIndent();
		} else {
			logfile().startIndent("Learning 3' G-A pattern:");
			const auto [G_A3, A_G3] = impl::makeFromTo<3, Base::G, Base::A>(table);
			_pmd3->learn(G_A3, A_G3);
			logfile().endIndent();
		}
	}

	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &data,
									 const TBaseLikelihoods &baseLikelihoodsNoPMD) const override {
		const auto pCT = _probCT(data);
		const auto pGA = _probGA(data);
		TBaseLikelihoods baseLikelihoods(baseLikelihoodsNoPMD);

		baseLikelihoods[Base::C] = (1.0 - pCT) * baseLikelihoodsNoPMD[Base::C] + pCT * baseLikelihoodsNoPMD[Base::T];
		baseLikelihoods[Base::G] = (1.0 - pGA) * baseLikelihoodsNoPMD[Base::G] + pGA * baseLikelihoodsNoPMD[Base::A];
		return baseLikelihoods;
	}

	TBaseProbabilities massFunction(genometools::Base b, const BAM::TSequencedBase &data,
									const TBaseLikelihoods &baseLikelihoodsNoPMD) const override {
		const auto pCT = _probCT(data);
		const auto pGA = _probGA(data);

		switch (b) {
		case Base::A: return TBaseProbabilities::normalize({1., 0., 0., 0.});
		case Base::C: {
			return TBaseProbabilities::normalize(
				{0., (1. - pCT) * baseLikelihoodsNoPMD[Base::C], 0., pCT * baseLikelihoodsNoPMD[Base::T]});
		}
		case Base::G: {
			return TBaseProbabilities::normalize(
				{pGA * baseLikelihoodsNoPMD[Base::A], 0., (1. - pGA) * baseLikelihoodsNoPMD[Base::G], 0.});
		}
		default: return TBaseProbabilities::normalize({0., 0., 0., 1.}); // case Base::T
		}
	}

	TBaseProbabilities massFunction(genometools::Genotype g, const BAM::TSequencedBase &data,
									const TBaseLikelihoods &baseLikelihoodsNoPMD) const override {
		const auto pCT = _probCT(data);
		const auto pGA = _probGA(data);

		using namespace genometools;
		std::array<Base, 2> bases{first(g), second(g)};
		TBaseData mf{0};
		for (const auto a : bases) {
			switch (a) {
			case Base::A: mf[Base::A] += baseLikelihoodsNoPMD[Base::A]; break;
			case Base::C: {
				mf[Base::C] += (1. - pCT) * baseLikelihoodsNoPMD[Base::C];
				mf[Base::T] += pCT * baseLikelihoodsNoPMD[Base::T];
			} break;
			case Base::G: {
				mf[Base::A] += pGA * baseLikelihoodsNoPMD[Base::A];
				mf[Base::G] += (1. - pGA) * baseLikelihoodsNoPMD[Base::G];
			} break;
			default: mf[Base::T] += baseLikelihoodsNoPMD[Base::T];
			}
		}
		return TBaseProbabilities::normalize(mf);
	}

	virtual void simulate(BAM::TSequencedBase &data) const override {
		const auto pCT = _probCT(data);
		const auto pGA = _probGA(data);
		auto &base     = data.base;

		if (base == Base::C) {
			if (randomGenerator().getRand() < pCT) base = Base::T;
		} else if (base == Base::G) {
			if (randomGenerator().getRand() < pGA) base = Base::A;
		}
	}
};

TPMDType *makeType(std::string_view pmdString) {
	// split by ':'
	std::vector<std::string> details;
	fillContainerFromString(pmdString, details, ":");

	// switch type
	if (details[0] == TPMDTypeNone::name) return new TPMDTypeNone;
	if (details[0] == TPMDTypeStrand<Strand::Single>::name) return new TPMDTypeStrand<Strand::Single>(details);
	if (details[0] == TPMDTypeStrand<Strand::Double>::name) return new TPMDTypeStrand<Strand::Double>(details);

	UERROR("Cannot initialize PMD: unknown PMD type '", details[0], "'!\nUse ", TPMDTypeNone::name, " or ",
		   TPMDTypeStrand<Strand::Single>::name, " or ", TPMDTypeStrand<Strand::Double>::name, ".");
}
}
