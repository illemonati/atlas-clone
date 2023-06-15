/*
 * TReadGroupPMD.h
 *
 *  Created on: Mars 31, 2023
 *      Author: andreas
 */

#ifndef PMD_TPERREADGROUP_H_
#define PMD_TPERREADGROUP_H_

#include "PMD/TFunction.h"
#include "TModel.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/probability.h"
#include "coretools/devtools.h"
#include "genometools/GenotypeTypes.h"
#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>
namespace GenotypeLikelihoods::PMD {
namespace impl {

static constexpr size_t _N = coretools::index(genometools::Base::max) + 1;
using PMDTable             = std::vector<coretools::TStrongArray<
	coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base, _N>, genometools::Base, _N>, ReadEnd>>;
using TMu =
	coretools::TStrongArray<coretools::TStrongArray<coretools::Probability, genometools::Base>, genometools::Base>;

struct PMDTables {
	std::vector<PMDTable> tables;
	size_t minLength;
	size_t index(size_t length) const noexcept {
		return std::min(std::max<size_t>(minLength, length) - minLength, tables.size() - 1);
	}
};

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

constexpr ReadEnd readEnd(bool is3, bool isReversed) {
	if (is3) {
		if (isReversed) {
			return ReadEnd::reverse3;
		} else {
			return ReadEnd::forward3;
		}
	} else {
		if (isReversed) {
			return ReadEnd::reverse5;
		} else {
			return ReadEnd::forward5;
		}
	}
}

template<size_t End, genometools::Base From, genometools::Base To>
auto makeFromTo(const PMDTable &table) {
	using genometools::Base;
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
template<size_t End>
TMu makeMu(const PMDTable &table) {
	using genometools::Base;
	constexpr auto forward = makeForward<End>();
	constexpr auto reverse = makeReverse<End>();

	coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base>, genometools::Base> from_to{};
	for (size_t p = 0; p < table.size(); ++p) {
		for (auto from = Base::min; from < Base::max; ++from) {
			for (auto to = Base::min; to < Base::max; ++to) {
				from_to[from][to] += table[p][forward][from][to];
				from_to[from][to] += table[p][reverse][from][to];
			}
		}
	}

	TMu mu;
	for (auto from = Base::min; from < Base::max; ++from) {
		OUT(from);
		size_t s = 0;
		for (auto to = Base::min; to < Base::max; ++to) {
			s += from_to[from][to];
		}
		for (auto to = Base::min; to < Base::max; ++to) {
			OUT(to);
			mu[from][to] = static_cast<double>(from_to[from][to])/s;
			OUT(mu[from][to]);
		}
	}
	return mu;
}

template<size_t End, genometools::Base From>
double ll_noPMD(const PMDTable &table, const TMu& mu) {
	using genometools::Base;
	constexpr auto forward = makeForward<End>();
	constexpr auto reverse = makeReverse<End>();

	double ll = 0.;
	for (size_t p = 0; p < table.size(); ++p) {
		double mu_FF = 1.;
		for (auto to = Base::min; to < Base::max; ++to) {
			// take "reverse" mu value, i.e. instead of CT, take TC
			// assuming these values are the same
			if (to == From) continue;
			ll += (table[p][forward][From][to] + table[p][reverse][From][to])*log(mu[to][From]);
			mu_FF -= mu[to][From];
		}
			ll += (table[p][forward][From][From] + table[p][reverse][From][From])*log(mu_FF);
	}
	return ll;
}

template<size_t End, genometools::Base From, genometools::Base To>
double ll_withPMD(const PMDTable &table, const TMu &mu, const TFunction *fun) {
	using genometools::Base;
	constexpr auto forward = makeForward<End>();
	constexpr auto reverse = makeReverse<End>();

	double ll = 0.;
	for (size_t p = 0; p < table.size(); ++p) {
			double mu_FF = 1.;

			for (auto to = Base::min; to < Base::max; ++to) {
			// take "reverse" mu value, i.e. instead of CT, take TC
			// assuming these values are the same
			if (to == From || to == To) continue;
			ll += (table[p][forward][From][to] + table[p][reverse][From][to]) * log(mu[to][From]);
			mu_FF -= mu[to][From];
			}

			const auto mu_FT = mu[To][From] + (1. - mu[To][From]) * fun->prob(p);
			ll += (table[p][forward][From][To] + table[p][reverse][From][To]) * log(mu_FT);
			mu_FF -= mu_FT;

			ll += (table[p][forward][From][From] + table[p][reverse][From][From]) * log(mu_FF);
	}
	return ll;
}

} // namespace impl


template<Strand strand, bool perLength> class TWithPMD final : public TModel {
private:
	static constexpr coretools::TStrongArray<std::string_view, Strand> _names{{"singleStrand", "doubleStrand"}};
	using Function = std::conditional_t<perLength, std::vector<std::unique_ptr<TFunction>>, std::unique_ptr<TFunction>>;
	using Table    = std::conditional_t<perLength, impl::PMDTables, impl::PMDTable>;

	Function _pmd5;
	Function _pmd3;
	Table _table;

	enum class From : bool {from5, from3};
	template<From from>
	coretools::Probability _pmd(const BAM::TSequencedBase &data) const noexcept {
		if constexpr (perLength) {
			if constexpr (from == From::from5) {
				return _pmd5[_table.index(data.fragmentLength)]->prob(data.distFrom5Prime);
			} else {
				return _pmd3[_table.index(data.fragmentLength)]->prob(data.distFrom3Prime);
			}
		} else {
			if constexpr (from == From::from5) {
				return _pmd5->prob(data.distFrom5Prime);
			} else {
				return _pmd3->prob(data.distFrom3Prime);
			}
		}
	}

	coretools::Probability _CT_single(const BAM::TSequencedBase &data) const noexcept {
		const bool isForward = !data.isReverseStrand();
		const bool from5     = data.distFrom5Prime < data.distFrom3Prime;

		if (!isForward) return coretools::Probability{}; 

		// forward:
		if (from5) {
			return _pmd<From::from5>(data);
		} else { // from3
			return _pmd<From::from3>(data);
		}
	}

	coretools::Probability _GA_single(const BAM::TSequencedBase &data) const noexcept {
		const bool isForward = !data.isReverseStrand();
		const bool from5     = data.distFrom5Prime < data.distFrom3Prime;

		if (isForward) return coretools::Probability{}; 

		// reversed:
		if (from5) {
			return _pmd<From::from5>(data);
		} else { // from3
			return _pmd<From::from3>(data);
		}
	}

	coretools::Probability _CT_double(const BAM::TSequencedBase &data) const noexcept {
		const bool isForward = !data.isReverseStrand();
		const bool from5     = data.distFrom5Prime < data.distFrom3Prime;

		if (isForward) {
			if (from5) return _pmd<From::from5>(data);
			return coretools::Probability{}; // from3
		} else { // reversed
			if (from5) return coretools::Probability{};
			return _pmd<From::from3>(data);
		}
	}

	coretools::Probability _GA_double(const BAM::TSequencedBase &data) const noexcept {
		const bool isForward = !data.isReverseStrand();
		const bool from5     = data.distFrom5Prime < data.distFrom3Prime;

		if (isForward) {
			if (from5) return coretools::Probability{};
			return _pmd<From::from3>(data);
		} else { // reversed
			if (from5) return _pmd<From::from5>(data);
			return coretools::Probability{}; // from3
		}
	}

	coretools::Probability _probCT(const BAM::TSequencedBase &data) const noexcept {
		if constexpr (strand == Strand::Single)
			return _CT_single(data);
		else
			return _CT_double(data);
	}

	coretools::Probability _probGA(const BAM::TSequencedBase &data) const noexcept {
		if constexpr (strand == Strand::Single)
			return _GA_single(data);
		else
			return _GA_double(data);
	}

public:
	static constexpr std::string_view name = _names[strand];
	template<typename... Ts>
	TWithPMD(std::string_view function5, std::string_view function3, Ts... ts) {
		constexpr auto N = sizeof...(ts);
		static_assert((perLength && (N == 1)) || (!perLength && !N));
		if constexpr (perLength) {
			_table.minLength = {ts...};
			if (function5.front() != '(') {
				_pmd5.emplace_back(makeFunction(function5));
			} else {
				function5.remove_prefix(1);
				function5.remove_suffix(1);
				coretools::str::TSplitter spl(function5, ';');
				for (auto f5: spl) {
					_pmd5.emplace_back(makeFunction(f5));
				}
			}
			if (function3.front() != '(') {
				_pmd3.emplace_back(makeFunction(function3));
			} else {
				function3.remove_prefix(1);
				function3.remove_suffix(1);
				coretools::str::TSplitter spl(function3, ';');
				for (auto f3: spl) {
					_pmd3.emplace_back(makeFunction(f3));
				}
			}
		} else {
		_pmd5.reset(makeFunction(function5));
		_pmd3.reset(makeFunction(function3));
		}
	}

	bool hasDamage() const noexcept override {
		if constexpr (perLength) {
			return _pmd5.front()->hasDamage() || _pmd3.front()->hasDamage();
		} else {
			return _pmd5->hasDamage() || _pmd3->hasDamage();
		}
	};

	std::string functionString() const noexcept override {
		if constexpr (perLength) {
			std::string s{name};
			s.append(":[");
			s.append(coretools::str::toString(_table.minLength));
			s.append("]:(");
			for (const auto &pmd5 : _pmd5) { s.append(pmd5->string()).append(1, ';'); }
			s.pop_back();
			s.append("):(");
			for (const auto &pmd3 : _pmd3) { s.append(pmd3->string()).append(1, ';'); }
			s.pop_back();
			s.append(")");
			return s;

		} else {
			return std::string(name).append(1, ':').append(_pmd5->string()).append(1, ':').append(_pmd3->string());
		}
	}

	std::string_view typeString() const noexcept override { return name; }

	void resize(size_t N) override {
		if constexpr (perLength) {
			_table.tables.clear();
			_table.tables.emplace_back();
			_table.tables.front().resize(N, {});
		} else {
			_table.resize(N, {});
		}
	}

	void writeTable(std::string_view name, std::array<coretools::TOutputFile, 2> &files) const noexcept override {
		using genometools::Base;
		constexpr auto directions = []() {
			coretools::TStrongArray<std::array<std::string_view, 2>, ReadEnd> ar{};
			ar[ReadEnd::forward3] = {"forward", "3'"};
			ar[ReadEnd::forward5] = {"forward", "5'"};
			ar[ReadEnd::reverse3] = {"reverse", "3'"};
			ar[ReadEnd::reverse5] = {"reverse", "5'"};
			return ar;
		}();

		if constexpr (perLength) {
			for (size_t l = 0; l < _table.tables.size(); ++l) {
				const auto &table = _table.tables[l];
				for (auto j = ReadEnd::min; j < ReadEnd::max; ++j) {
					for (Base f = Base::min; f <= Base::max; ++f) {
						std::vector<size_t> sums(table.size(), 0.);
						for (size_t pos = 0; pos < table.size(); ++pos) {
							for (Base t = Base::min; t < Base::max; ++t) { sums[pos] += table[pos][j][f][t]; }
						}

						for (Base t = Base::min; t <= Base::max; ++t) {
							files.front().writeNoDelim(name, "_");
							files.back().writeNoDelim(name, "_");
							files.front().write(_table.minLength + l, directions[j], f, t);
							files.back().write(_table.minLength + l, directions[j], f, t);
							for (size_t pos = 0; pos < table.size(); ++pos) {
								files.front().write(table[pos][j][f][t]);
								files.back().write(static_cast<double>(table[pos][j][f][t]) / sums[l]);
							}
							files.front().endln();
							files.back().endln();
						}
						files.front().writeNoDelim(name, "_");
						files.front().writeln(_table.minLength + l, directions[j], f, "sum", sums);
					}
				}
			}
		} else {
			for (auto j = ReadEnd::min; j < ReadEnd::max; ++j) {
				for (Base f = Base::min; f <= Base::max; ++f) {
					std::vector<size_t> sums(_table.size(), 0.);
					for (size_t i = 0; i < _table.size(); ++i) {
						for (Base t = Base::min; t < Base::max; ++t) { sums[i] += _table[i][j][f][t]; }
					}
					for (Base t = Base::min; t <= Base::max; ++t) {
						files.front().write(name, directions[j], f, t);
						files.back().write(name, directions[j], f, t);
						for (size_t i = 0; i < _table.size(); ++i) {
							files.front().write(_table[i][j][f][t]);
							files.back().write(static_cast<double>(_table[i][j][f][t]) / sums[i]);
						}
						files.front().endln();
						files.back().endln();
					}
					files.front().writeln(name, directions[j], f, "sum", sums);
				}
			}
		}
	}

	void add(genometools::Base from, BAM::TSequencedBase data) override {
		const auto to    = data.base;
		const auto from3 = data.distFrom3Prime < data.distFrom5Prime;
		if constexpr (perLength) {
			const auto index   = std::max<size_t>(_table.minLength, data.fragmentLength) - _table.minLength;
			const auto oldSize = _table.tables.size();
			if (oldSize <= index) {
				_table.tables.resize(index + 1);
				for (size_t i = oldSize; i < index + 1; ++i) _table.tables[i].resize(_table.tables.front().size(), {});
			}
			const auto pos = std::min<size_t>(_table.tables.size() - 1, from3 ? data.distFrom3Prime : data.distFrom5Prime);
			if (data.isReverseStrand()) {
				const auto readEnd = from3 ? ReadEnd::reverse3 : ReadEnd::reverse5;
				_table.tables[index][pos][readEnd][flipped(from)][flipped(to)]++;
			} else {
				const auto readEnd = from3 ? ReadEnd::forward3 : ReadEnd::forward5;
				_table.tables[index][pos][readEnd][from][to]++;
			}
		} else {
			const auto pos   = std::min<size_t>(_table.size() - 1, from3 ? data.distFrom3Prime : data.distFrom5Prime);
			if (data.isReverseStrand()) {
				const auto readEnd = from3 ? ReadEnd::reverse3 : ReadEnd::reverse5;
				_table[pos][readEnd][flipped(from)][flipped(to)]++;
			} else {
				const auto readEnd = from3 ? ReadEnd::forward3 : ReadEnd::forward5;
				_table[pos][readEnd][from][to]++;
			}
		}
	}

	void estimate() override {
		using coretools::instances::logfile;
		using genometools::Base;
		logfile().startIndent("Learning 5' C-T pattern:");
		if constexpr (perLength) {
			const std::string fun = _pmd5.front()->string();
			_pmd5.clear();
			for (size_t i = 0; i < _table.tables.size(); ++i) {
				const auto [C_T5, T_C5] = impl::makeFromTo<5, Base::C, Base::T>(_table.tables[i]);
				_pmd5.emplace_back(makeFunction(fun));
				_pmd5.back()->learn(C_T5, T_C5);
			}
		} else {
			const auto [C_T5, T_C5] = impl::makeFromTo<5, Base::C, Base::T>(_table);
			_pmd5->learn(C_T5, T_C5);
			const auto mu_5        = impl::makeMu<5>(_table);
			const double ll_no_5   = impl::ll_noPMD<5, Base::C>(_table, mu_5);
			const double ll_with_5 = impl::ll_withPMD<5, Base::C, Base::T>(_table, mu_5, _pmd5.get());
			OUT(ll_no_5);
			OUT(ll_with_5);
		}
		logfile().endIndent();

		constexpr auto from = [](){
			if constexpr (strand == Strand::Single) return Base::C;
			else return Base::G;
		}();
		constexpr auto to = []() {
			if constexpr (strand == Strand::Single)
				return Base::T;
			else
				return Base::A;
		}();

		logfile().startIndent("Learning 3' ", from, "-", to, " pattern:");
		if constexpr (perLength) {
			const std::string fun = _pmd3.front()->string();
			_pmd3.clear();
			for (size_t i = 0; i < _table.tables.size(); ++i) {
				const auto [from_to, to_from] = impl::makeFromTo<3, from, to>(_table.tables[i]);
				_pmd3.emplace_back(makeFunction(fun));
				_pmd3[i]->learn(from_to, to_from);
			}
		} else {
			const auto [from_to, to_from] = impl::makeFromTo<3, from, to>(_table);
			_pmd3->learn(from_to, to_from);

			const auto mu_3        = impl::makeMu<3>(_table);
			const double ll_no_3   = impl::ll_noPMD<3, from>(_table, mu_3);
			const double ll_with_3 = impl::ll_withPMD<3, from, to>(_table, mu_3, _pmd3.get());
			OUT(ll_no_3);
			OUT(ll_with_3);
		}
		logfile().endIndent();
	}

	TBaseLikelihoods baseLikelihoods(const BAM::TSequencedBase &data,
									 const TBaseLikelihoods &baseLikelihoodsNoPMD) const override {
		using genometools::Base;
		const auto pCT = _probCT(data);
		const auto pGA = _probGA(data);
		TBaseLikelihoods baseLikelihoods(baseLikelihoodsNoPMD);

		baseLikelihoods[Base::C] = (1.0 - pCT) * baseLikelihoodsNoPMD[Base::C] + pCT * baseLikelihoodsNoPMD[Base::T];
		baseLikelihoods[Base::G] = (1.0 - pGA) * baseLikelihoodsNoPMD[Base::G] + pGA * baseLikelihoodsNoPMD[Base::A];
		return baseLikelihoods;
	}

	TBaseProbabilities massFunction(genometools::Base b, const BAM::TSequencedBase &data,
									const TBaseLikelihoods &baseLikelihoodsNoPMD) const override {
		using genometools::Base;
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
		using genometools::Base;
		using coretools::instances::randomGenerator;
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

} // namespace GenotypeLikelihoods

#endif
