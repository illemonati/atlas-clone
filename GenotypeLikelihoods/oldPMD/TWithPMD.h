/*
 * TReadGroupPMD.h
 *
 *  Created on: Mars 31, 2023
 *      Author: andreas
 */

#ifndef PMD_TPERREADGROUP_H_
#define PMD_TPERREADGROUP_H_

#include <memory>
#include <type_traits>
#include <vector>

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"

#include "TFunction.h"
#include "TAlignment.h"
#include "TModel.h"

namespace GenotypeLikelihoods::oldPMD {
namespace impl {

static constexpr size_t _N = coretools::index(genometools::Base::max) + 1;
using PMDTable             = std::vector<coretools::TStrongArray<
	coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base, _N>, genometools::Base, _N>, ReadEnd>>;
using TMu =
	coretools::TStrongArray<coretools::TStrongArray<coretools::Probability, genometools::Base>, genometools::Base>;

struct PMDTables {
	std::vector<PMDTable> tables;
	size_t minLength;
	size_t maxSize;
	size_t index(size_t length) const noexcept {
		return std::min(std::max<size_t>(minLength, length) - minLength, tables.size() - 1);
	}
};

template<size_t End>
constexpr ReadEnd makeForward() {
	if constexpr (End == 3) return ReadEnd::forward3;
	if constexpr (End == 5) return ReadEnd::forward5;
	else static_assert(End == 3, "compile-time error");
}

template<size_t End>
constexpr ReadEnd makeReverse() {
	if constexpr (End == 3) return ReadEnd::reverse3;
	if constexpr (End == 5) return ReadEnd::reverse5;
	else static_assert(End == 3, "compile-time error");
}

template<size_t End>
auto makeFromTo(const PMDTable &table, genometools::Base _from, genometools::Base _to) {
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
		from_to.push_back(table[i][forward][_from][_to]);
		from_to.back() += table[i][reverse][_from][_to];
		to_from.push_back(table[i][forward][_to][_from]);
		to_from.back() += table[i][reverse][_to][_from];
		size_t s_from = 0;
		size_t s_to   = 0;
		for (auto b = Base::min; b < Base::max; ++b) {
			s_from += table[i][forward][_from][b];
			s_from += table[i][reverse][_from][b];
			s_to += table[i][forward][_to][b];
			s_to += table[i][reverse][_to][b];
		}
		from_to.back() /= s_from;
		to_from.back() /= s_to;

		if (s_from < 1000 || (from_to.back() - to_from.back()) < 0.01) {
			from_to.pop_back();
			to_from.pop_back();
			break;
		}

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
		size_t s = 0;
		for (auto to = Base::min; to < Base::max; ++to) {
			s += from_to[from][to];
		}
		for (auto to = Base::min; to < Base::max; ++to) {
			mu[from][to] = static_cast<double>(from_to[from][to])/s;
		}
	}
	return mu;
}

template<size_t End>
double ll_noPMD(const PMDTable &table, const TMu& mu, genometools::Base From) {
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

template<size_t End>
double ll_withPMD(const PMDTable &table, const TMu &mu, const TFunction *fun, genometools::Base From, genometools::Base To) {
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


template<bool perLength> class TWithPMD final : public TModel {
private:
	static constexpr coretools::TStrongArray<std::string_view, Strand, coretools::index(Strand::Unknown) + 1> _names{{"singleStrand", "doubleStrand", "unknownStrand"}};
	using Function = std::conditional_t<perLength, std::vector<std::unique_ptr<TFunction>>, std::unique_ptr<TFunction>>;
	using Table    = std::conditional_t<perLength, impl::PMDTables, impl::PMDTable>;

	Strand _strand;
	Function _pmd5;
	Function _pmd3;
	Table _table;

	template<size_t End>
	coretools::Probability _pmd(const BAM::TSequencedBase &data) const noexcept {
		if constexpr (perLength) {
			if constexpr (End == 5) {
				return _pmd5[_table.index(data.fragmentLength)]->prob(data.distFrom5Prime);
			} else if constexpr (End == 3){
				return _pmd3[_table.index(data.fragmentLength)]->prob(data.distFrom3Prime);
			} else {
				static_assert(End == 3, "compile-time error");
			}
		} else {
			if constexpr (End == 5) {
				return _pmd5->prob(data.distFrom5Prime);
			} else if constexpr (End == 3){
				return _pmd3->prob(data.distFrom3Prime);
			} else {
				static_assert(End == 3, "compile-time error");
			}
		}
	}

	template<Strand strand> coretools::Probability _CT(const BAM::TSequencedBase &data) const noexcept {
		if constexpr (strand == Strand::Single) {
			const bool isForward = !data.isReverseStrand();
			if (!isForward) return coretools::Probability{};

			// forward:
			const bool from5 = data.distFrom5Prime < data.distFrom3Prime;
			return from5 ? _pmd<5>(data) : _pmd<3>(data);

		} else if constexpr (strand == Strand::Double) {
			const bool isForward = !data.isReverseStrand();
			const bool from5     = data.distFrom5Prime < data.distFrom3Prime;

			if (isForward) {
				if (from5) return _pmd<5>(data);
				return coretools::Probability{}; // from3
			} else {                             // reversed
				if (from5) return coretools::Probability{};
				return _pmd<3>(data);
			}
		} else {
			return coretools::Probability{};
		}
	}

	template<Strand strand> coretools::Probability _GA(const BAM::TSequencedBase &data) const noexcept {
		if constexpr (strand == Strand::Single) {
			const bool isForward = !data.isReverseStrand();
			if (isForward) return coretools::Probability{};

			// reversed:
			const bool from5 = data.distFrom5Prime < data.distFrom3Prime;
			return from5 ? _pmd<5>(data) : _pmd<3>(data);

		} else if constexpr (strand == Strand::Double) {
			const bool isForward = !data.isReverseStrand();
			const bool from5     = data.distFrom5Prime < data.distFrom3Prime;

			if (isForward) {
				if (from5) return coretools::Probability{};
				return _pmd<3>(data);
			} else { // reversed
				if (from5) return _pmd<5>(data);
				return coretools::Probability{}; // from3
			}

		} else {
			return coretools::Probability{};
		}
	}

	coretools::Probability _probCT(const BAM::TSequencedBase &data) const noexcept {
		switch (_strand) {
		case Strand::Single: return _CT<Strand::Single>(data);
		case Strand::Double: return _CT<Strand::Double>(data);
		default: return _CT<Strand::Unknown>(data);
		}
	}

	coretools::Probability _probGA(const BAM::TSequencedBase &data) const noexcept {
		switch (_strand) {
		case Strand::Single: return _GA<Strand::Single>(data);
		case Strand::Double: return _GA<Strand::Double>(data);
		default: return _GA<Strand::Unknown>(data);
		}
	}

public:
	template<typename... Ts>
	TWithPMD(std::string_view function5, std::string_view function3, Strand strand, Ts... ts) : _strand(strand) {
		constexpr auto N = sizeof...(ts);
		if constexpr (perLength) {
			static_assert(N == 1);
			static_assert(std::is_same_v<Ts..., size_t>);

			_table.minLength = {ts...};
			if (function5.front() != '(') {
				_pmd5.emplace_back(makeFunction(function5));
			} else {
				function5.remove_prefix(1);
				function5.remove_suffix(1);
				coretools::str::TSplitter spl(function5, ';');
				for (auto f5 : spl) { _pmd5.emplace_back(makeFunction(f5)); }
			}
			if (function3.front() != '(') {
				_pmd3.emplace_back(makeFunction(function3));
			} else {
				function3.remove_prefix(1);
				function3.remove_suffix(1);
				coretools::str::TSplitter spl(function3, ';');
				for (auto f3 : spl) { _pmd3.emplace_back(makeFunction(f3)); }
			}
		} else {
			static_assert(N == 0);
			_pmd5.reset(makeFunction(function5));
			_pmd3.reset(makeFunction(function3));
		}
	}

	template<typename... Ts>
	TWithPMD(std::string_view function5, Ts... ts): TWithPMD(function5, function5, Strand::Unknown, ts...) {}


	bool hasDamage() const noexcept override {
		if constexpr (perLength) {
			return _pmd5.front()->hasDamage() || _pmd3.front()->hasDamage();
		} else {
			return _pmd5->hasDamage() || _pmd3->hasDamage();
		}
	};

	std::string functionString() const noexcept override {
		if constexpr (perLength) {
			std::string s{_names[_strand]};
			s.append("[");
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
			return std::string(_names[_strand]).append(1, ':').append(_pmd5->string()).append(1, ':').append(_pmd3->string());
		}
	}

	std::string_view typeString() const noexcept override { return _names[_strand]; }

	void resize(size_t N) override {
		if constexpr (perLength) {
			_table.maxSize = N;
			for (size_t i = 0; i < _table.tables.size(); ++i) {
				const auto sz = std::min(N, (_table.minLength + i - 1)/2 + 1); // this always ceils
				_table.tables[i].resize(sz, {});
			}
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
						std::vector<size_t> sums(_table.maxSize, 0.);
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
								files.back().write(static_cast<double>(table[pos][j][f][t]) / sums[pos]);
							}
							for (size_t pos = table.size(); pos < _table.maxSize; ++ pos) {
								files.front().write(0.);
								files.back().write(0.);
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
		const bool from5 = data.distFrom5Prime < data.distFrom3Prime;
		if constexpr (perLength) {
			const auto index   = std::max<size_t>(_table.minLength, data.fragmentLength) - _table.minLength;
			const auto oldSize = _table.tables.size();
			if (oldSize <= index) {
				_table.tables.resize(index + 1);
				for (size_t i = oldSize; i < index + 1; ++i) {
					const auto sz = (_table.minLength + i - 1) / 2 + 1; // this always ceils
					_table.tables[i].resize(sz, {});
				}
			}
			const auto pos = std::min<size_t>(_table.tables.size() - 1, from5 ? data.distFrom5Prime : data.distFrom3Prime);
			if (data.isReverseStrand()) {
				const auto readEnd = from5 ? ReadEnd::reverse5 : ReadEnd::reverse3;
				_table.tables[index][pos][readEnd][flipped(from)][flipped(to)]++;
			} else {
				const auto readEnd = from5 ? ReadEnd::forward5 : ReadEnd::forward3;
				_table.tables[index][pos][readEnd][from][to]++;
			}
		} else {
			const auto pos   = std::min<size_t>(_table.size() - 1, from5 ? data.distFrom5Prime : data.distFrom3Prime);
			if (data.isReverseStrand()) {
				const auto readEnd = from5 ? ReadEnd::reverse5 : ReadEnd::reverse3;
				_table[pos][readEnd][flipped(from)][flipped(to)]++;
			} else {
				const auto readEnd = from5 ? ReadEnd::forward5 : ReadEnd::forward3;
				_table[pos][readEnd][from][to]++;
			}
		}
	}

	void estimate() override {
		using coretools::instances::logfile;
		using genometools::Base;
		logfile().startIndent("Learning 5' C-T pattern:");
		if constexpr (perLength) {
			std::unique_ptr<TFunction> fun{_pmd5.front()->clone()};
			_pmd5.clear();
			for (size_t i = 0; i < _table.tables.size(); ++i) {
				const auto [C_T5, T_C5] = impl::makeFromTo<5>(_table.tables[i], Base::C, Base::T);
				fun->learn(C_T5, T_C5);
				_pmd5.emplace_back(fun->clone());
			}
		} else {
			const auto [C_T5, T_C5] = impl::makeFromTo<5>(_table, Base::C, Base::T);
			_pmd5->learn(C_T5, T_C5);
		}
		logfile().endIndent();

		constexpr auto from_tos = []() {
			coretools::TStrongArray<std::array<Base, 2>, Strand> as{};
			as[Strand::Single] = {Base::C, Base::T};
			as[Strand::Double] = {Base::G, Base::A};
			return as;
		}();
		coretools::TStrongArray<Function, Strand> funs;
		coretools::TStrongArray<double, Strand> ll_no;
		coretools::TStrongArray<double, Strand> ll_with;
		for (auto s = Strand::min; s < Strand::max; ++s) {
			const auto from = from_tos[s].front();
			const auto to   = from_tos[s].back();

			logfile().startIndent("Learning 3' ", from, "-", to, " pattern:");
			if constexpr (perLength) {
				for (size_t i = 0; i < _table.tables.size(); ++i) {
					const auto [from_to, to_from] = impl::makeFromTo<3>(_table.tables[i], from, to);
					funs[s].emplace_back(_pmd3.front()->clone());
					funs[s][i]->learn(from_to, to_from);
				}
			} else {
				const auto mu_3 = impl::makeMu<3>(_table);
				funs[s].reset(_pmd3->clone());
				const auto [from_to, to_from] = impl::makeFromTo<3>(_table, from, to);
				funs[s]->learn(from_to, to_from);

				ll_no[s]   = impl::ll_noPMD<3>(_table, mu_3, from);
				ll_with[s] = impl::ll_withPMD<3>(_table, mu_3, funs[s].get(), from, to);
			}
			logfile().endIndent();
		}
		coretools::TStrongArray<double, Strand> ll_tot;
		ll_tot[Strand::Single] = ll_with[Strand::Single] + ll_no[Strand::Double];
		ll_tot[Strand::Double] = ll_with[Strand::Double] + ll_no[Strand::Single];
		logfile().list("Log-Likelihood of single-strand library: ", ll_tot[Strand::Single]);
		logfile().list("Log-Likelihood of double-strand library: ", ll_tot[Strand::Double]);
		const auto delta = ll_tot[Strand::Single] - ll_tot[Strand::Double];
		if (delta > 0) {
			_strand = Strand::Single;
			logfile().list("Assuming single-strand library with delta-log-likelihood: ", delta);
		} else {
			_strand = Strand::Double;
			logfile().list("Assuming double-strand library with delta-log-likelihood: ", -delta);
		}
		std::swap(_pmd3, funs[_strand]);
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

	void simulate(BAM::TAlignment &aln) const override {
		using coretools::instances::randomGenerator;
		using genometools::Base;
		for (auto &d : aln) {
			const auto pCT = _probCT(d);
			const auto pGA = _probGA(d);
			auto &base     = d.base;

			if (base == Base::C) {
				if (randomGenerator().getRand() < pCT) base = Base::T;
			} else if (base == Base::G) {
				if (randomGenerator().getRand() < pGA) base = Base::A;
			}
		}
	}
};

} // namespace GenotypeLikelihoods

#endif
