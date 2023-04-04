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
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
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

struct PMDTables {
	std::vector<PMDTable> tables;
	size_t minLength;
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

} // namespace impl


template<Strand strand, bool perLength> class TWithPMD final : public TModel {
private:
	static constexpr coretools::TStrongArray<std::string_view, Strand> _names{{"singleStrand", "doubleStrand"}};
	using Function = std::conditional_t<perLength, std::vector<std::unique_ptr<TFunction>>, std::unique_ptr<TFunction>>;
	using Table    = std::conditional_t<perLength, impl::PMDTables, impl::PMDTable>;
	Function _pmd5;
	Function _pmd3;
	Table _table;

	coretools::Probability _probCT(const BAM::TSequencedBase &data) const noexcept {
		using coretools::Probability;
		if constexpr (strand == Strand::Single) {
			if (data.isReverseStrand()) return coretools::Probability{};
			if constexpr (perLength) {
				return data.distFrom3Prime < data.distFrom5Prime
						   ? _pmd3[data.fragmentLength]->prob(data.distFrom3Prime)
						   : _pmd5[data.fragmentLength]->prob(data.distFrom5Prime);

			} else {
				return data.distFrom3Prime < data.distFrom5Prime ? _pmd3->prob(data.distFrom3Prime)
																 : _pmd5->prob(data.distFrom5Prime);
			}
		} else {
			if (data.distFrom3Prime < data.distFrom5Prime) { // from 3
				if constexpr (perLength) {
					return !data.isReverseStrand() ? Probability{}
												   : _pmd3[data.fragmentLength]->prob(data.distFrom3Prime);
				} else {
					return !data.isReverseStrand() ? Probability{} : _pmd3->prob(data.distFrom3Prime);
				}
			} else { // from 5
				if constexpr (perLength) {
					return !data.isReverseStrand() ? _pmd5[data.fragmentLength]->prob(data.distFrom5Prime)
												   : Probability{};
				} else {
					return !data.isReverseStrand() ? _pmd5->prob(data.distFrom5Prime) : Probability{};
				}
			}
		}
	}

	coretools::Probability _probGA(const BAM::TSequencedBase &data) const noexcept {
		using coretools::Probability;
		if constexpr (strand == Strand::Single) {
			if (!data.isReverseStrand()) return coretools::Probability{};
			if constexpr (perLength) {
				return data.distFrom3Prime < data.distFrom5Prime
						   ? _pmd3[data.fragmentLength]->prob(data.distFrom3Prime)
						   : _pmd5[data.fragmentLength]->prob(data.distFrom5Prime);
			} else {
				return data.distFrom3Prime < data.distFrom5Prime ? _pmd3->prob(data.distFrom3Prime)
																 : _pmd5->prob(data.distFrom5Prime);
			}
		} else {
			if (data.distFrom3Prime < data.distFrom5Prime) { // from 3
				if constexpr (perLength) {
					return !data.isReverseStrand() ? _pmd3[data.fragmentLength]->prob(data.distFrom3Prime)
												   : Probability{};
				} else {
					return !data.isReverseStrand() ? _pmd3->prob(data.distFrom3Prime) : Probability{};
				}
			} else { // from 5
				if constexpr (perLength) {
					return !data.isReverseStrand() ? Probability{}
												   : _pmd5[data.fragmentLength]->prob(data.distFrom5Prime);
				} else {
					return !data.isReverseStrand() ? Probability{} : _pmd5->prob(data.distFrom5Prime);
				}
			}
		}
	}

public:
	static constexpr std::string_view name = _names[strand];
	template<typename... Ts>
	TWithPMD(std::string_view function5, std::string_view function3, Ts... ts) {
		constexpr auto N = sizeof...(ts);
		static_assert((perLength && (N == 2)) || (!perLength && !N));
		if constexpr (perLength) {
			auto&& tpl = std::forward_as_tuple(ts...);
			_table.minLength = std::get<0>(tpl);
			for (size_t i = _table.minLength; i <= std::get<1>(tpl); ++i) {
				_pmd5.emplace_back(makeFunction(function5));
				_pmd3.emplace_back(makeFunction(function3));
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
			for (const auto &pmd5 : _pmd5) { s.append(pmd5->string()).append(1, ','); }
			s.back() = ':';
			for (const auto &pmd3 : _pmd3) { s.append(pmd3->string()).append(1, ','); }
			s.back() = ']';
			return s;

		} else {
			return std::string(name).append(1, ':').append(_pmd5->string()).append(1, ':').append(_pmd3->string());
		}
	}

	std::string_view typeString() const noexcept override { return name; }

	void resize(size_t N) override {
		if constexpr (perLength) {
			_table.tables.resize(_pmd5.size());
			for (auto &table : _table.tables) table.resize(N, {});
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
			const auto len =
				std::clamp<size_t>(data.fragmentLength, _table.minLength, _table.minLength + _table.tables.size() - 1) -
				_table.minLength;
			const auto pos   = std::min<size_t>(_table.tables.size() - 1, from3 ? data.distFrom3Prime : data.distFrom5Prime);
			if (data.isReverseStrand()) {
				const auto readEnd = from3 ? ReadEnd::reverse3 : ReadEnd::reverse5;
				_table.tables[len][pos][readEnd][flipped(from)][flipped(to)]++;
			} else {
				const auto readEnd = from3 ? ReadEnd::forward3 : ReadEnd::forward5;
				_table.tables[len][pos][readEnd][from][to]++;
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
			for (size_t i = 0; i < _pmd5.size(); ++i) {
				const auto [C_T5, T_C5] = impl::makeFromTo<5, Base::C, Base::T>(_table.tables[i]);
				_pmd5[i]->learn(C_T5, T_C5);
			}
		} else {
			const auto [C_T5, T_C5] = impl::makeFromTo<5, Base::C, Base::T>(_table);
			_pmd5->learn(C_T5, T_C5);
		}
		logfile().endIndent();
		if constexpr (strand == Strand::Single) {
			logfile().startIndent("Learning 3' C-T pattern:");
			if constexpr (perLength) {
				for (size_t i = 0; i < _pmd3.size(); ++i) {
					const auto [C_T3, T_C3] = impl::makeFromTo<3, Base::C, Base::T>(_table.tables[i]);
					_pmd3[i]->learn(C_T3, T_C3);
				}
			} else {
				const auto [C_T3, T_C3] = impl::makeFromTo<3, Base::C, Base::T>(_table);
				_pmd3->learn(C_T3, T_C3);
			}
			logfile().endIndent();
		} else {
			logfile().startIndent("Learning 3' G-A pattern:");
			if constexpr (perLength) {
				for (size_t i = 0; i < _pmd3.size(); ++i) {
					const auto [G_A3, A_G3] = impl::makeFromTo<3, Base::G, Base::A>(_table.tables[i]);
					_pmd3[i]->learn(G_A3, A_G3);
				}
			} else {
				const auto [G_A3, A_G3] = impl::makeFromTo<3, Base::G, Base::A>(_table);
				_pmd3->learn(G_A3, A_G3);
			}
			logfile().endIndent();
		}
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
