#include "TPsi.h"

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Strings/stringManipulations.h"
#include "coretools/Types/probability.h"
#include <cmath>
#include <limits>

namespace GenotypeLikelihoods::PMD {
using coretools::Probability;
using coretools::P;
using coretools::index;
using BAM::End;
using coretools::str::fromString;
using coretools::str::toString;
using coretools::str::readAfter;
using coretools::str::readBefore;

namespace impl {
std::vector<Probability> exp(std::string_view sFunction, size_t N = 50) {
	const auto a    = fromString<double, true>(readBefore(sFunction, '*'));
	const auto rest = readAfter(sFunction, "exp(");
	const auto b    = fromString<double, true>(readBefore(rest, '*'));
	const auto c    = fromString<double, true>(readAfter(rest, '+'));
	std::vector<Probability> probs;
	for (size_t i = 0; i < N; ++i) {
		probs.emplace_back(a * std::exp(b * i) + c);
	}
	return probs;
}

std::vector<Probability> empiric(std::string_view sFunction) {
	coretools::str::TSplitter spl(sFunction, ',');
	if (spl.empty()) return {P(0.)};

	std::vector<Probability> probs;
	for (auto s: spl) {
		probs.emplace_back(fromString<double, true>(s));
	}
	return probs;
}

Type type(std::string_view sType) {
	if (sType == "CT") return Type::CT;
	if (sType == "GA") return Type::GA;
	UERROR("Type ", sType, " is not a recognized token!");
}

End end(std::string_view sEnd) {
	if (sEnd == "5") return End::from5;
	if (sEnd == "3") return End::from3;
	UERROR("End ", sEnd, " is not a recognized token!");
}

std::string toString(Type type) {
	if (type == Type::CT) return "CT";
	else return "GA";
}

} // namespace impl

void TPsi::_fromString(std::string_view Psi) {
	// CT5:0.1,0.09;GA3:0.3*exp()
	coretools::str::TSplitter semi(Psi, ';');
	for (auto s: semi) {
		const auto sType     = readBefore(s, ':');
		if (sType.size() != 3) UERROR(sType, " is not a recognized token!");

		const auto type      = impl::type(sType.substr(0, 2));
		const auto end       = index(impl::end(sType.substr(2, 1)));
		const auto sFunction = readAfter(s, ':');
		if (sFunction.find('*') != sFunction.npos) _tables[end][type] = impl::exp(sFunction);
		else _tables[end][type] = impl::empiric(sFunction);
	}
	for (auto& t: _tables)
		for (auto &v: t)
			if (v.empty()) v = {P(0.)};
}

	TPsi::TPsi(const BAM::RGInfo::TInfo &Info): _tables(2) {
	_parse(Info);
}

	TPsi::TPsi(): _tables(2), _tableSums(2) {
	_tables.resize(2);
	_tableSums.resize(2);
	for (auto &t : _tables)
		for (auto &v : t) v.emplace_back(0.);
}

coretools::Probability TPsi::prob(Type From_To, const BAM::TSequencedData &data) const noexcept {
	using BAM::End;
	const auto realType = data.get<BAM::Flags::ReversedStrand>() ? _flip[From_To] : From_To;
	const auto end      = data.get<BAM::Flags::Paired>() ? End::from5 : data.end();
	const auto pos      = data.dist(end).pseudo();

	const auto &table = _tables[index(end)][realType];
	if (pos >= table.size()) return table.back();
	return table[pos];
}

void TPsi::add(Type From_To, const BAM::TSequencedData &data, coretools::Probability P_g_I_di,
			   const TBaseBaseProbabilities &P_b_bbar_I_gdij) {
	using BAM::End;
	using genometools::Base;

	const auto end      = paired() ? End::from5 : data.end();
	const auto realType = data.get<BAM::Flags::ReversedStrand>() ? _flip[From_To] : From_To;
	auto &tSum          = _tableSums[index(end)][realType];
	if (tSum.empty()) return; // wrong pattern

	const auto pos  = std::min<size_t>(tSum.size() - 1, data.dist(end).pseudo());
	const auto from = _from[From_To];
	const auto to   = _to[From_To];
	tSum[pos].numDenom.num += P_b_bbar_I_gdij[from][to] * P_g_I_di;
	tSum[pos].numDenom.denom += (P_b_bbar_I_gdij[from][to] + P_b_bbar_I_gdij[from][from]) * P_g_I_di;
}

void TPsi::_add(Type From_To, const BAM::TSequencedData &data, genometools::Base ref) {
	using BAM::End;
	using genometools::Base;

	const auto end  = data.end();
	const auto pos  = data.dist(end).pseudo();
	const auto base = data.base;

	const auto realType = data.get<BAM::Flags::ReversedStrand>() ? _flip[From_To] : From_To;
	auto &tSum          = _tableSums[index(end)][realType];
	if (tSum.size() <= pos) tSum.resize(pos + 1, {{{0, 0, 0, 0}}});

	const auto From = _from[From_To];
	const auto To   = _to[From_To];
	if (ref == From) {
		if (base == To) tSum[pos].fromTo.fromTo += 1.;
		tSum[pos].fromTo.fromSum += 1.;
	} else if (ref == To) {
		if (base == From) tSum[pos].fromTo.toFrom += 1.;
		tSum[pos].fromTo.toSum += 1.;
	}
}

void TPsi::_parse(const BAM::RGInfo::TInfo & Info) {
	if (Info.is_string()) {
		_fromString(Info.get<std::string_view>());
	} else {
		for (auto e = End::min; e < End::max; ++e) {
			for (auto t = Type::min; t < Type::max; ++t) {
				auto &v        = _tables[index(e)][t];
				const auto key = impl::toString(t) + toString(e);
				if (!Info.contains(key) || Info[key].empty()) {
					v = {P(0.)};
				} else {
					if (Info[key].is_string()) {
						const auto sFunction = Info[key].get<std::string_view>();
						if (sFunction.find('*') != sFunction.npos)
							v = impl::exp(sFunction);
						else
							v = impl::empiric(sFunction);
					} else if (Info[key].is_array()) {
						for (auto [k, d] : Info[key].items()) {
							v.emplace_back(d);
						}
					} else {
						UERROR("Cannot parse json-token ", Info, "!");
					}
				}
			}
		}
	}
}

BAM::RGInfo::TInfo TPsi::info() const {
	BAM::RGInfo::TInfo info;
	for (auto e = End::min; e < End::max; ++e) {
		for (auto t = Type::min; t < Type::max; ++t) {
			const auto &v = _tables[index(e)][t];
			if (v.size() > 1 || v.front() > 0.) {
				const auto key = impl::toString(t) + toString(e);
				for (auto p : v) info[key].push_back(p.get()); // explicit conversion from probability to double
			}
		}
	}
	return info;
}

void TPsi::estimate() noexcept {
	constexpr double PMDmin = 1e-9;
	for (auto e = End::min; e < End::max; ++e) {
		//const auto t = type(e);
		for (auto t = Type::min; t < Type::max; ++t) {
			auto &table  = _tables[index(e)][t];
			auto &tSum   = _tableSums[index(e)][t];

			table.clear();
			for (size_t i = 0; i < tSum.size(); ++i) {
				auto &ts       = tSum[i];
				const auto PMD = std::max(PMDmin, ts.numDenom.num / ts.numDenom.denom);
				table.emplace_back(PMD);
				ts.numDenom.num   = 0.;
				ts.numDenom.denom = std::numeric_limits<double>::min(); // preventing any division by 0
			}
		}
	}
}

void TPsi::_printTable(std::string_view FName) {
	coretools::TOutputFile oFile(FName);
	for (auto e = End::min; e < End::max; ++e) {
		for (auto t = Type::min; t < Type::max; ++t) {
			oFile.write(impl::toString(t) + toString(e));
			auto &tSum = _tableSums[index(e)][t];
			for (const auto ts: tSum) {
				const auto fromTo = double(ts.fromTo.fromTo) / ts.fromTo.fromSum;
				const auto toFrom = double(ts.fromTo.toFrom) / ts.fromTo.toSum;
				const auto PMD    = std::max(0., (fromTo - toFrom) / (1.0 - toFrom));
				oFile.writeNoDelim(ts.fromTo.fromTo, "/", ts.fromTo.fromSum, ";", ts.fromTo.toFrom, "/", ts.fromTo.toSum, ":", PMD).writeDelim();
			}
			oFile.endln();
		}
	}
}

void TPsi::_initEnd(End e, int32_t MinData) {
	constexpr double PMDmin = 1e-9;

	coretools::TStrongArray<double, Type> sums{};
	for (auto t = Type::min; t < Type::max; ++t) {
		auto &table = _tables[index(e)][t];
		table.clear();

		auto &tSum = _tableSums[index(e)][t];
		if (tSum.empty()) continue;

		while (tSum.size() > 1) {
			const auto &ts = tSum.back();
			auto &ts_m     = *(tSum.end() - 2);

			if ((ts.fromTo.fromSum < MinData || ts.fromTo.toSum < MinData) ||
				(ts_m.fromTo.fromSum < MinData || ts_m.fromTo.toSum < MinData))
			{
				ts_m.fromTo.fromSum += ts.fromTo.fromSum;
				ts_m.fromTo.fromTo += ts.fromTo.fromTo;
				ts_m.fromTo.toSum += ts.fromTo.toSum;
				ts_m.fromTo.toFrom += ts.fromTo.toFrom;
				tSum.pop_back();
			} else {
				break;
			}
		}

		for (auto &ts : tSum) {
			const auto fromTo = double(ts.fromTo.fromTo) / ts.fromTo.fromSum;
			const auto toFrom = double(ts.fromTo.toFrom) / ts.fromTo.toSum;
			const auto PMD    = std::max(PMDmin, (fromTo - toFrom) / (1.0 - toFrom));
			// once 0, always 0, so add something small

			table.emplace_back(PMD);
			sums[t] += table.back();

			// Prepare ts for probabilistic sum
			ts.numDenom.num   = 0.;
			ts.numDenom.denom = std::numeric_limits<double>::min(); // preventing any division by 0
		}
	}
}

void TPsi::_joinTables() noexcept {
	for (auto t = Type::min; t < Type::max; ++t) {
		auto &ts5 = _tableSums[0][t];
		auto &ts3 = _tableSums[1][t];

		if (ts3.size() > ts5.size())
			ts5.resize(ts3.size());

		for (size_t i = 0; i < ts3.size(); ++i) {
			// if ts5.size() > ts3.size(), we only go to ts3.size()
			ts5[i].fromTo.fromTo += ts3[i].fromTo.fromTo;
			ts5[i].fromTo.fromSum += ts3[i].fromTo.fromSum;
			ts5[i].fromTo.toFrom += ts3[i].fromTo.toFrom;
			ts5[i].fromTo.toSum += ts3[i].fromTo.toSum;
		}
		ts3.clear();
	}
}

void TPsi::estimateInit(std::string_view OutputName, size_t MinData) noexcept {
	using coretools::instances::logfile;
	const auto fn = toString(OutputName, "_PsiTable.txt.gz");
	logfile().list("Writing countTable '", fn, "'.");
	_printTable(toString(OutputName, "_PsiTable.txt.gz"));

	if (_nPaired > 0 && _nSingle > 0) logfile().warning("Readgroup contains both single- and paired-ended reads!");
	if (paired()) {
		logfile().list("Assuming paired-ended reads.");

		// only 5' end
		_joinTables();
		_tables[1][Type::CT] = {P(0.)};
		_tables[1][Type::GA] = {P(0.)};
	} else {
		logfile().list("Assuming single-ended reads.");
		_initEnd(End::from3, MinData);
	}
	_initEnd(End::from5, MinData);
}

void TPsi::log() const noexcept {
	using coretools::instances::logfile;
	constexpr size_t Nmax = 3;

	bool hasAny = false;
	for (auto e = End::min; e < End::max; ++e) {
		if (paired() && e == End::from3) continue;

		for (auto t = Type::min; t < Type::max; ++t) {
			const auto &v = _tables[index(e)][t];
			if (v.size() > 1 || v.front() > 0.) {
				hasAny = true;
				auto ret = impl::toString(t) + toString(e) + ": [";
				if (v.size() <= Nmax*2) {
					for (auto p : v) ret.append(toString(p, ", "));
				} else {
					for (size_t i = 0; i < Nmax; ++i) ret.append(toString(v[i] , ", "));
					ret.append("..., ");
					const auto iStart = v.size() - 1 - Nmax;
					for (size_t i = 0; i < Nmax; ++i) ret.append(toString(v[iStart + i], ", "));
				}
				ret.pop_back();
				ret.back() = ']';
				logfile().list(ret);
			}
		}
	}
	if (!hasAny) logfile().list("[]");
}

} // namespace GenotypeLikelihoods::PMD
