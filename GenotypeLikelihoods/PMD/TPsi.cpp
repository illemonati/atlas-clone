#include "TPsi.h"

#include <cmath>
#include <limits>

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Strings/stringManipulations.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/probability.h"

#include "TReadGroupInfo.h"
#include "TSequencedData.h"

#include "coretools/devtools.h"


namespace GenotypeLikelihoods::PMD {

using BAM::End;
using coretools::P;
using coretools::Probability;
using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::user_assert;
using coretools::str::fromString;
using coretools::str::readAfter;
using coretools::str::readBefore;
using coretools::str::toString;

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
	throw coretools::TUserError("Type ", sType, " is not a recognized token!");
}

size_t end(std::string_view sEnd) {
	if (sEnd == "5") return 0;
	if (sEnd == "3") return 1;
	throw coretools::TUserError("End ", sEnd, " is not a recognized token!");
}

std::string toString(Type type) {
	if (type == Type::CT) return "CT";
	else return "GA";
}

} // namespace impl

TPsi::TPsi() : _S(parameters().get("S", 0)), _CMax(parameters().get("CMax", 0)) {
	if (_S > 0) {
		logfile().list("Estimating ", _S, " fragment-length specific CT- and GA-PMD patters at 3'-end. (Parameter 'S')");
	} else {
		logfile().list("Estimating a single CT and GA-PMD pattern at 3'-end. (Use Paramter 'S')");
	}

	if (_CMax > 0) {
		logfile().list("Assuming a maximum sequencing cycles of ", _CMax, ". (parameter 'CMax')");
	} else {
		logfile().list("Will assume longest read-length = maximum sequencing cycles.");
	}
}

std::pair<BAM::End, size_t> TPsi::_end_index(const BAM::TSequencedData &data) const noexcept {
	const auto end = data.get<BAM::Flags::Paired>() ? End::from5 : data.end();
	// 5'-end
	if (end == BAM::End::from5) return {end, 0};

	// 3'-end
	if (data.readLength() <= _CMax - _S) return {end, 1};

	const auto tIndex = data.readLength() - (_CMax - _S) + 1;
	if (tIndex >= _tables.size()) return {end ,_tables.size() - 1};

	return {end, tIndex};
}

coretools::Probability TPsi::prob(Type From_To, const BAM::TSequencedData &data) const noexcept {
	using BAM::End;
	const auto realType   = data.get<BAM::Flags::ReversedStrand>() ? _flip[From_To] : From_To;
	const auto [end, idx] = _end_index(data);
	const auto pos        = data.dist(end).pseudo();
	const auto &table     = _tables[idx][realType];

	if (pos >= table.size()) return table.back();
	return table[pos];
}

void TPsi::add(Type From_To, const BAM::TSequencedData &data, coretools::Probability P_g_I_di,
			   const TBaseBaseProbabilities &P_b_bbar_I_gdij) {
	assert(!_tableSums.empty());
	assert(!_tables.empty());
	using BAM::End;
	using genometools::Base;

	const auto realType   = data.get<BAM::Flags::ReversedStrand>() ? _flip[From_To] : From_To;
	const auto [end, idx] = _end_index(data);
	auto &tSum            = _tableSums[idx][realType];
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

	const auto end      = data.end();
	const auto idx      = end == End::from5 ? 0 : data.readLength();
	const auto pos      = data.dist(end).pseudo();
	const auto base     = data.base;
	const auto realType = data.get<BAM::Flags::ReversedStrand>() ? _flip[From_To] : From_To;

	if (_tableSums.size() <= idx) _tableSums.resize(idx + 1);
	auto &tSum = _tableSums[idx][realType];
	if (tSum.size() <= pos) tSum.resize(pos + 1, {{{0, 0, 0, 0}}});

	const auto from = _from[From_To];
	const auto to   = _to[From_To];
	if (ref == from) {
		if (base == to) tSum[pos].fromTo.fromTo += 1.;
		tSum[pos].fromTo.fromSum += 1.;
	} else if (ref == to) {
		if (base == from) tSum[pos].fromTo.toFrom += 1.;
		tSum[pos].fromTo.toSum += 1.;
	}
}

void TPsi::estimate() noexcept {
	constexpr double PMDmin = 1e-9;
	for (size_t i = 0; i < _tables.size(); ++i) {
		//const auto t = type(e);
		for (auto t = Type::min; t < Type::max; ++t) {
			auto &table  = _tables[i][t];
			auto &tSum   = _tableSums[i][t];

			table.clear();
			for (size_t i = 0; i < tSum.size(); ++i) {
				auto &ts       = tSum[i];
				const auto PMD = std::max(PMDmin, ts.numDenom.num / ts.numDenom.denom);
				table.emplace_back(PMD);
				ts.numDenom.num   = 0.;
				ts.numDenom.denom = std::numeric_limits<double>::min(); // preventing any division by 0
			}
			if (table.empty()) table.emplace_back(PMDmin);
		}
	}
}

void TPsi::_initEnd(size_t Idx, int32_t MinData) {
	constexpr size_t minSize = 15; // minFragmentLenght / 2
	constexpr double PMDmin = 1e-9;

	coretools::TStrongArray<double, Type> sums{};
	for (auto t = Type::min; t < Type::max; ++t) {
		auto &table = _tables[Idx][t];
		table.clear();

		auto &tSum = _tableSums[Idx][t];
		if (tSum.empty()) {
			table.emplace_back(PMDmin);
			continue;
		}

		while (tSum.size() > minSize) {
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

void TPsi::_joinTables(size_t From, size_t To) noexcept {
	for (auto t = Type::min; t < Type::max; ++t) {
		auto &to         = _tableSums[To][t];
		const auto &from = _tableSums[From][t];

		if (from.size() > to.size()) to.resize(from.size());

		for (size_t i = 0; i < from.size(); ++i) {
			// if ts5.size() > ts3.size(), we only go to ts3.size()
			to[i].fromTo.fromTo += from[i].fromTo.fromTo;
			to[i].fromTo.fromSum += from[i].fromTo.fromSum;
			to[i].fromTo.toFrom += from[i].fromTo.toFrom;
			to[i].fromTo.toSum += from[i].fromTo.toSum;
		}
	}
}

void TPsi::estimateInit(std::string_view OutputName, size_t MinData) noexcept {
	DEBUG_ASSERT(_nPaired > 0 || _nSingle > 0);

	_printTable(OutputName);

	if (_nPaired > 0 && _nSingle > 0) logfile().warning("Readgroup contains both single- and paired-ended reads!");
	if (paired()) {
		// only 5' end
		logfile().list("Assuming paired-ended reads.");

		for (size_t i = 1; i < _tableSums.size(); ++i) _joinTables(i, 0);
		_tableSums.resize(1);
		_tables.resize(_tableSums.size());
	} else {
		logfile().list("Assuming single-ended reads.");
		if (_CMax == 0) _CMax = _tableSums.size() - 1;
		logfile().list("Setting sequencing cycles to ", _CMax, ".");

		// If CMax was set by user and is smaller than max readlength
		for (size_t i = _tableSums.size() - 1; i > _CMax; --i) {
			_joinTables(i, _CMax);
			_tableSums.pop_back();
		}

		// CMax + 1 -_S may be larger than size()
		const auto maxJoin = std::min(_tableSums.size(), _CMax + 1 - _S);
		for (size_t i = 2; i < maxJoin; ++i) _joinTables(i, 1);

		_tableSums.erase(_tableSums.begin() + 2, _tableSums.begin() + maxJoin);
		_tables.resize(_tableSums.size());

		for (size_t i = 1; i < _tableSums.size(); ++i) _initEnd(i, MinData);
	}
	_initEnd(0, MinData);
}

BAM::RGInfo::TInfo TPsi::info() const {
	BAM::RGInfo::TInfo info;

	// CT5, GA5
	for (auto t = Type::min; t < Type::max; ++t) {
		const auto &v  = _tables[0][t];
		const auto key = impl::toString(t) + "5";
		for (auto p : v) info[key].push_back(p.get()); // explicit conversion from probability to double
	}

	// CT3, GA3
	if (_tables.size() > 1) {
		info["CMax"] = _CMax;
		for (auto t = Type::min; t < Type::max; ++t) {
			const auto key = impl::toString(t) + "3";
			for (size_t i = 1; i < _tables.size(); ++i) {
				const auto name = i == 0 ? "5" : toString("3_", _CMax - _S - 1 + i);
				const auto &v   = _tables[i][t];
				for (auto p : v) info[key][i - 1].push_back(p.get()); // explicit conversion from probability to double
			}
		}
	}
	return info;
}

TPsi::TPsi(const BAM::RGInfo::TInfo &Info) {
	using BAM::RGInfo::toString;
	_tables.resize(2); // at least
	if (Info.is_string()) {
		// CT5:0.1,0.09;GA3:0.3*exp()
		coretools::str::TSplitter semi(Info.get<std::string_view>(), ';');
		for (auto s : semi) {
			const auto sType = readBefore(s, ':');
			user_assert(sType.size() == 3, sType, " is not a recognized token!");

			const auto type = impl::type(sType.substr(0, 2));
			const auto end  = impl::end(sType.substr(2, 1));

			const auto sFunction = readAfter(s, ':');
			if (sFunction.find('*') != sFunction.npos)
				_tables[end][type] = impl::exp(sFunction);
			else
				_tables[end][type] = impl::empiric(sFunction);
		}
	} else {
		_CMax = Info.value("CMax", 0);
		_S    = 0;
		_tables.resize(2); // at least

		// from5
		for (auto t = Type::min; t < Type::max; ++t) {
			const auto k = impl::toString(t) + "5";
			if (Info.contains(k)) {
				if (Info[k].empty()) {
					_tables[0][t] = {P(0.)};
				} else if (Info[k].is_number()) {
					const auto d = Info[k].get<double>();
					user_assert(d >= 0 && d <= 1, "PMD must be between 0 and 1!");
					_tables[0][t] = {P(d)};
				} else if (Info[k].is_array()) {
					for (const auto &d : Info[k]) {
						user_assert(d >= 0 && d <= 1, "PMD must be between 0 and 1!");
						_tables[0][t].emplace_back(d);
					}
				} else {
					throw coretools::TUserError("Cannot parse json-token ", Info[k], "!");
				}
			}
		}

		// from3
		for (auto t = Type::min; t < Type::max; ++t) {
			const auto k = impl::toString(t) + "3";
			if (Info.contains(k)) {
				if (Info[k].empty()) {
					_tables[1][t] = {P(0.)};
				} else if (Info[k].is_number()) {
					const auto d = Info[k].get<double>();
					user_assert(d >= 0 && d <= 1, "PMD must be between 0 and 1!");

					_tables[1][t] = {P(d)};
				} else if (Info[k].is_array()) {
					if (Info[k][0].is_number()) {
						for (const auto &d : Info[k]) {
							user_assert(d >= 0 && d <= 1, "PMD must be between 0 and 1!");
							_tables[1][t].emplace_back(d);
						}
					} else if (Info[k][0].is_array()) {
						if (_CMax == 0 && Info[k].size() > 1)
							throw coretools::TUserError("several ", k + "-PMD on 3'-end needs a Cmax value > 0!");
						if (_S == 0) _S = Info[k].size() - 1;

						if (_tables.size() < 2 + _S) _tables.resize(2 + _S);
						for (size_t i = 1; i < _tables.size(); ++i) {
							for (const auto &d : Info[k][i - 1]) {
								user_assert(d >= 0 && d <= 1, "PMD must be between 0 and 1!");
								_tables[i][t].emplace_back(d);
							}
						}
					} else {
						throw coretools::TUserError("Cannot parse json-token ", Info[k], "!");
					}
				} else {
					throw coretools::TUserError("Cannot parse json-token ", Info[k], "!");
				}
			}
		}
	}
	// at least one value
	for (auto &t : _tables)
		for (auto &v : t)
			if (v.empty()) v = {P(0.)};
}

void TPsi::log() const noexcept {
	constexpr size_t Nmax = 3;

	bool hasAny = false;
	for (size_t i = 0; i < _tables.size(); ++i) {
		const auto name = i == 0
			? "5"
			: _CMax == 0 ? "3" : toString("3_", _CMax - _S - 1 + i);
		for (auto t = Type::min; t < Type::max; ++t) {
			const auto &v = _tables[i][t];
			hasAny        = true;
			auto ret      = impl::toString(t) + name + ": [";
			if (v.size() <= Nmax * 2) {
				for (auto p : v) ret.append(toString(p, ", "));
			} else {
				for (size_t i = 0; i < Nmax; ++i) ret.append(toString(v[i], ", "));
				ret.append("..., ");
				const auto iStart = v.size() - 1 - Nmax;
				for (size_t i = 0; i < Nmax; ++i) ret.append(toString(v[iStart + i], ", "));
			}
			ret.pop_back();
			ret.back() = ']';
			logfile().list(ret);
		}
	}
	if (!hasAny) logfile().list("[]");
}

void TPsi::_printTable(std::string_view OutputName) {
	size_t max = 0;
	for (const auto& outer: _tableSums) {
		for (const auto& inner: outer) {
			max = std::max(inner.size(), max);
		}
	}
	if (max == 0) return;
	std::vector<std::string> header{"End_FragLength"};
	for (size_t i = 0; i < max; ++i) {
		header.push_back(toString(i, "_from5"));
	}

	coretools::TOutputFile oFile(toString(OutputName, "_PsiTable.txt.gz"));
	oFile.writeln("#Format: [fromTo/fromSum:toFrom/toSum:PMDinit]");
	oFile.writeln(header);
	logfile().list("Writing countTable '", oFile.name(), "'.");

	for (size_t i = 0; i < _tableSums.size(); ++i) {
		const auto name = i == 0 ? "5" : toString("3_", i);
		for (auto t = Type::min; t < Type::max; ++t) {
			auto &tSum = _tableSums[i][t];
			if (tSum.empty()) continue;

			oFile.writeNoDelim(impl::toString(t), name).writeDelim();
			for (const auto ts : tSum) {
				const auto fromTo = double(ts.fromTo.fromTo) / ts.fromTo.fromSum;
				const auto toFrom = double(ts.fromTo.toFrom) / ts.fromTo.toSum;
				const auto PMD    = std::max(0., (fromTo - toFrom) / (1.0 - toFrom));
				oFile.writeNoDelim(ts.fromTo.fromTo, "/", ts.fromTo.fromSum, ":", ts.fromTo.toFrom, "/",
								   ts.fromTo.toSum, ":", PMD).writeDelim();
			}
			for (size_t i = tSum.size(); i < max; ++i) {
				oFile.write("0/0:0/0:0");
			}
			oFile.endln();
		}
	}
}

} // namespace GenotypeLikelihoods::PMD
