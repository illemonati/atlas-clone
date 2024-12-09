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

std::string toString(End end) {
	if (end == End::from5) return "5";
	else return "3";
}



} // namespace impl

void TPsi::_fromString(std::string_view Psi) {
	// CT5:0.1,0.09;GA3:0.3*exp()
	coretools::str::TSplitter semi(Psi, ';');
	for (auto s: semi) {
		const auto sType     = readBefore(s, ':');
		if (sType.size() != 3) UERROR(sType, " is not a recognized token!");

		const auto type      = impl::type(sType.substr(0, 2));
		const auto end       = impl::end(sType.substr(2, 1));
		const auto sFunction = readAfter(s, ':');
		if (sFunction.find('*') != sFunction.npos) _tables[end][type] = impl::exp(sFunction);
		else _tables[end][type] = impl::empiric(sFunction);
	}
	for (auto& t: _tables)
		for (auto &v: t)
			if (v.empty()) v = {P(0.)};
}

TPsi::TPsi(const BAM::RGInfo::TInfo &Info) {
	_parse(Info);
}

TPsi::TPsi() {
	for (auto &t : _tables)
		for (auto &v : t) v.emplace_back(0.);
}

void TPsi::_parse(const BAM::RGInfo::TInfo & Info) {
	if (Info.is_string()) {
		_fromString(Info.get<std::string_view>());
	} else {
		for (auto e = End::min; e < End::max; ++e) {
			for (auto t = Type::min; t < Type::max; ++t) {
				auto &v        = _tables[e][t];
				const auto key = impl::toString(t) + impl::toString(e);
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
			const auto &v = _tables[e][t];
			if (v.size() > 1 || v.front() > 0.) {
				const auto key = impl::toString(t) + impl::toString(e);
				for (auto p : v) info[key].push_back(p.get()); // explicit conversion from probability to double
			}
		}
	}
	return info;
}

void TPsi::estimate() noexcept {
	constexpr double PMDmin = 1e-9;
	for (auto e = End::min; e < End::max; ++e) {
		const auto t = type(e);
		auto &table  = _tables[e][t];
		auto &tSum   = _tableSums[e][t];

		table.clear();
		//for (auto &ts : tSum) {
		for (size_t i = 0; i < tSum.size(); ++i) {
			auto& ts = tSum[i];
			const auto PMD = std::max(PMDmin, ts.numDenom.num / ts.numDenom.denom);
			table.emplace_back(PMD);
			ts.numDenom.num   = 0.;
			ts.numDenom.denom = std::numeric_limits<double>::min(); // preventing any division by 0
		}
	}
}

void TPsi::_printTable(std::string_view FName) {
	coretools::TOutputFile oFile(FName);
	for (auto e = End::min; e < End::max; ++e) {
		for (auto t = Type::min; t < Type::max; ++t) {
			oFile.write(impl::toString(t) + impl::toString(e));
			auto &tSum = _tableSums[e][t];
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

void TPsi::_initEnd(End e) {
	constexpr int Nmin      = 100;
	constexpr double PMDmin = 1e-9;

	coretools::TStrongArray<double, Type> sums{};
	for (auto t = Type::min; t < Type::max; ++t) {
		auto &table = _tables[e][t];
		table.clear();

		auto &tSum = _tableSums[e][t];
		if (tSum.empty()) continue;

		while (tSum.size() > 1) {
			const auto &ts = tSum.back();
			auto &ts_m     = *(tSum.end() - 2);

			if (ts.fromTo.fromSum < Nmin || ts.fromTo.toSum < Nmin) {
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
	// Either CT or GA
	const auto worseType = sums[Type::CT] >= sums[Type::GA] ? Type::GA : Type::CT;
	_tableSums[e][worseType].clear();
	_tables[e][worseType] = {P(0.)};
}

void TPsi::_joinTables() noexcept {
	for (auto t = Type::min; t < Type::max; ++t) {
		const auto &ts3 = _tableSums[End::from3][t];
		auto &ts5       = _tableSums[End::from5][t];

		if (ts3.size() > ts5.size())
			ts5.resize(ts3.size());

		for (size_t i = 0; i < ts3.size(); ++i) {
			// if ts5.size() > ts3.size(), we only go to ts3.size()
			ts5[i].fromTo.fromTo += ts3[i].fromTo.fromTo;
			ts5[i].fromTo.fromSum += ts3[i].fromTo.fromSum;
			ts5[i].fromTo.toFrom += ts3[i].fromTo.toFrom;
			ts5[i].fromTo.toSum += ts3[i].fromTo.toSum;
		}
	}
}

void TPsi::estimateInit(std::string_view OutputName) noexcept {
	using coretools::instances::logfile;
	const auto fn = toString(OutputName, "_PsiTable.txt.gz");
	logfile().list("Writing countTable '", fn, "'.");
	_printTable(toString(OutputName, "_PsiTable.txt.gz"));

	if (paired()) {
		logfile().list("Assuming paired-ended read.");

		// only 5' end
		_joinTables();
		_tables[End::from3][Type::CT] = {P(0.)};
		_tables[End::from3][Type::GA] = {P(0.)};
	} else {
		logfile().list("Assuming single-ended read.");
		_initEnd(End::from3);
	}
	_initEnd(End::from5);
}

void TPsi::log() const noexcept {
	using coretools::instances::logfile;
	constexpr size_t Nmax = 3;

	bool hasAny = false;
	for (auto e = End::min; e < End::max; ++e) {
		if (paired() && e == End::from3) continue;

		for (auto t = Type::min; t < Type::max; ++t) {
			const auto &v = _tables[e][t];
			if (v.size() > 1 || v.front() > 0.) {
				hasAny = true;
				auto ret = impl::toString(t) + impl::toString(e) + ": [";
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

void TPsi::reset(const BAM::RGInfo::TInfo &info) {
	for (auto e = End::min; e < End::max; ++e) {
		for (auto t = Type::min; t < Type::max; ++t) {
			_tables[e][t].clear();
		}
	}
	_parse(info);

	const auto stInit = []() {
		SumType stInit{};
		stInit.numDenom.num   = 0.;
		stInit.numDenom.denom = std::numeric_limits<double>::min(); // preventing any division by 0
		return stInit;
	}();

	for (auto e = End::min; e < End::max; ++e) {
		if ((_tables[e][Type::CT].size() > _tables[e][Type::GA].size())) {
			_tableSums[e][Type::CT].assign(_tables[e][Type::CT].size(), stInit);
		} else {
			_tableSums[e][Type::GA].assign(_tables[e][Type::GA].size(), stInit);
		}
	}
}

} // namespace GenotypeLikelihoods::PMD
