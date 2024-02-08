#include "TPsi.h"

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Strings/stringManipulations.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"
#include <cmath>
#include <limits>

namespace GenotypeLikelihoods::PMD {
using coretools::Probability;
using namespace coretools::str;

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
	TSplitter spl(sFunction, ',');
	if (spl.empty()) return {0.};

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
	TSplitter semi(Psi, ';');
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
			if (v.empty()) v = {0.};
}

TPsi::TPsi(std::string_view Psi) {
	_fromString(Psi);
}

TPsi::TPsi(const BAM::RGInfo::TInfo &info) {
	if (info.is_string()) {
		_fromString(info.get<std::string_view>());
	} else {
		for (auto e = End::min; e < End::max; ++e) {
			for (auto t = Type::min; t < Type::max; ++t) {
				auto &v        = _tables[e][t];
				const auto key = impl::toString(t) + impl::toString(e);
				if (!info.contains(key) || info[key].empty()) {
					v = {0.};
				} else {
					if (info[key].is_string()) {
						const auto sFunction = info[key].get<std::string_view>();
						if (sFunction.find('*') != sFunction.npos)
							v = impl::exp(sFunction);
						else
							v = impl::empiric(sFunction);
					} else if (info[key].is_array()) {
						for (auto [k, d] : info[key].items()) {
							v.emplace_back(d);
						}
					} else {
						UERROR("Cannot parse json-token ", info, "!");
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
	for (auto e = End::min; e < End::max; ++e) {
		const auto t = _tableSums[e][Type::CT].empty() ? Type::GA : Type::CT;
		auto &table  = _tables[e][t];
		auto &tSum   = _tableSums[e][t];

		table.clear();
		for (auto &ts : tSum) {
			table.push_back(ts.numDenom.num / ts.numDenom.denom);
			ts.numDenom.num   = 0.;
			ts.numDenom.denom = std::numeric_limits<double>::min(); // preventing any division by 0
		}
	}
}

void TPsi::estimateInit() noexcept {
	constexpr int Nmin = 100;
	constexpr double psiMin = 1./Nmin/10;
	for (auto e = End::min; e < End::max; ++e) {
		coretools::TStrongArray<size_t, Type> sums{};
		for (auto t = Type::min; t < Type::max; ++t) {
			auto &table = _tables[e][t];
			auto &tSum  = _tableSums[e][t];

			// add up end to have enough data
			for (size_t i = tSum.size() - 1; i > 0; --i) {
				const auto &ts      = tSum[i];
				auto& ts_m = tSum[i - 1]; // i > 0 -> always ok

				auto merge = [](auto ts) {
					return (ts.fromTo.fromTo <= ts.fromTo.toFrom) || (ts.fromTo.fromSum + ts.fromTo.toSum < Nmin);
				};

				if (merge(ts) || merge(ts_m)) {
					ts_m.fromTo.fromSum += ts.fromTo.fromSum;
					ts_m.fromTo.fromTo += ts.fromTo.fromTo;
					ts_m.fromTo.toSum += ts.fromTo.toSum;
					ts_m.fromTo.toFrom += ts.fromTo.toFrom;
					tSum.pop_back();
				} else {
					break;
				}
			}
			table.clear();
			for (auto &ts : tSum) {
				const auto fromTo = double(ts.fromTo.fromTo) / ts.fromTo.fromSum;
				const auto toFrom = double(ts.fromTo.toFrom) / ts.fromTo.toSum;
				table.push_back(std::max(psiMin, (fromTo - toFrom) / (1.0 - toFrom))); // once 0, always 0, so add something small
				sums[t] += std::max(0, ts.fromTo.fromTo - ts.fromTo.toFrom);
				ts.numDenom.num   = 0.;
				ts.numDenom.denom = std::numeric_limits<double>::min(); // preventing any division by 0
			}
		}
		// Either CT or GA
		const auto worseType = sums[Type::CT] >= sums[Type::GA] ? Type::GA : Type::CT;
		_tableSums[e][worseType].clear();
		_tables[e][worseType] = {0.};
	}
}

void TPsi::log() const noexcept {
	using coretools::instances::logfile;
	constexpr size_t Nmax = 3;

	bool hasAny = false;
	for (auto e = End::min; e < End::max; ++e) {
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

} // namespace GenotypeLikelihoods::PMD
