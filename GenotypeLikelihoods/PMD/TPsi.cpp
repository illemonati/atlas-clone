#include "TPsi.h"
#include "coretools/Strings/fromString.h"
#include "coretools/Strings/splitters.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/probability.h"
#include "coretools/devtools.h"
#include "genometools/GenotypeTypes.h"

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
		if (sFunction.find('*') != sFunction.npos) _tables[type][end] = impl::exp(sFunction);
		else _tables[type][end] = impl::empiric(sFunction);
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
		for (auto t = Type::min; t < Type::max; ++t) {
			for (auto e = End::min; e < End::max; ++e) {
				auto &v        = _tables[t][e];
				const auto key = impl::toString(t) + impl::toString(e);
				if (info[key].empty()) {
					v = {0.};
				} else {
					if (info[key].is_string()) {
						const auto sFunction = info[key].get<std::string_view>();
						if (sFunction.find('*') != sFunction.npos)
							v = impl::exp(sFunction);
						else
							v = impl::empiric(sFunction);
					} else if (info[key].is_array()) {
						for (auto d : info[key]) { v.emplace_back(d); }
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
	for (auto t = Type::min; t < Type::max; ++t) {
		for (auto e = End::min; e < End::max; ++e) {
			const auto & v = _tables[t][e];
			if (v.size() > 1 || v.front() > 0.) {
				const auto key = impl::toString(t) + impl::toString(e);
				for (auto p : v) info[key].push_back(p.get()); // explicit conversion from probability to double
			}
		}
	}
	return info;
}

void TPsi::estimate() noexcept {
	for (auto t = Type::min; t < Type::max; ++t) {
		for (auto e = End::min; e < End::max; ++e) {
			auto &table = _tables[t][e];
			auto &tSum  = _tableSums[t][e];

			table.clear();
			for (const auto& ts: tSum) {
				table.push_back(ts.num/ts.denom);
			}
			tSum.clear();
		}
	}
	OUT(_tables);
}

} // namespace GenotypeLikelihoods::PMD
