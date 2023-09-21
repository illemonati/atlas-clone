#include "TRho.h"
#include "TGenotypeData.h"

namespace GenotypeLikelihoods::SequencingError {
using namespace coretools::str;
using genometools::Base;
using namespace std::literals;

TRho::TRho(std::string_view Def) {
	using coretools::str::toString;
	//"default" implies default rho

	if (Def.empty() || Def == "default" || Def == "-") {
		return;
	}

	TSplitter spl(Def, ';');
	size_t i = 0;
	for (auto s: spl) {
		if (i >= 4) UERROR("Too many rows given for rho, needed only 4!");

		std::array<double, 4> ar;
		TSplitter spl2(strip(s, "[]"), ',');
		size_t j = 0;

		for (auto ss : spl2) {
			if (j >= ar.size()) UERROR("Too many rho values given for row ", i, ", needed ", ar.size(), "!");

			if (strip(ss) == "-") {
				ar[j] = 0.;
			} else {
				coretools::str::fromString<true>(strip(ss), ar[j]);
			}
			++j;
		}
		if (j < ar.size()) UERROR("Too few(", j, ") rho values given, needed ", ar.size(), "!");

		ar[i]         = 0.;
		_rho[Base(i)] = TBaseProbabilities::normalize(ar);
		++i;
	}
	if (i < 4) UERROR("Too few rows given for rho, needed 4, not ", i, "!");
}

TRho::TRho(const BAM::RGInfo::TInfo &info) {
	auto b = Base::min;
	for (const std::array<double, 4> line: info) {
		_rho[b] = TBaseProbabilities::normalize(line);
		++b;
	}
}

void TRho::log() const {
	using coretools::str::toString;
	using coretools::instances::logfile;
	logfile().list("[[   -    , ", _rho[Base::A][Base::C], ", ", _rho[Base::A][Base::G], ", ", _rho[Base::A][Base::T], "]");
	logfile().list(" [", _rho[Base::C][Base::A], ",    -    , ", _rho[Base::C][Base::G], ", ", _rho[Base::C][Base::T], "]");
	logfile().list(" [", _rho[Base::G][Base::A], ", ", _rho[Base::G][Base::C], ",    -    , ", _rho[Base::G][Base::T], "]");
	logfile().list(" [", _rho[Base::T][Base::A], ", ", _rho[Base::T][Base::C], ", ", _rho[Base::T][Base::G], ",    -    ]]");
}

BAM::RGInfo::TInfo TRho::info() const {
	return {
		{_rho[Base::A][Base::A].get(), _rho[Base::A][Base::C].get(), _rho[Base::A][Base::G].get(), _rho[Base::A][Base::T].get()},
		{_rho[Base::C][Base::A].get(), _rho[Base::C][Base::C].get(), _rho[Base::C][Base::G].get(), _rho[Base::C][Base::T].get()},
		{_rho[Base::G][Base::A].get(), _rho[Base::G][Base::C].get(), _rho[Base::G][Base::G].get(), _rho[Base::G][Base::T].get()},
		{_rho[Base::T][Base::A].get(), _rho[Base::T][Base::C].get(), _rho[Base::T][Base::G].get(), _rho[Base::T][Base::T].get()}
	};
}

void TRho::add(genometools::Base l, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept {
	for (auto k = Base::min; k < Base::max; ++k) {
		_rhoSum[k][l] += P_g_I_d*P_bbar_I_d[k];
	}
}

void TRho::estimate() noexcept {
	for (Base k = Base::min; k < Base::max; ++k) {
		_rhoSum[k][k] = 0.0;
		_rho[k] = TBaseProbabilities::normalize(_rhoSum[k]);
	}
	// reset
	_rhoSum.fill({});
}
}

