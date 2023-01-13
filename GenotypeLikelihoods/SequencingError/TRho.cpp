#include "TRho.h"

namespace GenotypeLikelihoods::SequencingError {
using namespace coretools::str;
using genometools::Base;
using namespace std::literals;

TRho::TRho(std::string_view Def) {
	using coretools::str::toString;
	//"default" implies default rho

	if (Def == "default" || Def == "-") {
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

TRho::TRho(const BAM::RGInfo::TInfo &info) { ECHO(info.dump()); }

std::string TRho::definition() const noexcept {
	using coretools::str::toString;
	return "[[-,"s + toString(_rho[Base::A][Base::C]) + ',' + toString(_rho[Base::A][Base::G]) + ',' + toString(_rho[Base::A][Base::T]) + "];["
		+ toString(_rho[Base::C][Base::A]) + ",-," + toString(_rho[Base::C][Base::G]) + ',' + toString(_rho[Base::C][Base::T]) + "];["
		+ toString(_rho[Base::G][Base::A]) + ',' + toString(_rho[Base::G][Base::C]) + ",-," + toString(_rho[Base::G][Base::T]) + "];["
		+ toString(_rho[Base::T][Base::A]) + ',' + toString(_rho[Base::T][Base::C]) + ',' + toString(_rho[Base::T][Base::G]) + ",-]]";
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

