#include "TSimulatorMutationtable.h"
#include <numeric>

namespace Simulations {

using genometools::Base;

TSimulatorMutationtable::TSimulatorMutationtable(const genometools::TBaseProbabilities &baseFreq) {
	for (Base a = Base::min; a < Base::max; ++a) {
		for (Base b = Base::min; b < Base::max; ++b) {
			_mutTable[a][b] = baseFreq[a] * baseFreq[b];
		}
		_mutTable[a][a] = 0.0;
	}
	_normalizeAndMakeCumulative();
}

TSimulatorMutationtable::TSimulatorMutationtable(const genometools::TBaseProbabilities &baseFreq,
						 double theta) {
	const double exp_t   = exp(-theta);
	const double o_exp_t = 1 - exp_t;
	for (Base a = Base::min; a < Base::max; ++a) {
		for (Base b = Base::min; b < Base::max; ++b) {
			_mutTable[a][b] = baseFreq[a] * baseFreq[b] * o_exp_t;
		}
		_mutTable[a][a] += baseFreq[a]*exp_t;
	}
	_normalizeAndMakeCumulative();
}

void TSimulatorMutationtable::_normalizeAndMakeCumulative() {
	// normalize within row
	for (auto & mu : _mutTable) {
		const double sum = std::accumulate(mu.begin(), mu.end(), 0.);
		mu[Base::A] /= sum;
		mu[Base::C] = mu[Base::C] / sum + mu[Base::A];
		mu[Base::G] = mu[Base::G] / sum + mu[Base::C];
		mu[Base::T] = 1.0;
	}
}
	
}
