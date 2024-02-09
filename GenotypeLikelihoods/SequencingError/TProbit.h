#ifndef SEQUENCINGERROR_TPROBIT_H_
#define SEQUENCINGERROR_TPROBIT_H_

#include "SequencingError/TFunction.h"
#include "coretools/Distributions/TNormalDistr.h"

namespace GenotypeLikelihoods::SequencingError {
	
template<typename Covariate> class TProbit final : public TFunction {
private:
	struct TProbitTmpStorage {
		double phi;
		double phiCumul;

		TProbitTmpStorage(const std::array<double, 3> &betas, uint16_t q) {
			const double z = betas[1] + betas[2] * q;
			phiCumul       = coretools::probdist::TNormalDistr::cumulativeDensity(z, 0, 1);
			phi            = coretools::probdist::TNormalDistr::density(z, 0, 1);
		}
	};
	std::array<double, 3> _betas{1., 0, 1.};    // betas of the model

	// tmp storage
	mutable std::vector<TProbitTmpStorage> _tmpStorage;

	void _expandTmpStorage(size_t MaxValue) const {
		for (size_t q = _tmpStorage.size(); q <= MaxValue; ++q) { _tmpStorage.emplace_back(_betas, q); }
	}

public:
	static constexpr std::string_view name = "probit";
	TProbit(size_t FirstParameterIndex) : TFunction(FirstParameterIndex) {}

	size_t numParameters() const noexcept override { return 3; }

	double *begin() noexcept override { return _betas.data(); }
	double *end() noexcept override { return _betas.data() + _betas.size(); }
	const double *begin() const noexcept override { return _betas.data(); }
	const double *end() const noexcept override { return _betas.data() + _betas.size(); }

	void init(const RecalEstimatorTools::TRecalDataTable &, size_t FirstParameterIndex) noexcept override {
		_firstParameterIndex = FirstParameterIndex;
	}

	double getEta(const BAM::TSequencedBase &base) const noexcept override {
		const auto q = Covariate::extract(base);
		if (q >= _tmpStorage.size()) { _expandTmpStorage(q); }
		return _tmpStorage[q].phiCumul*_betas.front();
	}

	double getEta(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
				  std::vector<T2ndDerivative> &der2) const noexcept override {
		const auto q = Covariate::extract(base);
		if (q >= _tmpStorage.size()) { _expandTmpStorage(q); }

		const auto b1 = firstParameterIndex();
		const auto b2 = b1 + 1;
		const auto b3 = b1 + 2;

		const auto z        = _betas[1] + _betas[2] * q;
		const auto phi      = _tmpStorage[q].phi;
		const auto phiCumul = _tmpStorage[q].phiCumul;
		const double bPhi   = phi * _betas.front();
		const double bPhiQ  = bPhi * q;
		const double bPhiQZ = bPhiQ * z;

		der1.emplace_back(b1, phiCumul);
		der1.emplace_back(b2, bPhi);
		der1.emplace_back(b3, bPhiQ);

		//der2.emplace_back(b1, b1, 0.); this is zero
		der2.emplace_back(b1, b2, phi);
		der2.emplace_back(b1, b3, phi*q);

		//der2.emplace_back(b2, b1, phi); only upper half
		der2.emplace_back(b2, b2, -bPhi*z);
		der2.emplace_back(b2, b3, -bPhiQZ*q);

		//der2.emplace_back(b3, b1, phi*p); only upper half
		//der2.emplace_back(b3, b2, -bPhi*z); only upper half
		der2.emplace_back(b3, b3, -bPhiQZ);

		return phiCumul *_betas.front();
	}

	double adjust() noexcept override { return 0.; }
	std::string typeString() const noexcept override { return std::string(Covariate::name).append(1, ':').append(name); }

	void addInfo(BAM::RGInfo::TInfo& info) const override {
		info[Covariate::name] = {{name, _betas}};
	}
};
}
#endif
