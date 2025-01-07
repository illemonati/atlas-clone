#ifndef SEQUENCINGERROR_TEMPIRIC_H_
#define SEQUENCINGERROR_TEMPIRIC_H_

#include "SequencingError/TCovariate.h"
#include "SequencingError/TFunction.h"
#include "SequencingError/TRecalDataTable.h"
#include "coretools/Main/TLog.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Types/TPseudoInt.h"

namespace GenotypeLikelihoods::SequencingError {

template<typename Covariate> class TEmpiric final : public TFunction {
private:
	std::vector<double> _betas;    // betas of the model

public:
	static constexpr std::string_view name = "empiric";

	TEmpiric(size_t FirstParameterIndex) : TFunction(FirstParameterIndex) {}

	size_t numParameters() const noexcept override { return _betas.size(); }

	double *begin() noexcept override { return _betas.data(); }
	double *end() noexcept override { return _betas.data() + _betas.size(); }
	const double *begin() const noexcept override { return _betas.data(); }
	const double *end() const noexcept override { return _betas.data() + _betas.size(); }

	void push_back(double val) noexcept {_betas.push_back(val);}

	void init(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t FirstParameterIndex) override {
		_firstParameterIndex = FirstParameterIndex;
		_betas.assign(dataTable[Covariate::index].size(), NAN);
		for (size_t i = 0; i < _betas.size(); ++i) {
			if (dataTable[Covariate::index][i]) {
				if constexpr (std::is_same_v<Covariate, TCovariate_quality>) {
					const coretools::Probability p = coretools::Probability(coretools::PhredInt(i));
					_betas[i] = coretools::logit(p);
				} else {
					_betas[i] = 0.;
				}
			}
		}
	}

	double adjust() noexcept override {
		double mean = 0.;
		size_t N    = 0;
		for (auto bi : _betas) {
			if (!std::isnan(bi)) {
				++N;
				mean += bi;
			}
		}
		if (N > 1) mean /= N;

		for (auto &bi : _betas) {
			if (!std::isnan(bi)) { bi -= mean; }
		}
		return mean;
	}

	double getEta(const BAM::TSequencedBase &base) const noexcept override {
		const auto val = Covariate::extract(base);
		if (val < _betas.size()) return _betas[val];

		return _betas.back();
	}

	double getEta(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
				  std::vector<T2ndDerivative> &) const noexcept override {
		const auto val = std::min<size_t>(Covariate::extract(base), _betas.size() - 1);

		const size_t der_index = firstParameterIndex() + static_cast<size_t>(val);
		der1.emplace_back(der_index, 1.0);
		return _betas[val];
	}

	std::string typeString() const noexcept override {
		return std::string(Covariate::name).append(1, ':').append(name);
	}

	void addInfo(BAM::RGInfo::TInfo &info) const override {
		BAM::RGInfo::TInfo ar = nlohmann::json::array();
		for (size_t i = 0; i < _betas.size(); ++i) {
			if (!std::isnan(_betas[i])) {
				if constexpr (std::is_same_v<Covariate, TCovariate_position>) {
					ar += {coretools::TPseudoInt::fromPseudo(i).linear(), _betas[i]};
				} else if constexpr (std::is_same_v<Covariate, TCovariate_fragmentLength>) {
					ar += {coretools::TLogInt::fromLog(i).linear(), _betas[i]};
				} else {
					ar += {i, _betas[i]};
				}
			}
		}
		info[Covariate::name] = {{name, ar}};
	}

	void log() const override {
		using coretools::str::toString;
		using coretools::instances::logfile;
		constexpr size_t Nmax = 3;

		std::vector<size_t> iis;
		iis.reserve(_betas.size());
		for (size_t i = 0; i < _betas.size(); ++i) {
			if (!std::isnan(_betas[i])) iis.push_back(i);
		}
		if (iis.empty()) {
			logfile().list(typeString(), ": []");
			return;
		}

		std::string ret = "[";
		if (iis.size() <= 2 * Nmax) {
			for (auto i : iis) ret.append(toString(i, ": ", _betas[i], ", "));
		} else {
			for (size_t j = 0; j < Nmax; ++j)
				ret.append(toString(iis[j], ": ", _betas[iis[j]], ", "));
			ret.append("..., ");
			const auto jStart = iis.size() - Nmax;
			for (size_t j = 0; j < Nmax; ++j)
				ret.append(toString(iis[jStart + j], ": ", _betas[iis[jStart + j]], ", "));
		}
		ret.pop_back();
		ret.back() = ']';
		logfile().list(typeString(), ": ", ret);
	}
};
}

#endif
