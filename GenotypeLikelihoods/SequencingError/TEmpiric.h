#ifndef SEQUENCINGERROR_TEMPIRIC_H_
#define SEQUENCINGERROR_TEMPIRIC_H_

#include "SequencingError/TCovariate.h"
#include "SequencingError/TFunction.h"
#include "SequencingError/TRecalDataTable.h"
#include "coretools/Main/TLog.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Types/TPseudoInt.h"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>

namespace GenotypeLikelihoods::SequencingError {

template<typename Covariate> class TEmpiric final : public TFunction {
private:
	static constexpr size_t _N = [](){
		if constexpr (std::is_same_v<Covariate, TCovariate_context>) {
			return 5;
		} else if constexpr (std::is_same_v<Covariate, TCovariate_fragmentLength>) {
			return 32; //2**32 > length of largest chromosome
		} else {
			return 256;
		}
	}();
	static constexpr uint8_t _nope = -1;

	std::vector<double> _vals;       // betas of the model
	std::array<uint8_t, _N> _iis{};  // index to betas or _nope

	double _beta(size_t i) const noexcept {return _vals[_iis[i]];}
	double* _beta(size_t i) noexcept {return &_vals[_iis[i]];}

public:
	static constexpr std::string_view name = "empiric";

	TEmpiric(size_t FirstParameterIndex) : TFunction(FirstParameterIndex) {
		_iis.fill(_nope);
	}

	size_t numParameters() const noexcept override { return _vals.size(); }

	double *begin() noexcept override { return _vals.data(); }
	double *end() noexcept override { return _vals.data() + _vals.size(); }
	const double *begin() const noexcept override { return _vals.data(); }
	const double *end() const noexcept override { return _vals.data() + _vals.size(); }

	void setData(std::vector<std::pair<size_t, double>> Data) {
		// make sure
		std::sort(Data.begin(), Data.end(), [](const auto& p1, const auto& p2) {
			return p1.first < p2.first;
		});
		// everything before 1st value = 1st value
		size_t i = 0;
		for (const auto &p : Data) {
			_vals.push_back(p.second);
			// unlearned indexes get pooled with next higher up
			for (;i <= p.first; ++i) {
				_iis[i] = _vals.size() - 1;
			}
		}
		// everithing after last value = last value
		for (; i < _iis.size(); ++i) {
			_iis[i] = _vals.size() - 1;
		}
	}

	void init(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t FirstParameterIndex, size_t MinData) override {
		_vals.clear();
		_iis.fill(_nope);

		_firstParameterIndex = FirstParameterIndex;

		if constexpr (std::is_same_v<Covariate, TCovariate_quality>) {
			// don't pool qualities
			const auto& table = dataTable[TCovariate_quality::index];
			for (size_t i = 0; i < table.size(); ++i) {
				if (table[i]) {
					const coretools::Probability p = coretools::Probability(coretools::PhredInt(i));
					_vals.push_back(coretools::logit(p));
					_iis[i] = _vals.size() - 1;
				}
			}
		} else if constexpr (std::is_same_v<Covariate, TCovariate_context>) {
			// don't pool context
			const auto& table = dataTable[TCovariate_context::index];
			for (size_t i = 0; i < table.size(); ++i) {
				if (table[i]) {
					_vals.push_back(0.);
					_iis[i] = _vals.size() - 1;
				}
			}
		} else {
			const auto& table = dataTable[Covariate::index];
			std::vector<std::vector<size_t>> pool;
			for (size_t i = 0; i < table.size(); ++i) {
				if (table[i]) {
					_vals.push_back(table[i]);
					pool.push_back({i});
				}
			}
			assert(_vals.size() == pool.size());

			// Lower
			while (_vals.size() > 2 && _vals.front() < MinData) {
				_vals[1] += _vals[0];
				pool[1].insert(pool[1].end(), pool[0].begin(), pool[0].end());
				_vals.erase(_vals.begin());
				pool.erase(pool.begin());
			}
			// Upper
			while (_vals.size() > 2 && _vals.back() < MinData) {
				const auto i = _vals.size() - 1;
				_vals[i - 1] += _vals[i];
				pool[i - 1].insert(pool[i - 1].end(), pool[i].begin(), pool[i].end());
				_vals.pop_back();
				pool.pop_back();
			}

			// Middle
			size_t iMin = std::distance(_vals.begin(), std::min_element(_vals.begin(), _vals.end()));
			while (_vals[iMin] < MinData) {
				size_t dir = _vals[iMin - 1] < _vals[iMin + 1] ? -1 : 1;
				_vals[iMin + dir] += _vals[iMin];
				pool[iMin + dir].insert(pool[iMin + dir].end(), pool[iMin].begin(), pool[iMin].end());
				_vals.erase(_vals.begin() + iMin);
				pool.erase(pool.begin() + iMin);
				iMin = std::distance(_vals.begin(), std::min_element(_vals.begin(), _vals.end()));
			}

			for (size_t i = 0; i < pool.size(); ++i) {
				assert(!pool[i].empty());
				for (auto j: pool[i]) {
					_iis[j] = i;
				}
				_vals[i] = 0;
			}
		}
		assert(_vals.size() <= size_t(_nope));
	}

	double adjust() noexcept override {
		double mean = 0.;
		for (const auto bi: _vals) { 
			mean += bi;
		}
		if (mean != 0.) mean /= numParameters();

		for (auto& bi: _vals) { 
			bi -= mean;
		}

		return mean;
	}

	double getEta(const BAM::TSequencedBase &base) const noexcept override {
		const auto val = Covariate::extract(base);
		return _beta(val);
	}

	double getEta(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
				  std::vector<T2ndDerivative> &) const noexcept override {
		const auto val = Covariate::extract(base);

		const size_t der_index = firstParameterIndex() + _iis[val];
		der1.emplace_back(der_index, 1.0);
		return _beta(val);
	}

	std::string typeString() const noexcept override {
		return std::string(Covariate::name).append(1, ':').append(name);
	}

	void addInfo(BAM::RGInfo::TInfo &info) const override {
		BAM::RGInfo::TInfo ar = nlohmann::json::array();
		for (size_t i = 0; i < _iis.size(); ++i) {
			if (_iis[i] != _nope) {
				if constexpr (std::is_same_v<Covariate, TCovariate_position>) {
					ar += {coretools::TPseudoInt::fromPseudo(i).linear(), _beta(i)};
				} else if constexpr (std::is_same_v<Covariate, TCovariate_fragmentLength>) {
					ar += {coretools::TLogInt::fromLog(i).linear(), _beta(i)};
				} else {
					ar += {i, _beta(i)};
				}
			}
		}
		info[Covariate::name] = {{name, ar}};
	}

	void log() const override {
		using coretools::str::toString;
		using coretools::instances::logfile;
		constexpr size_t Nmax = 3;

		std::string ret = "[";
		if (numParameters() <= 2 * Nmax) {
			// write all parameters
			for (size_t i = 0; i < _iis.size(); ++i) {
				if (_iis[i] != _nope) {
					ret.append(toString(i, ": ", _beta(i), ", "));
				}
			}
		} else {
			// write first Nmax parameters
			for (size_t i = 0, j=0; j < Nmax; ++i) {
				if (_iis[i] != _nope) {
					ret.append(toString(i, ": ", _beta(i), ", "));
					++j;
				}
			}
			ret.append("..., ");
			// find last Nmax parameters
			std::array<size_t, Nmax> ilast;
			size_t i = _iis.size() - 1;
			for (size_t j = 0; j < Nmax; --i) {
				if (_iis[i] != _nope) {
					++j;
					ilast[Nmax - j] = i;
				}
			}
			for (const auto& i: ilast) {
				ret.append(toString(i, ": ", _beta(i), ", "));
			}
		}
		ret.pop_back();
		ret.back() = ']';
		logfile().list(typeString(), ": ", ret);
	}
};
}

#endif
