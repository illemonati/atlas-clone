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
#include <sys/types.h>

namespace GenotypeLikelihoods::SequencingError {

template<typename Covariate> class TEmpiric final : public TFunction {
private:
	static constexpr size_t _N = [](){
		if constexpr (std::is_same_v<Covariate, TCovariate_fragmentLength>) {
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
	double *end() noexcept override { return _vals.data() + numParameters(); }
	const double *begin() const noexcept override { return _vals.data(); }
	const double *end() const noexcept override { return _vals.data() + numParameters(); }

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

		const auto &table = dataTable[Covariate::index];
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
			if constexpr (std::is_same_v<Covariate, TCovariate_quality>) {
				size_t iTot = 0;
				for (auto j : pool[i]) {
					_iis[j] = i;
					iTot += j;
				}
				iTot /= pool[i].size(); // take the mean as initial value
				const auto p = coretools::Probability(coretools::PhredInt(iTot));
				_vals[i]     = coretools::logit(p);
			} else {
				for (auto j : pool[i]) { _iis[j] = i; }
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

	double getEta(const BAM::TSequencedData &data) const noexcept override {
		const auto val = Covariate::extract(data);
		return _beta(val);
	}

	double getEta(const BAM::TSequencedData &data, std::vector<T1stDerivative> &der1,
				  std::vector<T2ndDerivative> &) const noexcept override {
		const auto val = Covariate::extract(data);

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

template<>
class TEmpiric<TCovariate_context> final : public TFunction {
private:
	static constexpr size_t _N = 5;
	std::array<double, _N> _betas;

public:
	static constexpr std::string_view name = "empiric";

	TEmpiric(size_t FirstParameterIndex) : TFunction(FirstParameterIndex) {
		_betas.fill(0.);
	}

	size_t numParameters() const noexcept override { return _N - 1; } // don't count Base::N

	double *begin() noexcept override { return _betas.data(); }
	double *end() noexcept override { return _betas.data() + numParameters(); }
	const double *begin() const noexcept override { return _betas.data(); }
	const double *end() const noexcept override { return _betas.data() + numParameters(); }

	void setData(std::vector<std::pair<size_t, double>> Data) {
		using coretools::instances::logfile;

		if (Data.size() > numParameters() + 1) {
			UERROR("A Can only have  4 contexts: A, C, G, T, but your input has ", Data.size(), " contexts!");
		} else  if (Data.size() == numParameters() + 1) { // Base::N was measured
			logfile().warning("A value for context = N was estimated, this is outdated!");
		} else if  (Data.size() == numParameters() + 1) { // too little data
			logfile().warning("Not all contexts: A, C, G, T were estimated!");
		}
		_betas.fill(0.);
		for (const auto &p : Data) {
			_betas[p.first] = p.second;
		}
	}

	void init(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t FirstParameterIndex, size_t ) override {
		using coretools::instances::logfile;
		_betas.fill(0.);
		_firstParameterIndex = FirstParameterIndex;
		
		const auto& table = dataTable[TCovariate_context::index];
		if (table.size() > numParameters() + 1) { // data table counts Base::N
			UERROR("B Can only have  4 contexts: A, C, G, T, but your data has ", table.size(), " contexts!");
		}
	}

	double adjust() noexcept override {
		double mean = 0.;
		for (size_t i = 0; i < numParameters(); ++i) {
			mean += _betas[i];
		}
		if (mean != 0.) mean /= numParameters();

		for (size_t i = 0; i < numParameters(); ++i) {
			_betas[i] -= mean;
		}

		return mean;
	}

	double getEta(const BAM::TSequencedData &data) const noexcept override {
		const auto val = TCovariate_context::extract(data);
		return _betas[val];
	}

	double getEta(const BAM::TSequencedData &data, std::vector<T1stDerivative> &der1,
				  std::vector<T2ndDerivative> &) const noexcept override {
		const auto val = TCovariate_context::extract(data);

		if (val < numParameters()) {
			const size_t der_index = firstParameterIndex() + val;
			der1.emplace_back(der_index, 1.0);
		}
		return _betas[val];
	}

	std::string typeString() const noexcept override {
		return std::string(TCovariate_context::name).append(1, ':').append(name);
	}

	void addInfo(BAM::RGInfo::TInfo &info) const override {
		BAM::RGInfo::TInfo ar = nlohmann::json::array();
		for (size_t i = 0; i < numParameters(); ++i) {
			ar += {i, _betas[i]};
		}
		info[TCovariate_context::name] = {{name, ar}};
	}

	void log() const override {
		using coretools::str::toString;
		using coretools::instances::logfile;

		std::string ret = "[";
		for (size_t i = 0; i < numParameters(); ++i) {
			ret.append(toString(i, ": ", _betas[i], ", "));
		}
		ret.pop_back();
		ret.back() = ']';
		logfile().list(typeString(), ": ", ret);
	}
};
}

#endif
