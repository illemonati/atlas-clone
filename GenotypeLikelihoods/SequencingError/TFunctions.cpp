#include "TFunctions.h"
#include "SequencingError/TCovariate.h"
#include "TFunction.h"
#include "TReadGroupInfo.h"
#include "coretools/Main/TError.h"
#include "coretools/enum.h"
#include <memory>
#include <type_traits>

namespace GenotypeLikelihoods::SequencingError {
using namespace coretools::str;
namespace impl {

auto parseFunctions(std::string_view Defs) {
	struct Voldemort {
		std::string_view covariate;
		std::string_view function;
		Voldemort(std::string_view Covariate, std::string_view Function): covariate(Covariate), function(Function) {}
	};
	std::vector<Voldemort> functions;

	for (auto def: TSplitter(Defs, ';')) {
		const auto pos  = def.find(':');
		if (pos == std::string_view::npos) functions.emplace_back(def, "");
		else functions.emplace_back(def.substr(0, pos), def.substr(pos + 1));
	}
	return functions;
}

auto parseFunction(std::string_view str) {
	const auto beg = str.find('[');
	if (beg == std::string_view::npos) { // default arguments
		return std::make_pair(str, TSplitter<>("", ','));
	}
	const auto end = str.find(']', beg);
	if (end == std::string_view::npos) {
		UERROR("Wrong format for recal function '", str, "': missing ']'! ",
			   "Expected format is TYPE[BETAS], where [BETAS] is optional.");
	}
	return std::make_pair(str.substr(0, beg), TSplitter<>(str.substr(beg + 1, end - beg - 1), ','));
}

template<typename Covariate> TFunction *makeCovFunction(std::string_view Function, size_t index) {
	auto [type, Spl] = parseFunction(Function);

	if (stringStartsWith(type, TPolynomial<1, Covariate>::name)) {
		const bool isDefault = Spl.empty();
		size_t o = 0;
		if (type.size() == TPolynomial<1, Covariate>::name.size()) {
			if (isDefault) UERROR("You must specify the order of the polynomial function or give initial values!");
		} else if (type.size() == TPolynomial<1, Covariate>::name.size() + 1) {
			o = type.back() - '0';
		} else {
			UERROR("Unknow function: ", type, "!");
		}

		std::array<double, 9> bs{};
		size_t i = 0;
		while(!Spl.empty()) {
			fromString<true>(Spl.front(), bs[i]);
			Spl.popFront();
			++i;
		}
		if (i > 0) {
			if (o == 0) o = i;
			else if (o != i) UERROR("You specified a polynomial order of ", o, " but gave ", i, " arguments!");
		}
		TFunction* fn;
		switch (o) {
		case 1: fn = new TPolynomial<1, Covariate>(index); break;
		case 2: fn = new TPolynomial<2, Covariate>(index); break;
		case 3: fn = new TPolynomial<3, Covariate>(index); break;
		case 4: fn = new TPolynomial<4, Covariate>(index); break;
		case 5: fn = new TPolynomial<5, Covariate>(index); break;
		case 6: fn = new TPolynomial<6, Covariate>(index); break;
		case 7: fn = new TPolynomial<7, Covariate>(index); break;
		case 8: fn = new TPolynomial<8, Covariate>(index); break;
		case 9: fn = new TPolynomial<9, Covariate>(index); break;
		default: UERROR("Only Polynomials from order 1 to 9 can be used!");
		}
		if (isDefault) return fn;
		i = 0;
		for (auto &beta : *fn) {
			beta = bs[i];
			++i;
		}
		return fn;
	}
	if (type == TProbit<Covariate>::name) {
		TFunction* fn = new TProbit<Covariate>(index);
		if (Spl.empty()) {
			return fn;
		}
		size_t i = 0;
		for (auto &beta : *fn) {
			if (Spl.empty()) UERROR("Not enough parameters given for function ", fn->typeString(), ". Expected ", fn->numParameters(), " got ", i, " !");
			fromString<true>(Spl.front(), beta);
			Spl.popFront();
			++i;
		}
		return fn;
	}
	if (type == TEmpiric<Covariate>::name || type=="") {
			auto fn = new TEmpiric<Covariate>(index);
			size_t size = 0;
			double back = 0.;
			for (auto s : Spl) {
				TSplitter ss(s, ':');
				const auto i = fromString<size_t, true>(strip(ss.front()));
				ss.popFront();
				const auto v = fromString<double, true>(strip(ss.front()));

				if (size == 0) { // fill first i positions with v
					for (size_t j = 0; j <= i; ++j) {
						fn->push_back(v);
					}
				}
				else { // interpolate
					const auto di   = i - (size - 1); // 1
					const auto dv   = v - back;       // v
					const auto dvdi = dv / di;
				    for (size_t j = 1; j <= di; ++j) fn->push_back(back + j * dvdi);
			    }
				back = v;
				size = i + 1;
		    }
		    return fn;
	}
	UERROR("Function '", type, "' does not exist!");
}

template<typename Covariate> TFunction *makeCovFunction(const BAM::RGInfo::TInfo& info, size_t index) {
	if (info.contains(TPolynomial<1, Covariate>::name)) {
		TFunction* fn;
		const std::vector<double>& betas = info[TPolynomial<1, Covariate>::name];
		const size_t o = betas.size();
		switch (o) {
		case 1: fn = new TPolynomial<1, Covariate>(index); break;
		case 2: fn = new TPolynomial<2, Covariate>(index); break;
		case 3: fn = new TPolynomial<3, Covariate>(index); break;
		case 4: fn = new TPolynomial<4, Covariate>(index); break;
		case 5: fn = new TPolynomial<5, Covariate>(index); break;
		case 6: fn = new TPolynomial<6, Covariate>(index); break;
		case 7: fn = new TPolynomial<7, Covariate>(index); break;
		case 8: fn = new TPolynomial<8, Covariate>(index); break;
		case 9: fn = new TPolynomial<9, Covariate>(index); break;
		default: UERROR("Only Polynomials from order 1 to 9 can be used!");
		}
		size_t i = 0;
		for (auto &beta : *fn) {
			beta = betas[i];
			++i;
		}
		return fn;
	}
	if (info.contains(TProbit<Covariate>::name)) {
		auto *fn                         = new TProbit<Covariate>(index);
		const std::vector<double> &betas = info[TProbit<Covariate>::name];
		if (betas.empty()) { return fn; }
		if (betas.size() != fn->numParameters()) {
			UERROR("Not enough parameters given for function ", fn->typeString(), ". Expected ", fn->numParameters(),
				   " got ", betas.size(), " !");
		}
		size_t i = 0;
		for (auto &beta : *fn) {
			beta = betas[i];
			++i;
		}
		return fn;
	}
	if (info.contains(TEmpiric<Covariate>::name)) {
		auto fn           = new TEmpiric<Covariate>(index);
		const auto &betas = info[TEmpiric<Covariate>::name];
		size_t size       = 0;
		double back       = 0.;
		for (const std::vector<double> &b : betas) {
			const size_t i = b.front();
			const auto v   = b.back();
			if (size == 0) { // fill first i positions with v
				for (size_t j = 0; j <= i; ++j) { fn->push_back(v); }
			} else {                              // interpolate
				const auto di   = i - (size - 1); // 1
				const auto dv   = v - back;       // v
				const auto dvdi = dv / di;
				for (size_t j = 1; j <= di; ++j) fn->push_back(back + j * dvdi);
			}
			back = v;
			size = i + 1;
		}
		return fn;
	}
	return new TNoFunction;
}

constexpr coretools::Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return coretools::Probability(0.9999999999);
	if (eta < -23.03) return coretools::Probability(0.0000000001);

	return coretools::logistic(eta);
}


class TFunctionsTpl final : public TFunctions {
	TIntercept _intercept;
	enum class Covariates : size_t { min = 0, Quality = min, Position, Context, FragmentLength, MappingQuality, max };

	coretools::TStrongArray<std::unique_ptr<TFunction>, Covariates> _covariates{{nullptr, nullptr, nullptr, nullptr, nullptr}};
	std::vector<double> _oldBetas;

	template<typename Cov, size_t I = 0> constexpr static Covariates _covIndex() noexcept {
		using Covs = std::tuple<TCovariate_quality, TCovariate_position, TCovariate_context, TCovariate_fragmentLength,
								TCovariate_mappingQuality>;
		if constexpr (std::is_same_v<Cov, std::tuple_element_t<I, Covs>>)
			return Covariates(I);
		else
			return _covIndex<Cov, I + 1>();
	}

	template<typename Cov, typename MD> size_t _makeFn(size_t index, const MD &md) {
		if (md.covariate == Cov::name) {
			constexpr auto c = _covIndex<Cov>();
			if (_covariates[c]) UERROR("Covariate ", Cov::name, " has two functions, only one is allowed!");
			_covariates[c].reset(impl::makeCovFunction<Cov>(md.function, index));
			index += _covariates[c]->numParameters();
		}
		return index;
	}

	template<typename Cov> size_t _makeFn(size_t index, const BAM::RGInfo::TInfo& info) {
		constexpr auto c = _covIndex<Cov>();
		if (info.contains(Cov::name)) {
			_covariates[c].reset(impl::makeCovFunction<Cov>(info[Cov::name], index));
		} else {
			_covariates[c].reset(new TNoFunction);
		}
		return index;
	}

public:
	TFunctionsTpl(const BAM::RGInfo::TInfo& info) {
		_intercept.beta() = info[TIntercept::name];
		size_t index        = _intercept.numParameters();
		index += _makeFn<TCovariate_quality>(index, info);
		index += _makeFn<TCovariate_context>(index, info);
		index += _makeFn<TCovariate_position>(index, info);
		index += _makeFn<TCovariate_fragmentLength>(index, info);
		index += _makeFn<TCovariate_mappingQuality>(index, info);

		for (auto &cov: _covariates) {
			if (!cov) cov.reset(new TNoFunction);
		}
	}
	TFunctionsTpl(std::string_view Def) {
		auto modelDef  = impl::parseFunctions(Def);

		// intercept
		for (const auto &md : modelDef) {
			if (md.covariate.empty()) {
				const auto [type, betas] = parseFunction(md.function);
				if (type != TIntercept::name) UERROR("You must specify a covariate for function ", type);
				if (!betas.empty()) fromString<true>(betas.front(), _intercept.beta());
			} else if (md.function.empty()) {
				// intercept can be parsed the wrong way around :-(
				const auto [type, betas] = parseFunction(md.covariate);
				if (type == TIntercept::name && !betas.empty()) fromString<true>(betas.front(), _intercept.beta());
			}
		}
		size_t index = _intercept.numParameters();

		for (const auto &md : modelDef) {
			index += _makeFn<TCovariate_quality>(index, md);
			index += _makeFn<TCovariate_context>(index, md);
			index += _makeFn<TCovariate_position>(index, md);
			index += _makeFn<TCovariate_fragmentLength>(index, md);
			index += _makeFn<TCovariate_mappingQuality>(index, md);
		}
		for (auto &cov: _covariates) {
			if (!cov) cov.reset(new TNoFunction);
		}
	}

	void checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable) override {
		size_t index = _intercept.numParameters();

		for (auto &cov : _covariates) {
			if (!cov->checkOrInitValueRange(DataTable, index)) {
				UERROR("Function ", cov->typeString(), " does not cover full range of data");
			}
			index += cov->numParameters();
		}
	}

	size_t numParameters() const noexcept override {
		size_t numParameters = _intercept.numParameters();
		for (const auto& cov: _covariates) {
			numParameters += cov->numParameters();
		}
		return numParameters;
	}

	coretools::Probability getEpsilon(const BAM::TSequencedBase &base) const override {
		double eta = _intercept.getEta();
		for (const auto &cov : _covariates) eta += cov->getEta(base);
		return impl::calcEpsilon(eta);
	}

	coretools::Probability getEpsilon(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
									  std::vector<T2ndDerivative> &der2) const noexcept override {
		double eta = _intercept.getEta(der1);
		for (const auto &cov : _covariates) eta += cov->getEta(base, der1, der2);
		return impl::calcEpsilon(eta);
	}

	void reject() noexcept override {
		auto old = _oldBetas.begin();
		_intercept.beta() = *old;
		++old;
		for (const auto &cov : _covariates) {
			for (auto &beta : *cov) {
				beta = *old;
				++old;
			}
		}
	}
	void propose(double lambda, const arma::mat &_JxF) noexcept override {
		size_t index = 0;
		_oldBetas.clear();

		_oldBetas.push_back(_intercept.beta());
		_intercept.beta() -= lambda * _JxF(index++);

		for (const auto &fn : _covariates) {
			for (auto &beta : *fn) {
				_oldBetas.push_back(beta);
				beta -= lambda * _JxF(index++);
			}
		}
	}
	std::string definition() const noexcept override {
		std::string ret = _intercept.modelString();
		for (const auto& cov: _covariates) {
			const auto ms = cov->modelString();
			if (!ms.empty()) ret.append(1, ';').append(ms);
		}
		return ret;
	}
	BAM::RGInfo::TInfo info() const override {
		BAM::RGInfo::TInfo in;
		_intercept.addInfo(in);
		for (const auto &cov : _covariates) { cov->addInfo(in); }
		return in;
	}
};

} // namespace impl

TFunctions *makeFunctions(std::string_view Def) {
	//return new impl::TFunctionsTpl(Def);
	return new impl::TFunctionsTpl(Def);
}
TFunctions *makeFunctions(const BAM::RGInfo::TInfo& info) {
	return new impl::TFunctionsTpl(info);
}
} // namespace GenotypeLikelihoods::SequencingError
