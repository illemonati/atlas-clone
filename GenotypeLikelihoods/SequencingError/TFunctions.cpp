#include "TFunctions.h"
#include "SequencingError/TCovariate.h"
#include "TFunction.h"
#include "coretools/Main/TError.h"
#include <memory>

namespace GenotypeLikelihoods::SequencingError {
using namespace coretools::str;
namespace impl {

constexpr coretools::Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return coretools::Probability(0.9999999999);
	if (eta < -23.03) return coretools::Probability(0.0000000001);

	return coretools::logistic(eta);
}

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

template<typename Covariate> TFunction *makeCovFunction(std::string_view Function, size_t FirstParameterIndex) {
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
		case 1: fn = new TPolynomial<1, Covariate>(FirstParameterIndex); break;
		case 2: fn = new TPolynomial<2, Covariate>(FirstParameterIndex); break;
		case 3: fn = new TPolynomial<3, Covariate>(FirstParameterIndex); break;
		case 4: fn = new TPolynomial<4, Covariate>(FirstParameterIndex); break;
		case 5: fn = new TPolynomial<5, Covariate>(FirstParameterIndex); break;
		case 6: fn = new TPolynomial<6, Covariate>(FirstParameterIndex); break;
		case 7: fn = new TPolynomial<7, Covariate>(FirstParameterIndex); break;
		case 8: fn = new TPolynomial<8, Covariate>(FirstParameterIndex); break;
		case 9: fn = new TPolynomial<9, Covariate>(FirstParameterIndex); break;
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
		TFunction* fn = new TProbit<Covariate>(FirstParameterIndex);
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
	if (type == TEmpiric<Covariate>::name) {
			auto fn = new TEmpiric<Covariate>(FirstParameterIndex);
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

TFunction *makeFunction(std::string_view Covariate, std::string_view Function, size_t FirstParameterIndex) {
	// No covariate
	if (Covariate.empty()) {
		const auto [type, betas] = parseFunction(Function);
		if (type != TIntercept::name) UERROR("You must specify a covariate");
		TFunction *in = new TIntercept(FirstParameterIndex);
		if (betas.empty()) return in;
		fromString<true>(betas.front(), *in->begin());
		return in;
	}

	if (Function.empty()) {
		if (stringStartsWith(Covariate, TIntercept::name)) { // intercept can be parsed as either covariate or function
			const auto [type, betas] = parseFunction(Covariate);
			if (type != TIntercept::name) UERROR("You must specify a covariate");
			TFunction *in = new TIntercept(FirstParameterIndex);
			if (betas.empty()) return in;
			fromString<true>(betas.front(), *in->begin());
			return in;
		}
		Function = TEmpiric<TCovariate_context>::name;
	}

	if (Covariate == TCovariate_quality::name)
		return makeCovFunction<TCovariate_quality>(Function, FirstParameterIndex);
	if (Covariate == TCovariate_position::name)
		return makeCovFunction<TCovariate_position>(Function, FirstParameterIndex);
	if (Covariate == TCovariate_context::name)
		return makeCovFunction<TCovariate_context>(Function, FirstParameterIndex);
	if (Covariate == TCovariate_fragmentLength::name)
		return makeCovFunction<TCovariate_fragmentLength>(Function, FirstParameterIndex);
	if (Covariate == TCovariate_mappingQuality::name)
		return makeCovFunction<TCovariate_mappingQuality>(Function, FirstParameterIndex);

	UERROR("Covariate '", Covariate, "' does not exist!");
}

class TFunctionsTpl final : public TFunctions {
	std::vector<std::unique_ptr<TFunction>> _functions;
	std::vector<double> _oldBetas;
public:
	TFunctionsTpl(std::string_view Def) {
		// create functions
		size_t numParameters = 0;
		const auto modelDef  = impl::parseFunctions(Def);
		for (const auto &cov : modelDef) {
			_functions.emplace_back(impl::makeFunction(cov.covariate, cov.function, numParameters));

			// add new parameters
			numParameters += _functions.back()->numParameters();
		}
	}

	void checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable) override {
		size_t numParameters = 0;

		for (auto &fn : _functions) {
			if (!fn->checkOrInitValueRange(DataTable, numParameters)) {
				UERROR("Function ", fn->typeString(), " does not cover full range of data");
			}
			numParameters += fn->numParameters();
		}
	}

	size_t numParameters() const noexcept override {
		size_t numParameters = 0;
		for (const auto& f: _functions) {
			numParameters += f->numParameters();
		}
		return numParameters;
	}

	coretools::Probability getEpsilon(const BAM::TSequencedBase &base) const override {
		double eta = 0.;
		for (const auto &fn : _functions) eta += fn->getEta(base);
		return impl::calcEpsilon(eta);
	}

	coretools::Probability getEpsilon(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
								  std::vector<T2ndDerivative> &der2) const noexcept override {
		double eta = 0.;
		for (const auto &fn : _functions) eta += fn->getEta(base, der1, der2);
		return impl::calcEpsilon(eta);
	}

	void reject() noexcept override {
		auto old = _oldBetas.begin();
		for (const auto &fn : _functions) {
			for (auto &beta : *fn) {
				beta = *old;
				++old;
			}
		}
	}
	void propose(double lambda, const arma::mat &_JxF) noexcept override {
		size_t index = 0;
		_oldBetas.clear();
		for (const auto &fn : _functions) {
			for (auto &beta : *fn) {
				_oldBetas.push_back(beta);
				beta -= lambda * _JxF(index++);
			}
		}
	}
	std::string definition() const noexcept override {
		return std::accumulate(_functions.begin() + 1, _functions.end(), _functions.front()->modelString(),
							   [](auto tot, auto &f) { return tot.append(1, ';').append(f->modelString()); });
	}
	BAM::RGInfo::TInfo info() const override {
		BAM::RGInfo::TInfo in;
		for (const auto &fn : _functions) { fn->addInfo(in); }
		return in;
	}
};

class TNewFunctionsTpl final : public TFunctions {
	TIntercept _intercept;
	enum class Covariates : size_t { min = 0, Quality = min, Position, Context, FragmentLength, MappingQuality, max };
	static constexpr coretools::TStrongArray<std::string_view, Covariates> _covNames{
		{TCovariate_quality::name, TCovariate_quality::name, TCovariate_quality::name, TCovariate_quality::name,
		 TCovariate_quality::name}};

	coretools::TStrongArray<std::unique_ptr<TFunction>, Covariates> _covariates{{nullptr, nullptr, nullptr, nullptr, nullptr}};
	std::vector<double> _oldBetas;

public:
	TNewFunctionsTpl(std::string_view Def) : _intercept(0) {
		OUT(Def);
		auto modelDef  = impl::parseFunctions(Def);

		// intercept
		for (const auto &md : modelDef) {
			if (md.covariate == TIntercept::name) {
				const auto [type, betas] = parseFunction(md.function);
				if (!betas.empty()) fromString<true>(betas.front(), *_intercept.begin());
			} else if (md.function == TIntercept::name) {
				// intercept can be parsed the wrong way around :-(
				const auto [type, betas] = parseFunction(md.covariate);
				if (!betas.empty()) fromString<true>(betas.front(), *_intercept.begin());
			}
		}

		size_t index = _intercept.numParameters();

		for (auto cov = Covariates::min; cov < Covariates::max; ++cov) {
			OUT(_covNames[cov]);
			for (const auto &md : modelDef) {
				OUT(md.covariate);
				OUT(md.function);
				if (md.covariate == _covNames[cov]) {
					if (_covariates[cov]) UERROR("Covariate ", _covNames[cov], " has two functions, only one is allowed!");
					_covariates[cov].reset(impl::makeCovFunction<TCovariate_quality>(md.function, index));
					index += _covariates[cov]->numParameters();
				}
			}
			if (!_covariates[cov]) _covariates[cov] = std::make_unique<TNoFunction>();
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
		*_intercept.begin() = *old;
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
		_oldBetas.push_back(*_intercept.begin());
		*_intercept.begin() -= lambda * _JxF(index++);
		for (const auto &fn : _covariates) {
			for (auto &beta : *fn) {
				_oldBetas.push_back(beta);
				beta -= lambda * _JxF(index++);
			}
		}
	}
	std::string definition() const noexcept override {
		return std::accumulate(_covariates.begin(), _covariates.end(), _intercept.modelString(),
							   [](auto tot, auto &cov) { return tot.append(1, ';').append(cov->modelString()); });
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
	return new impl::TNewFunctionsTpl(Def);
}
} // namespace GenotypeLikelihoods::SequencingError
