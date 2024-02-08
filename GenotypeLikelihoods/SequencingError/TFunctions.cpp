#include "TFunctions.h"
#include "coretools/Strings/stringManipulations.h"
#include "coretools/Strings/stringProperties.h"

namespace GenotypeLikelihoods::SequencingError {
using namespace coretools::str;
namespace impl {

auto parseFunctions(std::string_view Defs) {
	coretools::TStrongArray<std::string_view, Covariates, coretools::index(Covariates::max) + 1> functions{};

	for (auto def: TSplitter(Defs, ';')) {
		const auto cov = coretools::str::readBefore(def, ':');
		auto fn        = coretools::str::readAfter(def, ':');
		if (fn.empty()) fn = TEmpiric<TCovariate_context>::name;

		if (cov == TCovariate_context::name) {
			functions[TCovariate_context::index] = fn;
		} else if (cov == TCovariate_fragmentLength::name) {
			functions[TCovariate_fragmentLength::index] = fn;
		} else if (cov == TCovariate_mappingQuality::name) {
			functions[TCovariate_mappingQuality::index] = fn;
		} else if (cov == TCovariate_position::name) {
			functions[TCovariate_position::index] = fn;
		} else if (cov == TCovariate_quality::name) {
			functions[TCovariate_quality::index] = fn;
		} else if (cov == TIntercept::name) {
			functions[TCovariate_quality::index] = fn;
		} else if (fn == TIntercept::name) { // intercept can be wrong way around
			functions[TCovariate_quality::index] = cov;
		}
	}
	return functions;
}

template<typename Covariate> TFunction *makeCovFunction(const BAM::RGInfo::TInfo &RGinfo, size_t index) {
	if (!RGinfo.contains(Covariate::name)) {
			return new TNoFunction;
	}
	const auto info = RGinfo[Covariate::name];
	if (info.contains(TPolynomial<1, Covariate>::name)) {
		TFunction *fn;
		const std::vector<double> &betas = info[TPolynomial<1, Covariate>::name];
		const size_t o                   = betas.size();
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

inline auto parseFunction(std::string_view str) {
	const auto beg = str.find('[');
	if (beg == std::string_view::npos) { // default arguments
		return std::make_pair(str, coretools::str::TSplitter<>("", ','));
	}
	const auto end = str.find(']', beg);
	if (end == std::string_view::npos) {
		UERROR("Wrong format for recal function '", str, "': missing ']'! ",
			   "Expected format is TYPE[BETAS], where [BETAS] is optional.");
	}
	return std::make_pair(str.substr(0, beg), coretools::str::TSplitter<>(str.substr(beg + 1, end - beg - 1), ','));
}

template<typename Covariate> TFunction *makeCovFunction(std::string_view Function, size_t index) {
	using namespace coretools::str;

	if (Function.empty()) return new TNoFunction;

	auto [type, Spl] = parseFunction(Function);

	if (stringStartsWith(type, TPolynomial<1, Covariate>::name)) {
		const bool isDefault = Spl.empty();
		size_t o             = 0;
		if (type.size() == TPolynomial<1, Covariate>::name.size()) {
			if (isDefault) UERROR("You must specify the order of the polynomial function or give initial values!");
		} else if (type.size() == TPolynomial<1, Covariate>::name.size() + 1) {
			o = type.back() - '0';
		} else {
			UERROR("Unknow function: ", type, "!");
		}

		std::array<double, 9> bs{};
		size_t i = 0;
		while (!Spl.empty()) {
			coretools::str::fromString<true>(Spl.front(), bs[i]);
			Spl.popFront();
			++i;
		}
		if (i > 0) {
			if (o == 0)
				o = i;
			else if (o != i)
				UERROR("You specified a polynomial order of ", o, " but gave ", i, " arguments!");
		}
		TFunction *fn;
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
		TFunction *fn = new TProbit<Covariate>(index);
		if (Spl.empty()) { return fn; }
		size_t i = 0;
		for (auto &beta : *fn) {
			if (Spl.empty())
				UERROR("Not enough parameters given for function ", fn->typeString(), ". Expected ",
					   fn->numParameters(), " got ", i, " !");
			coretools::str::fromString<true>(Spl.front(), beta);
			Spl.popFront();
			++i;
		}
		return fn;
	}
	if (type == TEmpiric<Covariate>::name) {
		auto fn     = new TEmpiric<Covariate>(index);
		size_t size = 0;
		double back = 0.;
		if (Spl.empty()) fn->push_back(0.); // Initialize with at least 1 value
		for (auto s : Spl) {
			coretools::str::TSplitter ss(s, ':');
			const auto i = fromString<size_t, true>(strip(ss.front()));
			ss.popFront();
			const auto v = fromString<double, true>(strip(ss.front()));

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
	UERROR("Function '", type, "' does not exist!");
}

constexpr coretools::Probability calcEpsilon(double eta) noexcept {
	if (eta > 23.03) return coretools::Probability(0.9999999999);
	if (eta < -23.03) return coretools::Probability(0.0000000001);

	return coretools::logistic(eta);
}
} // namespace impl

TFunctions::TFunctions(const BAM::RGInfo::TInfo &info) {
	_intercept.beta() = info[TIntercept::name];
	size_t index      = _intercept.numParameters();

	using T0 = TCovariate_context;
	_covariates[T0::index].reset(impl::makeCovFunction<T0>(info, index)); 
	index += _covariates[T0::index]->numParameters();

	using T1 = TCovariate_fragmentLength;
	_covariates[T1::index].reset(impl::makeCovFunction<T1>(info, index)); 
	index += _covariates[T1::index]->numParameters();

	using T2 = TCovariate_mappingQuality;
	_covariates[T2::index].reset(impl::makeCovFunction<T2>(info, index)); 
	index += _covariates[T2::index]->numParameters();

	using T3 = TCovariate_position;
	_covariates[T3::index].reset(impl::makeCovFunction<T3>(info, index)); 
	index += _covariates[T3::index]->numParameters();

	using T4 = TCovariate_quality;
	_covariates[T4::index].reset(impl::makeCovFunction<T4>(info, index)); 
	index += _covariates[T4::index]->numParameters();
}

TFunctions::TFunctions(std::string_view Def) {
	auto fns = impl::parseFunctions(Def);
	if (fns[Covariates::max].empty()) {
		// No intercept
		_intercept.beta() = 0.;
	} else {
		fromString<true>(fns[Covariates::max], _intercept.beta());
	}
	size_t index = _intercept.numParameters();

	using T0 = TCovariate_context;
	_covariates[T0::index].reset(impl::makeCovFunction<T0>(fns[T0::index], index)); 
	index += _covariates[T0::index]->numParameters();

	using T1 = TCovariate_fragmentLength;
	_covariates[T1::index].reset(impl::makeCovFunction<T1>(fns[T1::index], index)); 
	index += _covariates[T1::index]->numParameters();

	using T2 = TCovariate_mappingQuality;
	_covariates[T2::index].reset(impl::makeCovFunction<T2>(fns[T2::index], index)); 
	index += _covariates[T2::index]->numParameters();

	using T3 = TCovariate_position;
	_covariates[T3::index].reset(impl::makeCovFunction<T3>(fns[T3::index], index)); 
	index += _covariates[T3::index]->numParameters();

	using T4 = TCovariate_quality;
	_covariates[T4::index].reset(impl::makeCovFunction<T4>(fns[T4::index], index)); 
	index += _covariates[T4::index]->numParameters();
}

void TFunctions::init(const RecalEstimatorTools::TRecalDataTable &DataTable) {
	size_t index = _intercept.numParameters();
	for (auto &cov : _covariates) {
		cov->init(DataTable, index);
		index             += cov->numParameters();
		_intercept.beta() += cov->adjust();
	}
}

size_t TFunctions::numParameters() const noexcept {
	size_t numParameters = _intercept.numParameters();
	for (const auto &cov : _covariates) { numParameters += cov->numParameters(); }
	return numParameters;
}

coretools::Probability TFunctions::getEpsilon(const BAM::TSequencedBase &base) const noexcept {
	double eta = _intercept.getEta();
	for (const auto &cov : _covariates) eta += cov->getEta(base);
	return impl::calcEpsilon(eta);
}

coretools::Probability TFunctions::getEpsilon(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
											  std::vector<T2ndDerivative> &der2) const noexcept {
	double eta = _intercept.getEta(der1);
	for (const auto &cov : _covariates) eta += cov->getEta(base, der1, der2);
	return impl::calcEpsilon(eta);
}

void TFunctions::reject() noexcept {
	auto old          = _oldBetas.begin();
	_intercept.beta() = *old;
	++old;
	for (const auto &cov : _covariates) {
		for (auto &beta : *cov) {
			beta = *old;
			++old;
		}
	}
}

void TFunctions::propose(double lambda, const arma::mat &_JxF) noexcept {
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

void TFunctions::adjust() noexcept {
	for (auto &fn : _covariates) {
		_intercept.beta() += fn->adjust();
	}
}

void TFunctions::log() const {
	_intercept.log();
	for (const auto &cov : _covariates) { cov->log(); }
}

BAM::RGInfo::TInfo TFunctions::info() const {
	BAM::RGInfo::TInfo in;
	_intercept.addInfo(in);
	for (const auto &cov : _covariates) { cov->addInfo(in); }
	return in;
}

} // namespace GenotypeLikelihoods::SequencingError
