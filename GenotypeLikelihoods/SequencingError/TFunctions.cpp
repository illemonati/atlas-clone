#include "TFunctions.h"
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
		if constexpr (Covariate::isIndexed) {
			auto fn = new TIndexedEmpiric<Covariate>(FirstParameterIndex);
			for (auto s : Spl) {
				TSplitter ss(s, ':');
				const auto i = fromString<size_t, true>(strip(ss.front()));
				ss.popFront();
				const auto v = fromString<double, true>(strip(ss.front()));
				ss.popFront();
				if (!ss.empty()) UERROR("Too many arguments in indexed data ", s, "!");
				fn->push_back(i, v);
			}
			return fn;
		}
		else {
			auto fn = new TEmpiric<Covariate>(FirstParameterIndex);
			for (auto s : Spl) { fn->push_back(fromString<double, true>(s)); }
			return fn;
		}
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
}

TFunctions::TFunctions(std::string_view Def) {
	// create functions
	size_t numParameters = 0;
	const auto modelDef  = impl::parseFunctions(Def);
	for (const auto &cov : modelDef) {
		_functions.emplace_back(impl::makeFunction(cov.covariate, cov.function, numParameters));

		// add new parameters
		numParameters += _functions.back()->numParameters();
	}
}
void TFunctions::checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable) {
	size_t numParameters = 0;

	for (auto &fn : _functions) {
		if (!fn->checkOrInitValueRange(DataTable, numParameters)) {
			throw "Function " + fn->typeString() + " does not cover full range of data";
		}
		numParameters += fn->numParameters();
	}
}
} // namespace GenotypeLikelihoods::SequencingError
