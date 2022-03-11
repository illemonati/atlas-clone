/*
 * TRecalibrationEMModules.cpp
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */

#include "TSequencingErrorCovariate.h"
#include "TSequencingErrorCovariateFunction.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

namespace /*Anonymous*/ {
	auto parseModuleString(const std::string &str) {
	constexpr auto format =
		"Expected format is TYPE(ARGS)[VALUES], where [VALUES] is optional and (ARGS) is only required for some TYPE.";
	std::string type;
	std::vector<std::string> args;
	std::vector<std::string> values;

	// split string into parameters and values
	if (const auto pos = str.find('['); pos != std::string::npos) {
		// extract type
		type = str.substr(0, pos);

		// extract parameters
		const auto pos2 = str.find(']');
		if (pos == std::string::npos) { throw "Wrong format for recal function '" + str + "': missing ']'! " + format; }
		coretools::str::fillContainerFromStringAny(str.substr(pos + 1, pos2 - pos - 1), values, ",", true);
	} else {
		type = str;
	}

	// name might contain args
	if (const auto pos = type.find('('); pos != std::string::npos) {
		// extract parameters
		const auto pos2 = str.find(')');
		if (pos == std::string::npos) { throw "Wrong format for recal function '" + str + "': missing ')'! " + format; }
		coretools::str::fillContainerFromStringAny(type.substr(pos + 1, pos2 - pos - 1), args, ",", true);

		// extract type
		type = type.substr(0, pos);
	}
	return std::make_tuple(type, args, values);
}

TCovariateFunction *polynomialFunction(size_t FirstParameterIndex, const std::string &functionString,
				       const std::vector<std::string> &args, const std::vector<std::string> &values,
				       TRecalibrationEMTransformationMap *transformationMap = nullptr) {
	if (values.empty()) {
		if (args.size() != 1) {
			throw "Wrong number of arguments for polynomial recal function '" + functionString +
				"': expect one argument (order).";
		}
		auto f = new TCovariateFunction_polynomial(FirstParameterIndex,
							 coretools::str::convertStringCheck<uint32_t>(args[0]), transformationMap);
		// if no values are provided, set first beta = 1
		f->beta(0) = 1;
		return f;
	}
	return new TCovariateFunction_polynomial(FirstParameterIndex, values, transformationMap);
}

template<bool poly, bool specific, bool specificMap>
TCovariateFunction *function(const size_t FirstParameterIndex, const std::string &functionString,
			     TRecalibrationEMTransformationMap *transformationMap = nullptr) {
	const auto [type, values, args] = parseModuleString(functionString);

	// are values provided?
	if (values.empty()) {
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}
	// create function
	if constexpr (poly) {
		if (type == TCovariateFunction_polynomial::name) {
			return new TCovariateFunction_polynomial(FirstParameterIndex, values, transformationMap);
		}
	}
	if constexpr (specific) {
		if (type == TCovariateFunction_specific::name) {
			return new TCovariateFunction_specific(FirstParameterIndex, values);
		}
	}
	if constexpr (specificMap) {
		if (type == TCovariateFunction_specificMap::name) {
			return new TCovariateFunction_specificMap(FirstParameterIndex, values);
		}
	}
	throw "Recalibration function '" + type + "' not valid for covariate!";
}

template<bool poly, bool specific, bool specificMap>
TCovariateFunction *function(const size_t FirstParameterIndex, const std::string &functionString,
			     const RecalEstimatorTools::TRecalDataTable &dataTable, int v = 0,
			     TRecalibrationEMTransformationMap *transformationMap = nullptr) {
	// parse
	const auto [type, values, args] = parseModuleString(functionString);

	// create function
	if constexpr (poly) {
		if (type == TCovariateFunction_polynomial::name) {
			if (type == TCovariateFunction_polynomial::name) {
				return polynomialFunction(FirstParameterIndex, functionString, args, values, transformationMap);
			}
		}
	}
	if constexpr (specific) {
		if (type == TCovariateFunction_specific::name) {
			if (values.empty()) return new TCovariateFunction_specific(FirstParameterIndex, v);
			return new TCovariateFunction_specific(FirstParameterIndex, values);
		}
	}
	if constexpr (specificMap) {
		if (type == TCovariateFunction_specificMap::name) {
			if (values.empty())
				return new TCovariateFunction_specificMap(FirstParameterIndex, dataTable.qualities().vectorOfUsed());
			return new TCovariateFunction_specificMap(FirstParameterIndex, values);
		}
	}
	throw "Recalibration function '" + type + "' not valid for covariate quality!";
}

} // namespace

//-------------------------------------------
// TCovariate_quality
//-------------------------------------------
TCovariate_quality::TCovariate_quality(const size_t FirstParameterIndex, const std::string &functionString,
				       const RecalEstimatorTools::TRecalDataTable &dataTable) {
	_function.reset(function<true, true, false>(FirstParameterIndex, functionString, dataTable, 0, &_qualityToLogit));
}

TCovariate_quality::TCovariate_quality(const size_t FirstParameterIndex, const std::string &functionString) {
		_function.reset(function<true, true, false>(FirstParameterIndex, functionString, &_qualityToLogit));
}

bool TCovariate_quality::checkParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept {
	return _function->checkValueRange(dataTable.qualities().vectorOfUsed());
}

bool TCovariate_quality::checkParameterRange(std::vector<uint16_t> &usedQualities, uint16_t) const noexcept {
	return _function->checkValueRange(usedQualities);
}

void TCovariate_quality::adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) {
	return _function->adjustValueRanges(dataTable.qualities().vectorOfUsed());
}

//-------------------------------------------
// TCovariate_position
//-------------------------------------------
TCovariate_position::TCovariate_position(size_t FirstParameterIndex, const std::string &functionString,
					 const RecalEstimatorTools::TRecalDataTable &dataTable) {
	_function.reset(function<true, true, false>(FirstParameterIndex, functionString, dataTable, dataTable.positions().max()));
}

TCovariate_position::TCovariate_position(size_t FirstParameterIndex, const std::string &functionString) {
	_function.reset(function<true, true, false>(FirstParameterIndex, functionString));
}

bool TCovariate_position::checkParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept {
	return _function->checkValueRange(dataTable.positions().max());
}

bool TCovariate_position::checkParameterRange(std::vector<uint16_t> &, uint16_t maxPos) const noexcept {
	return _function->checkValueRange(maxPos);
}

void TCovariate_position::adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) {
	return _function->adjustValueRanges(dataTable.positions().vectorOfUsed());
}

//-------------------------------------------
// TCovariate_context
//-------------------------------------------
TCovariate_context::TCovariate_context(size_t FirstParameterIndex, const std::string &functionString,
				       const RecalEstimatorTools::TRecalDataTable &dataTable) {
	_function.reset(function<false, true, false>(FirstParameterIndex, functionString, dataTable, numContext - 1));
}

TCovariate_context::TCovariate_context(size_t FirstParameterIndex, const std::string &functionString) {
	_function.reset(function<false, true, false>(FirstParameterIndex, functionString));
}

bool TCovariate_context::checkParameterRange(const RecalEstimatorTools::TRecalDataTable &) const noexcept {
	return _function->checkValueRange(numContext);
}

bool TCovariate_context::checkParameterRange(std::vector<uint16_t> &, uint16_t) const noexcept {
	return _function->checkValueRange(numContext);
}

//-------------------------------------------
// TCovariate_fragmentLength
//-------------------------------------------

TCovariate_fragmentLength::TCovariate_fragmentLength(size_t FirstParameterIndex, const std::string &functionString,
						     const RecalEstimatorTools::TRecalDataTable &dataTable) {
	_function.reset(function<true, false, true>(FirstParameterIndex, functionString, dataTable));
}

TCovariate_fragmentLength::TCovariate_fragmentLength(size_t FirstParameterIndex, const std::string &functionString) {
	_function.reset(function<true, false, true>(FirstParameterIndex, functionString));
}

bool TCovariate_fragmentLength::checkParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept {
	return _function->checkValueRange(dataTable.fragmentLengths().vectorOfUsed());
}

bool TCovariate_fragmentLength::checkParameterRange(std::vector<uint16_t> &usedLengths, uint16_t) const noexcept {
	return _function->checkValueRange(usedLengths);
}

void TCovariate_fragmentLength::adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) {
	return _function->adjustValueRanges(dataTable.fragmentLengths().vectorOfUsed());
}

//-------------------------------------------
// TCovariate_mappingQuality
//-------------------------------------------

TCovariate_mappingQuality::TCovariate_mappingQuality(const size_t FirstParameterIndex,
						     const std::string &functionString,
						     const RecalEstimatorTools::TRecalDataTable &dataTable) {
	_function.reset(function<true, false, true>(FirstParameterIndex, functionString, dataTable));
}

TCovariate_mappingQuality::TCovariate_mappingQuality(const size_t FirstParameterIndex,
						     const std::string &functionString) {
	_function.reset(function<true, false, true>(FirstParameterIndex, functionString));
}

bool TCovariate_mappingQuality::checkParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) const noexcept {
	return _function->checkValueRange(dataTable.mappingQualities().vectorOfUsed());
}

bool TCovariate_mappingQuality::checkParameterRange(std::vector<uint16_t> &usedLengths, uint16_t) const noexcept {
	return _function->checkValueRange(usedLengths);
}

void TCovariate_mappingQuality::adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &dataTable) {
	return _function->adjustValueRanges(dataTable.mappingQualities().vectorOfUsed());
}

} // namespace SequencingError
} // namespace GenotypeLikelihoods
