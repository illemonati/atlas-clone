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
void parseModuleString(const std::string &str, std::string &type, std::vector<std::string> &args,
		       std::vector<std::string> &values) {
	std::string format =
		"Expected format is TYPE(ARGS)[VALUES], where [VALUES] is optional and (ARGS) is only required for some TYPE.";

	// empty
	values.clear();
	args.clear();

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
}

TCovariateFunction *polynomialFunction(size_t FirstParameterIndex, const std::string &functionString,
					  std::vector<std::string> &args, std::vector<std::string> &values, TRecalibrationEMTransformationMap* transformationMap = nullptr) {
	if (values.empty()) {
		if (args.size() != 1) {
			throw "Wrong number of arguments for polynomial recal function '" + functionString +
				"': expect one argument (order).";
		}
		return new TCovariateFunction_polynomial(FirstParameterIndex,
							 coretools::str::convertStringCheck<uint32_t>(args[0]), transformationMap);
	}
	return new TCovariateFunction_polynomial(FirstParameterIndex, values, transformationMap);
}

template<size_t N>
TCovariateFunction *function(size_t FirstParameterIndex, const std::string &functionString, const std::array<std::string, N> & allowed) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// Is this function allowed?
	if (std::find(allowed.begin(), allowed.end(), type) == allowed.end())
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
		

	// are values provided?
	if (values.empty()) 
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";

	// create function
	if (type == TCovariateFunction_polynomial::name) {
		return polynomialFunction(FirstParameterIndex, functionString, args, values);
	}
	if (type == TCovariateFunction_specific::name) {
		return new TCovariateFunction_specific(FirstParameterIndex, values);
	} 

	throw "Recalibration function '" + type + "' not valid for covariate quality!";
}

} // namespace

//-------------------------------------------
// TCovariate_quality
//-------------------------------------------
TCovariate_quality::TCovariate_quality(const size_t FirstParameterIndex, const std::string &functionString,
				       const RecalEstimatorTools::TRecalDataTable &dataTable) {
		addFunction(FirstParameterIndex, functionString, dataTable);
}

TCovariate_quality::TCovariate_quality(const size_t FirstParameterIndex, const std::string &functionString) {
	_function.reset(function(FirstParameterIndex, functionString, std::array<std::string, 1>{TCovariateFunction_specific::name}));
}

void TCovariate_quality::addFunction(const size_t FirstParameterIndex, const std::string &functionString,
				     const RecalEstimatorTools::TRecalDataTable &dataTable) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// create function
	if (type == TCovariateFunction_polynomial::name) {
		_function.reset(polynomialFunction(FirstParameterIndex, functionString, args, values, &_qualityToLogit));

		// if no values are provided, set first beta = 1
		if (values.empty()) _function->beta(0) = 1.0;
	} else if (type == TCovariateFunction_specific::name) {
		if (values.empty()) {
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, dataTable.qualities().vectorOfUsed()));
		} else {
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
}

void TCovariate_quality::addFunction(const size_t FirstParameterIndex, const std::string &functionString) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// are values provided?
	if (values.empty()) {
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	// create function
	if (type == TCovariateFunction_polynomial::name) {
		_function.reset(polynomialFunction(FirstParameterIndex, functionString, args, values, &_qualityToLogit));
	} else if (type == TCovariateFunction_specific::name) {
		_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
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
	addFunction(FirstParameterIndex, functionString, dataTable);
}

TCovariate_position::TCovariate_position(size_t FirstParameterIndex, const std::string &functionString) {
	addFunction(FirstParameterIndex, functionString);
}

void TCovariate_position::addFunction(size_t FirstParameterIndex, const std::string &functionString,
				      const RecalEstimatorTools::TRecalDataTable &dataTable) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// create function
	if (type == TCovariateFunction_polynomial::name) {
		_function.reset(polynomialFunction(FirstParameterIndex, functionString, args, values));
	} else if (type == TCovariateFunction_specific::name) {
		if (values.empty()) {
			_function.reset(new TCovariateFunction_specific(FirstParameterIndex, dataTable.positions().max()));
		} else {
			_function.reset(new TCovariateFunction_specific(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
}

void TCovariate_position::addFunction(size_t FirstParameterIndex, const std::string &functionString) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// are values provided?
	if (values.empty()) {
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	// create function
	if (type == TCovariateFunction_polynomial::name) {
		_function.reset(polynomialFunction(FirstParameterIndex, functionString, args, values));
	} else if (type == TCovariateFunction_specific::name) {
		_function.reset(new TCovariateFunction_specific(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
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
	addFunction(FirstParameterIndex, functionString, dataTable);
}

TCovariate_context::TCovariate_context(size_t FirstParameterIndex, const std::string &functionString) {
	addFunction(FirstParameterIndex, functionString);
}

void TCovariate_context::addFunction(size_t FirstParameterIndex, const std::string &functionString,
				     const RecalEstimatorTools::TRecalDataTable &) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// create function
	if (type == TCovariateFunction_specific::name) {
		if (values.empty()) {
			_function.reset(new TCovariateFunction_specific(FirstParameterIndex, numContext - 1));
		} else {
			_function.reset(new TCovariateFunction_specific(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
}

void TCovariate_context::addFunction(size_t FirstParameterIndex, const std::string &functionString) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// are values provided?
	if (values.empty()) {
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	// create function
	if (type == TCovariateFunction_specific::name) {
		_function.reset(new TCovariateFunction_specific(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
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
	addFunction(FirstParameterIndex, functionString, dataTable);
}

TCovariate_fragmentLength::TCovariate_fragmentLength(size_t FirstParameterIndex, const std::string &functionString) {
	addFunction(FirstParameterIndex, functionString);
}

void TCovariate_fragmentLength::addFunction(const size_t FirstParameterIndex, const std::string &functionString,
					    const RecalEstimatorTools::TRecalDataTable &dataTable) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// create function
	if (type == TCovariateFunction_polynomial::name) {
		_function.reset(polynomialFunction(FirstParameterIndex, functionString, args, values));

		// if no values are provided, set first beta = 1
		if (values.empty()) { _function->beta(0) = 1.0; }
	} else if (type == TCovariateFunction_specific::name) {
		std::vector<uint16_t> usedLengths = dataTable.fragmentLengths().vectorOfUsed();
		if (values.empty()) {
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, usedLengths));
		} else {
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
}

void TCovariate_fragmentLength::addFunction(size_t FirstParameterIndex, const std::string &functionString) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// are values provided?
	if (values.empty()) {
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	// create function
	if (type == TCovariateFunction_polynomial::name) {
		_function.reset(polynomialFunction(FirstParameterIndex, functionString, args, values));
	} else if (type == TCovariateFunction_specific::name) {
		_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
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
	addFunction(FirstParameterIndex, functionString, dataTable);
}

TCovariate_mappingQuality::TCovariate_mappingQuality(const size_t FirstParameterIndex,
						     const std::string &functionString) {
	addFunction(FirstParameterIndex, functionString);
}

void TCovariate_mappingQuality::addFunction(const size_t FirstParameterIndex, const std::string &functionString,
					    const RecalEstimatorTools::TRecalDataTable &dataTable) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// create function
	if (type == TCovariateFunction_polynomial::name) {
		_function.reset(polynomialFunction(FirstParameterIndex, functionString, args, values));

		// if no values are provided, set first beta = 1
		if (values.empty()) _function->beta(0) = 1.0;
	} else if (type == TCovariateFunction_specific::name) {
		std::vector<uint16_t> usedMQ = dataTable.mappingQualities().vectorOfUsed();
		if (values.empty()) {
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, usedMQ));
		} else {
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
}

void TCovariate_mappingQuality::addFunction(const size_t FirstParameterIndex, const std::string &functionString) {
	// parse
	std::string type;
	std::vector<std::string> values, args;
	parseModuleString(functionString, type, args, values);

	// are values provided?
	if (values.empty()) {
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	// create function
	if (type == TCovariateFunction_polynomial::name) {
		_function.reset(polynomialFunction(FirstParameterIndex, functionString, args, values));
	} else if (type == TCovariateFunction_specific::name) {
		_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
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
