/*
 * TRecalibrationEMModules.cpp
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */


#include "TSequencingErrorCovariate.h"

namespace GenotypeLikelihoods{
namespace SequencingError {

//-------------------------------------------
// TCovariate
//-------------------------------------------
void TCovariate::_parseModuleString(const std::string & str, std::string & type, std::vector<std::string> & args, std::vector<std::string> & values){
	std::string format = "Expected format is TYPE(ARGS)[VALUES], where [VALUES] is optional and (ARGS) is only required for some TYPE.";

	//empty
	values.clear();
	args.clear();

	//split string into parameters and values
	size_t pos = str.find('[');
	if(pos != std::string::npos){
		//extract type
		type = str.substr(0, pos);

		//extract parameters
		size_t pos2 = str.find(']');
		if(pos == std::string::npos){
			throw "Wrong format for recal function '" + str + "': missing ']'! " + format;
		}
		coretools::str::fillContainerFromStringAny(str.substr(pos+1, pos2-pos-1), values, ",", true);
	} else {
		type = str;
	}

	//name might contain args
	pos = type.find('(');
	if(pos != std::string::npos){
		//extract parameters
		size_t pos2 = str.find(')');
		if(pos == std::string::npos){
			throw "Wrong format for recal function '" + str + "': missing ')'! " + format;
		}
		coretools::str::fillContainerFromStringAny(type.substr(pos+1, pos2-pos-1), args, ",", true);

		//extract type
		type = type.substr(0, pos);
	}
};

uint16_t TCovariate::numParameters(){
	if(_function){
		return _function->numParameters();
	} else {
		return 0;
	}
};

uint16_t TCovariate::numNonZeroFirstDerivatives(){
	if(_function){
		return _function->numNonZeroFirstDerivatives();
	} else {
		return 0;
	}
};

uint16_t TCovariate::numNonZeroSecondDerivatives(){
	if(_function){
		return _function->numNonZeroSecondDerivatives();
	} else {
		return 0;
	}
};

void TCovariate::_addPolynomialFunction(const size_t FirstParameterIndex, const std::string & functionString, std::vector<std::string> & args, std::vector<std::string> & values){
	if(values.empty()){
		if(args.size() != 1){
			throw "Wrong number of arguments for polynomial recal function '" + functionString + "': expect one argument (order).";
		}
		_function.reset( new TCovariateFunction_polynomial(FirstParameterIndex, coretools::str::convertStringCheck<uint32_t>(args[0])));
	} else {
		_function.reset( new TCovariateFunction_polynomial(FirstParameterIndex, values));
	}
};

std::string TCovariate::functionString(){
	if(_function){
		return _function->getModelString();
	} else return SequencingErrorCovariateFunction_none;
};

//-------------------------------------------
// TCovariate_quality
//-------------------------------------------
TCovariate_quality::TCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TCovariate_quality::TCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TCovariate_quality::addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == SequencingErrorCovariateFunction_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);

		//if no values are provided, set first beta = 1
		if(values.empty()){
			_function->setBeta(0, 1.0);
		}
	} else if(type == SequencingErrorCovariateFunction_specific){
		std::vector<uint16_t> usedQualities = dataTable.qualities().vectorOfUsed();
		if(values.empty()){
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, usedQualities));
		} else {
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	//add transformation
	_function->addTransformation(&qualityToLogit);
};

void TCovariate_quality::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//are values provided?
	if(values.empty()){
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	//create function
	if(type == SequencingErrorCovariateFunction_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == SequencingErrorCovariateFunction_specific){
		_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	//add transformation
	_function->addTransformation(&qualityToLogit);
};

bool TCovariate_quality::checkParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable){
	std::vector<uint16_t> usedQualities = dataTable.qualities().vectorOfUsed();
	return _function->checkValueRange(usedQualities);
};

bool TCovariate_quality::checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t){
	return _function->checkValueRange(usedQualities);
};

void TCovariate_quality::adjustParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable){
	return _function->adjustValueRanges(dataTable.qualities().vectorOfUsed());
};


//-------------------------------------------
// TCovariate_position
//-------------------------------------------
TCovariate_position::TCovariate_position(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TCovariate_position::TCovariate_position(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TCovariate_position::addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == SequencingErrorCovariateFunction_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == SequencingErrorCovariateFunction_specific){
		if(values.empty()){
			_function.reset(new TCovariateFunction_specific(FirstParameterIndex, dataTable.positions().max()));
		} else {
			_function.reset(new TCovariateFunction_specific(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

void TCovariate_position::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//are values provided?
	if(values.empty()){
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	//create function
	if(type == SequencingErrorCovariateFunction_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == SequencingErrorCovariateFunction_specific){
		_function.reset(new TCovariateFunction_specific(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};
bool TCovariate_position::checkParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable){
	return _function->checkValueRange(dataTable.positions().max());
};

bool TCovariate_position::checkParameterRange(std::vector<uint16_t> &, uint16_t maxPos){
	return _function->checkValueRange(maxPos);
};

void TCovariate_position::adjustParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable){
	return _function->adjustValueRanges(dataTable.positions().vectorOfUsed());
};

//-------------------------------------------
// TCovariate_context
//-------------------------------------------
TCovariate_context::TCovariate_context(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable){
	numContext = 20;
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TCovariate_context::TCovariate_context(const size_t FirstParameterIndex, const std::string & functionString){
	numContext = 20;
	addFunction(FirstParameterIndex, functionString);
};

void TCovariate_context::addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable &){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == SequencingErrorCovariateFunction_specific){
		if(values.empty()){
			_function.reset(new TCovariateFunction_specific(FirstParameterIndex, 19));
		} else {
			_function.reset(new TCovariateFunction_specific(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

void TCovariate_context::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//are values provided?
	if(values.empty()){
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	//create function
	if(type == SequencingErrorCovariateFunction_specific){
		_function.reset(new TCovariateFunction_specific(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

bool TCovariate_context::checkParameterRange(const RecalEstimatorTools::TRecalDataTable &){
	return _function->checkValueRange(20);
};

bool TCovariate_context::checkParameterRange(std::vector<uint16_t> &, uint16_t){
	return _function->checkValueRange(20);
};

//-------------------------------------------
// TCovariate_fragmentLength
//-------------------------------------------

TCovariate_fragmentLength::TCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TCovariate_fragmentLength::TCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TCovariate_fragmentLength::addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == SequencingErrorCovariateFunction_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);

		//if no values are provided, set first beta = 1
		if(values.empty()){
			_function->setBeta(0, 1.0);
		}
	} else if(type == SequencingErrorCovariateFunction_specific){
		std::vector<uint16_t> usedLengths = dataTable.fragmentLengths().vectorOfUsed();
		if(values.empty()){
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, usedLengths));
		} else {
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

void TCovariate_fragmentLength::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//are values provided?
	if(values.empty()){
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	//create function
	if(type == SequencingErrorCovariateFunction_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == SequencingErrorCovariateFunction_specific){
		_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

bool TCovariate_fragmentLength::checkParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable){
	std::vector<uint16_t> usedLengths = dataTable.fragmentLengths().vectorOfUsed();
	return _function->checkValueRange(usedLengths);
};

bool TCovariate_fragmentLength::checkParameterRange(std::vector<uint16_t> & usedLengths, uint16_t){
	return _function->checkValueRange(usedLengths);
};

void TCovariate_fragmentLength::adjustParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable){
	return _function->adjustValueRanges(dataTable.fragmentLengths().vectorOfUsed());
};

//-------------------------------------------
// TCovariate_mappingQuality
//-------------------------------------------

TCovariate_mappingQuality::TCovariate_mappingQuality(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TCovariate_mappingQuality::TCovariate_mappingQuality(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TCovariate_mappingQuality::addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == SequencingErrorCovariateFunction_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);

		//if no values are provided, set first beta = 1
		if(values.empty()){
			_function->setBeta(0, 1.0);
		}
	} else if(type == SequencingErrorCovariateFunction_specific){
		std::vector<uint16_t> usedMQ = dataTable.mappingQualities().vectorOfUsed();
		if(values.empty()){
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, usedMQ));
		} else {
			_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

void TCovariate_mappingQuality::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//are values provided?
	if(values.empty()){
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	//create function
	if(type == SequencingErrorCovariateFunction_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == SequencingErrorCovariateFunction_specific){
		_function.reset(new TCovariateFunction_specificMap(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

bool TCovariate_mappingQuality::checkParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable){
	std::vector<uint16_t> usedMQ = dataTable.mappingQualities().vectorOfUsed();
	return _function->checkValueRange(usedMQ);
};

bool TCovariate_mappingQuality::checkParameterRange(std::vector<uint16_t> & usedLengths, uint16_t){
	return _function->checkValueRange(usedLengths);
};

void TCovariate_mappingQuality::adjustParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable){
	return _function->adjustValueRanges(dataTable.mappingQualities().vectorOfUsed());
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods
