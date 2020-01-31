/*
 * TRecalibrationEMModules.cpp
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */


#include <TRecalibrationEMCovariate.h>

//-------------------------------------------
// TRecalibrationEMCovariate
//-------------------------------------------
void TRecalibrationEMCovariate::_parseModuleString(const std::string & str, std::string & type, std::vector<std::string> & args, std::vector<std::string> & values){
	std::string format = "Expected format is TYPE(ARGS)[VALUES], where [VALUES] is optional and (ARGS) is only required for TYPE polynomial.";

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
		fillVectorFromStringAnySkipEmpty(str.substr(pos+1, pos2-pos-1), values, ",");
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
		fillVectorFromStringAnySkipEmpty(type.substr(pos+1, pos2-pos-1), args, ",");

		//extract type
		type = type.substr(0, pos);
	}
};

size_t TRecalibrationEMCovariate::numParameters(){
	if(_initialized){
		return _function->numParameters();
	} else {
		return 0;
	}
};

size_t TRecalibrationEMCovariate::numNonZeroFirstDerivatives(){
	if(_initialized){
		return _function->numNonZeroFirstDerivatives();
	} else {
		return 0;
	}
};

size_t TRecalibrationEMCovariate::numNonZeroSecondDerivatives(){
	if(_initialized){
		return _function->numNonZeroSecondDerivatives();
	} else {
		return 0;
	}
};

void TRecalibrationEMCovariate::_addPolynomialFunction(const size_t FirstParameterIndex, const std::string & functionString, std::vector<std::string> & args, std::vector<std::string> & values){
	if(values.empty()){
		if(args.size() != 1){
			throw "Wrong number of arguments for polynomial recal function '" + functionString + "': expect one argument (order).";
		}
		_function = new TRecalibrationEMCovariateFunction_polynomial(FirstParameterIndex, stringToIntCheck(args[0]));
	} else {
		_function = new TRecalibrationEMCovariateFunction_polynomial(FirstParameterIndex, values);
	}
};

std::string TRecalibrationEMCovariate::functionString(){
	if(_initialized){
		return _function->getModelString();
	} else return RecalModuleFunctionName_none;
};

//-------------------------------------------
// TRecalibrationEMCovariate_quality
//-------------------------------------------
TRecalibrationEMCovariate_quality::TRecalibrationEMCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TRecalibrationEMCovariate_quality::TRecalibrationEMCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TRecalibrationEMCovariate_quality::addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	//clear if already initialized
	clear();

	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == RecalModuleFunctionName_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == RecalModuleFunctionName_specific){
		std::vector<uint16_t> usedQualities;
		dataTable->fillVectorWithUsedQualities(usedQualities);
		if(values.empty()){
			_function = new TRecalibrationEMCovariateFunction_specificMap(FirstParameterIndex, usedQualities);
		} else {
			_function = new TRecalibrationEMCovariateFunction_specificMap(FirstParameterIndex, values);
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	_initialized = true;
};

void TRecalibrationEMCovariate_quality::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
	//clear if already initialized
	clear();

	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//are values provided?
	if(values.empty()){
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	//create function
	if(type == RecalModuleFunctionName_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == RecalModuleFunctionName_specific){
		_function = new TRecalibrationEMCovariateFunction_specificMap(FirstParameterIndex, values);
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	_initialized = true;
};

bool TRecalibrationEMCovariate_quality::checkParameterRange(TRecalibrationEMDataTable* dataTable){
	std::vector<uint16_t> usedQualities;
	dataTable->fillVectorWithUsedQualities(usedQualities);

	return _function->checkValueRange(usedQualities);
};

bool TRecalibrationEMCovariate_quality::checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos){
	return _function->checkValueRange(usedQualities);
};

//-------------------------------------------
// TRecalibrationEMCovariate_position
//-------------------------------------------
TRecalibrationEMCovariate_position::TRecalibrationEMCovariate_position(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TRecalibrationEMCovariate_position::TRecalibrationEMCovariate_position(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TRecalibrationEMCovariate_position::addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	//clear if already initialized
	clear();

	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == RecalModuleFunctionName_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == RecalModuleFunctionName_specific){
		if(values.empty()){
			_function = new TRecalibrationEMCovariateFunction_specific(FirstParameterIndex, dataTable->maxPos);
		} else {
			_function = new TRecalibrationEMCovariateFunction_specific(FirstParameterIndex, values);
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	_initialized = true;
};

void TRecalibrationEMCovariate_position::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
	//clear if already initialized
	clear();

	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//are values provided?
	if(values.empty()){
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	//create function
	if(type == RecalModuleFunctionName_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == RecalModuleFunctionName_specific){
		_function = new TRecalibrationEMCovariateFunction_specific(FirstParameterIndex, values);
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	_initialized = true;
};
bool TRecalibrationEMCovariate_position::checkParameterRange(TRecalibrationEMDataTable* dataTable){
	return _function->checkValueRange(dataTable->maxPos);
};

bool TRecalibrationEMCovariate_position::checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos){
	return _function->checkValueRange(maxPos);
};

//-------------------------------------------
// TRecalibrationEMCovariate_context
//-------------------------------------------
TRecalibrationEMCovariate_context::TRecalibrationEMCovariate_context(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	numContext = 20;
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TRecalibrationEMCovariate_context::TRecalibrationEMCovariate_context(const size_t FirstParameterIndex, const std::string & functionString){
	numContext = 20;
	addFunction(FirstParameterIndex, functionString);
};

void TRecalibrationEMCovariate_context::addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	//clear if already initialized
	clear();

	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == RecalModuleFunctionName_specific){
		if(values.empty()){
			_function = new TRecalibrationEMCovariateFunction_specific(FirstParameterIndex, 20);
		} else {
			_function = new TRecalibrationEMCovariateFunction_specific(FirstParameterIndex, values);
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	_initialized = true;
};

void TRecalibrationEMCovariate_context::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
	//clear if already initialized
	clear();

	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//are values provided?
	if(values.empty()){
		throw "Failed to initialize recalibration covariate: missing [VALUES] in '" + functionString + "'!";
	}

	//create function
	if(type == RecalModuleFunctionName_specific){
		_function = new TRecalibrationEMCovariateFunction_specific(FirstParameterIndex, values);
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	_initialized = true;
};

bool TRecalibrationEMCovariate_context::checkParameterRange(TRecalibrationEMDataTable* dataTable){
	return _function->checkValueRange(20);
};

bool TRecalibrationEMCovariate_context::checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos){
	return _function->checkValueRange(20);
};

