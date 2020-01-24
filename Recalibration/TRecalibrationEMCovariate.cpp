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
		fillVectorFromStringAnySkipEmpty(str.substr(pos, pos2), values, ",");
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
		fillVectorFromStringAnySkipEmpty(type.substr(pos, pos2), args, ",");

		//extract type
		type = type.substr(0, pos);
	}
};

void TRecalibrationEMCovariate::_addPolynomialFunction(const std::string & functionString, std::vector<std::string> & args, std::vector<std::string> & values){
	if(values.empty()){
		if(args.size() != 1){
			throw "Wrong number of arguments for polynomial recal function '" + functionString + "': expect one argument (order).";
		}
		_function = new TRecalibrationEMCovariateFunction_polynomial(args[0]);
	} else {
		_function = new TRecalibrationEMCovariateFunction_polynomial(values);
	}
};

//-------------------------------------------
// TRecalibrationEMCovariate_quality
//-------------------------------------------
void TRecalibrationEMCovariate_quality::addFunction(const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == RecalModuleFunctionName_polynomial){
		_addPolynomialFunction(functionString, args, values);
	} else if(type == RecalModuleFunctionName_specific){
		std::vector<uint16_t> usedQualities;
		dataTable->fillVectorWithUsedQualities(usedQualities);
		if(values.empty()){
			_function = new TRecalibrationEMCovariateFunction_specificMap(usedQualities);
		} else {
			_function = new TRecalibrationEMCovariateFunction_specificMap(values);
		}
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

double TRecalibrationEMCovariate_quality::getEtaTerm(const TBase & base){
	return _function->getEtaTerm(base.qualityAsPhredInt);
};

double TRecalibrationEMCovariate_quality::getEtaTerm(const TRecalibrationEMReadData & data){
	return _function->getEtaTerm(data.quality);
};

//-------------------------------------------
// TRecalibrationEMCovariate_position
//-------------------------------------------
void TRecalibrationEMCovariate_position::addFunction(const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == RecalModuleFunctionName_polynomial){
		_addPolynomialFunction(functionString, args, values);
	} else if(type == RecalModuleFunctionName_specific){
		if(values.empty()){
			_function = new TRecalibrationEMCovariateFunction_specific(dataTable->maxPos);
		} else {
			_function = new TRecalibrationEMCovariateFunction_specific(values);
		}
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

double TRecalibrationEMCovariate_position::getEtaTerm(const TBase & base){
	return _function->getEtaTerm(base.distFrom5Prime);
};

double TRecalibrationEMCovariate_position::getEtaTerm(const TRecalibrationEMReadData & data){
	return _function->getEtaTerm(data.position);
};

//-------------------------------------------
// TRecalibrationEMCovariate_context
//-------------------------------------------
void TRecalibrationEMCovariate_context::addFunction(const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == RecalModuleFunctionName_specific){
		if(values.empty()){
			_function = new TRecalibrationEMCovariateFunction_specific(20);
		} else {
			_function = new TRecalibrationEMCovariateFunction_specific(values);
		}
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

double TRecalibrationEMCovariate_context::getEtaTerm(const TBase & base){
	return _function->getEtaTerm(base.context);
};

double TRecalibrationEMCovariate_context::getEtaTerm(const TRecalibrationEMReadData & data){
	return _function->getEtaTerm(data.context);
};

