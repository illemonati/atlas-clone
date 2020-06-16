/*
 * TRecalibrationEMModules.cpp
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */


#include "../GenotypeLikelihoods/TSequencingErrorCovariate.h"

namespace GenotypeLikelihoods{

//-------------------------------------------
// TSequencingErrorCovariate
//-------------------------------------------
void TSequencingErrorCovariate::_parseModuleString(const std::string & str, std::string & type, std::vector<std::string> & args, std::vector<std::string> & values){
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

uint16_t TSequencingErrorCovariate::numParameters(){
	if(_function){
		return _function->numParameters();
	} else {
		return 0;
	}
};

uint16_t TSequencingErrorCovariate::numNonZeroFirstDerivatives(){
	if(_function){
		return _function->numNonZeroFirstDerivatives();
	} else {
		return 0;
	}
};

uint16_t TSequencingErrorCovariate::numNonZeroSecondDerivatives(){
	if(_function){
		return _function->numNonZeroSecondDerivatives();
	} else {
		return 0;
	}
};

void TSequencingErrorCovariate::_addPolynomialFunction(const size_t FirstParameterIndex, const std::string & functionString, std::vector<std::string> & args, std::vector<std::string> & values){
	if(values.empty()){
		if(args.size() != 1){
			throw "Wrong number of arguments for polynomial recal function '" + functionString + "': expect one argument (order).";
		}
		_function.reset( new TSequencingErrorCovariateFunction_polynomial(FirstParameterIndex, convertStringCheck<uint32_t>(args[0])));
	} else {
		_function.reset( new TSequencingErrorCovariateFunction_polynomial(FirstParameterIndex, values));
	}
};

std::string TSequencingErrorCovariate::functionString(){
	if(_function){
		return _function->getModelString();
	} else return SequencingErrorCovariateFunction_none;
};

//-------------------------------------------
// TSequencingErrorCovariate_quality
//-------------------------------------------
TSequencingErrorCovariate_quality::TSequencingErrorCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TSequencingErrorCovariate_quality::TSequencingErrorCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TSequencingErrorCovariate_quality::addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
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
		std::vector<uint16_t> usedQualities;
		dataTable->fillVectorWithUsedQualities(usedQualities);
		if(values.empty()){
			_function.reset(new TSequencingErrorCovariateFunction_specificMap(FirstParameterIndex, usedQualities));
		} else {
			_function.reset(new TSequencingErrorCovariateFunction_specificMap(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	//add transformation
	_function->addTransformation(&qualityToLogit);
};

void TSequencingErrorCovariate_quality::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
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
		_function.reset(new TSequencingErrorCovariateFunction_specificMap(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}

	//add transformation
	_function->addTransformation(&qualityToLogit);
};

bool TSequencingErrorCovariate_quality::checkParameterRange(TRecalibrationEMDataTable* dataTable){
	std::vector<uint16_t> usedQualities;
	dataTable->fillVectorWithUsedQualities(usedQualities);

	return _function->checkValueRange(usedQualities);
};

bool TSequencingErrorCovariate_quality::checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos){
	return _function->checkValueRange(usedQualities);
};

//-------------------------------------------
// TSequencingErrorCovariate_position
//-------------------------------------------
TSequencingErrorCovariate_position::TSequencingErrorCovariate_position(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TSequencingErrorCovariate_position::TSequencingErrorCovariate_position(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TSequencingErrorCovariate_position::addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == SequencingErrorCovariateFunction_polynomial){
		_addPolynomialFunction(FirstParameterIndex, functionString, args, values);
	} else if(type == SequencingErrorCovariateFunction_specific){
		if(values.empty()){
			_function.reset(new TSequencingErrorCovariateFunction_specific(FirstParameterIndex, dataTable->maxPos));
		} else {
			_function.reset(new TSequencingErrorCovariateFunction_specific(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

void TSequencingErrorCovariate_position::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
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
		_function.reset(new TSequencingErrorCovariateFunction_specific(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};
bool TSequencingErrorCovariate_position::checkParameterRange(TRecalibrationEMDataTable* dataTable){
	return _function->checkValueRange(dataTable->maxPos);
};

bool TSequencingErrorCovariate_position::checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos){
	return _function->checkValueRange(maxPos);
};

//-------------------------------------------
// TSequencingErrorCovariate_context
//-------------------------------------------
TSequencingErrorCovariate_context::TSequencingErrorCovariate_context(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	numContext = 20;
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TSequencingErrorCovariate_context::TSequencingErrorCovariate_context(const size_t FirstParameterIndex, const std::string & functionString){
	numContext = 20;
	addFunction(FirstParameterIndex, functionString);
};

void TSequencingErrorCovariate_context::addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	//parse
	std::string type;
	std::vector<std::string> values, args;
	_parseModuleString(functionString, type, args, values);

	//create function
	if(type == SequencingErrorCovariateFunction_specific){
		if(values.empty()){
			_function.reset(new TSequencingErrorCovariateFunction_specific(FirstParameterIndex, 19));
		} else {
			_function.reset(new TSequencingErrorCovariateFunction_specific(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

void TSequencingErrorCovariate_context::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
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
		_function.reset(new TSequencingErrorCovariateFunction_specific(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

bool TSequencingErrorCovariate_context::checkParameterRange(TRecalibrationEMDataTable* dataTable){
	return _function->checkValueRange(20);
};

bool TSequencingErrorCovariate_context::checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos){
	return _function->checkValueRange(20);
};

//-------------------------------------------
// TSequencingErrorCovariate_fragmentLength
//-------------------------------------------

TSequencingErrorCovariate_fragmentLength::TSequencingErrorCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TSequencingErrorCovariate_fragmentLength::TSequencingErrorCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TSequencingErrorCovariate_fragmentLength::addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
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
		std::vector<uint16_t> usedLengths;
		dataTable->fillVectorWithUsedFragmentLengths(usedLengths);
		if(values.empty()){
			_function.reset(new TSequencingErrorCovariateFunction_specificMap(FirstParameterIndex, usedLengths));
		} else {
			_function.reset(new TSequencingErrorCovariateFunction_specificMap(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

void TSequencingErrorCovariate_fragmentLength::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
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
		_function.reset(new TSequencingErrorCovariateFunction_specificMap(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

bool TSequencingErrorCovariate_fragmentLength::checkParameterRange(TRecalibrationEMDataTable* dataTable){
	std::vector<uint16_t> usedLengths;
	dataTable->fillVectorWithUsedFragmentLengths(usedLengths);

	return _function->checkValueRange(usedLengths);
};

bool TSequencingErrorCovariate_fragmentLength::checkParameterRange(std::vector<uint16_t> & usedLengths, uint16_t maxPos){
	return _function->checkValueRange(usedLengths);
};

//-------------------------------------------
// TSequencingErrorCovariate_mappingQuality
//-------------------------------------------

TSequencingErrorCovariate_mappingQuality::TSequencingErrorCovariate_mappingQuality(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
	addFunction(FirstParameterIndex, functionString, dataTable);
};

TSequencingErrorCovariate_mappingQuality::TSequencingErrorCovariate_mappingQuality(const size_t FirstParameterIndex, const std::string & functionString){
	addFunction(FirstParameterIndex, functionString);
};

void TSequencingErrorCovariate_mappingQuality::addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){
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
		std::vector<uint16_t> usedMQ;
		dataTable->fillVectorWithUsedFragmentLengths(usedMQ);
		if(values.empty()){
			_function.reset(new TSequencingErrorCovariateFunction_specificMap(FirstParameterIndex, usedMQ));
		} else {
			_function.reset(new TSequencingErrorCovariateFunction_specificMap(FirstParameterIndex, values));
		}
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

void TSequencingErrorCovariate_mappingQuality::addFunction(const size_t FirstParameterIndex, const std::string & functionString){
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
		_function.reset(new TSequencingErrorCovariateFunction_specificMap(FirstParameterIndex, values));
	} else {
		throw "Recalibration function '" + type + "' not valid for covariate quality!";
	}
};

bool TSequencingErrorCovariate_mappingQuality::checkParameterRange(TRecalibrationEMDataTable* dataTable){
	std::vector<uint16_t> usedMQ;
	dataTable->fillVectorWithUsedMQ(usedMQ);

	return _function->checkValueRange(usedMQ);
};

bool TSequencingErrorCovariate_mappingQuality::checkParameterRange(std::vector<uint16_t> & usedLengths, uint16_t maxPos){
	return _function->checkValueRange(usedLengths);
};

//end namespace recal
};
