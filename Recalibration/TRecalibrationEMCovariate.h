/*
 * TRecalibrationEMModules.h
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */

#ifndef RECALIBRATION_TRECALIBRATIONEMCOVARIATE_H_
#define RECALIBRATION_TRECALIBRATIONEMCOVARIATE_H_

#include <memory>
#include <TRecalibrationEMCovariateFunction.h>

//define covariate names
#define RecalCovariateName_none "none"
#define RecalCovariateName_quality "quality"
#define RecalCovariateName_position "position"
#define RecalCovariateName_context "context"
#define RecalCovariateName_fragmentLength "fragmentLength"

//------------------------------------------------------------------------------------
// TRecalibrationEMCovariate
// This is the base class without any covariate. Not intended to be used in models!
//------------------------------------------------------------------------------------
class TRecalibrationEMCovariate{
protected:
	std::unique_ptr<TRecalibrationEMCovariateFunction> _function;

	void _parseModuleString(const std::string & str, std::string & type, std::vector<std::string> & args, std::vector<std::string> & values);
	void _addPolynomialFunction(const size_t FirstParameterIndex, const std::string & functionString, std::vector<std::string> & args, std::vector<std::string> & values);

	//extract
	virtual uint16_t _extractCovariate(const TBase & base){
		throw "No covariate defined for base class TRecalibrationEMCovariate!";
	};
	virtual uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		throw "No covariate defined for base class TRecalibrationEMCovariate!";
	};

public:
	TRecalibrationEMCovariate(){};
	virtual ~TRecalibrationEMCovariate(){};

	uint16_t numParameters();
	uint16_t numNonZeroFirstDerivatives();
	uint16_t numNonZeroSecondDerivatives();

	virtual std::string name(){ return RecalCovariateName_none; };

	//covariate function
	virtual void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){};
	virtual void addFunction(const size_t FirstParameterIndex, const std::string & functionString){};
	std::string functionString();
	virtual bool checkParameterRange(TRecalibrationEMDataTable* dataTable){ return true; };
	virtual bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos){ return true; };
	TRecalibrationEMCovariateFunction* getPointerToFunction(){ return _function.get(); };

	//calculate terms
	double getEtaTerm(const TBase & base){
		return _function->getEtaTerm( _extractCovariate(base) );
	};
	double getEtaTerm(const TRecalibrationEMReadData & data){
		return _function->getEtaTerm( _extractCovariate(data) );
	};

	void fillDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second){
		_function->fillDerivatives(_extractCovariate(data), first, second);
	};

	double adjustParametersPostEstimation(){
		return _function->adjustParametersPostEstimation();
	};
};

//-------------------------------------------
// TRecalibrationEMCovariate_quality
//-------------------------------------------
class TRecalibrationEMCovariate_quality:public TRecalibrationEMCovariate{
private:
	TRecalibrationEMQualityTransformationMap qualityToLogit;

	uint16_t _extractCovariate(const TBase & base){
		return base.qualityAsPhredInt;
	};
	uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		return data.qualityAsPhredInt;
	};

public:
	TRecalibrationEMCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	TRecalibrationEMCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString);

	std::string name(){ return RecalCovariateName_quality; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
	void fillFirstDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMFirstDerivatives & first, size_t & index);
	void fillSecondDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMSecondDerivatives & second, size_t & index);
};

//-------------------------------------------
// TRecalibrationEMCovariate_position
//-------------------------------------------
class TRecalibrationEMCovariate_position:public TRecalibrationEMCovariate{
private:
	uint16_t _extractCovariate(const TBase & base){
		return base.distFrom5Prime;
	};
	uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		return data.positionFrom5Prime;
	};

public:
	TRecalibrationEMCovariate_position(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	TRecalibrationEMCovariate_position(const size_t FirstParameterIndex, const std::string & functionString);

	std::string name(){ return RecalCovariateName_position; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
};

//-------------------------------------------
// TRecalibrationEMCovariate_context
//-------------------------------------------
class TRecalibrationEMCovariate_context:public TRecalibrationEMCovariate{
private:
	int numContext;

	uint16_t _extractCovariate(const TBase & base){
		return base.context;
	};
	uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		return data.context;
	};

public:
	TRecalibrationEMCovariate_context(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	TRecalibrationEMCovariate_context(const size_t FirstParameterIndex, const std::string & functionString);

	std::string name(){ return RecalCovariateName_context; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
};

//-------------------------------------------
// TRecalibrationEMCovariate_fragmentLength
//-------------------------------------------

class TRecalibrationEMCovariate_fragmentLength:public TRecalibrationEMCovariate{
private:

	uint16_t _extractCovariate(const TBase & base){
		return base.qualityAsPhredInt;
	};
	uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		return data.qualityAsPhredInt;
	};

public:
	TRecalibrationEMCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	TRecalibrationEMCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString);

	std::string name(){ return RecalCovariateName_fragmentLength; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
	void fillFirstDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMFirstDerivatives & first, size_t & index);
	void fillSecondDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMSecondDerivatives & second, size_t & index);

};

#endif /* RECALIBRATION_TRECALIBRATIONEMCOVARIATE_H_ */
