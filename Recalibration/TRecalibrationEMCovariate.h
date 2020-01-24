/*
 * TRecalibrationEMModules.h
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */

#ifndef RECALIBRATION_TRECALIBRATIONEMCOVARIATE_H_
#define RECALIBRATION_TRECALIBRATIONEMCOVARIATE_H_

#include <TRecalibrationEMCovariateFunction.h>

//define covariate names
#define RecalCovariateName_none "none"
#define RecalCovariateName_quality "quality"
#define RecalCovariateName_position "position"
#define RecalCovariateName_context "context"
#define RecalCovariateName_fragmentLength "fragmentLength"

//-------------------------------------------
// TRecalibrationEMCovariate
//-------------------------------------------
class TRecalibrationEMCovariate{
protected:
	TRecalibrationEMCovariateFunction* _function;
	bool _initialized;

	void _parseModuleString(const std::string & str, std::string & type, std::vector<std::string> & args, std::vector<std::string> & values);
	void _addPolynomialFunction(const std::string & functionString, std::vector<std::string> & args, std::vector<std::string> & values);

public:
	TRecalibrationEMCovariate(){
		_initialized = false;
		_function = nullptr;
	};

	virtual ~TRecalibrationEMCovariate(){
		if(_initialized)
			delete _function;
	};

	int numParameters(){
		if(_initialized){
			return _function->numParameters();
		} else {
			return 0;
		}
	};

	virtual void addFunction(const std::string & functionString, TRecalibrationEMDataTable* dataTable){};
	virtual bool checkParameterRange(TRecalibrationEMDataTable* dataTable){ return true; };
	virtual bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos){ return true; };
	TRecalibrationEMCovariateFunction* getPointerToFunction(){ return _function;};

	virtual double getEtaTerm(const TBase & base){ return 0.0; };
	virtual double getEtaTerm(const TRecalibrationEMReadData & data){ return 0.0; };
};

//-------------------------------------------
// TRecalibrationEMCovariate_quality
//-------------------------------------------
class TRecalibrationEMCovariate_quality:public TRecalibrationEMCovariate{
public:
	TRecalibrationEMCovariate_quality(){};

	void addFunction(const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
	double getEtaTerm(const TBase & base);
	double getEtaTerm(const TRecalibrationEMReadData & data);
};

//-------------------------------------------
// TRecalibrationEMCovariate_position
//-------------------------------------------
class TRecalibrationEMCovariate_position:public TRecalibrationEMCovariate{
public:
	TRecalibrationEMCovariate_position(){};

	void addFunction(const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
	double getEtaTerm(const TBase & base);
	double getEtaTerm(const TRecalibrationEMReadData & data);
};

//-------------------------------------------
// TRecalibrationEMCovariate_context
//-------------------------------------------
class TRecalibrationEMCovariate_context:public TRecalibrationEMCovariate{
private:
	int numContext;
public:
	TRecalibrationEMCovariate_context(){
		numContext = 20;
	};

	void addFunction(const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
	double getEtaTerm(const TBase & base);
	double getEtaTerm(const TRecalibrationEMReadData & data);
};

#endif /* RECALIBRATION_TRECALIBRATIONEMCOVARIATE_H_ */
