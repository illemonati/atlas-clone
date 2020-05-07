/*
 * TRecalibrationEMModules.h
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_

#include <memory>
#include "../GenotypeLikelihoods/TSequencingErrorCovariateFunction.h"

namespace GenotypeLikelihoods{

//define covariate names
#define SequencingErrorCovariateName_none "none"
#define SequencingErrorCovariateName_quality "quality"
#define SequencingErrorCovariateName_position "position"
#define SequencingErrorCovariateName_context "context"
#define SequencingErrorCovariateName_fragmentLength "fragmentLength"
#define SequencingErrorCovariateName_mappingQuality "mappingQuality"

//------------------------------------------------------------------------------------
// TSequencingErrorCovariate
// This is the base class without any covariate. Not intended to be used in models!
//------------------------------------------------------------------------------------
class TSequencingErrorCovariate{
protected:
	std::unique_ptr<TSequencingErrorCovariateFunction> _function;

	void _parseModuleString(const std::string & str, std::string & type, std::vector<std::string> & args, std::vector<std::string> & values);
	void _addPolynomialFunction(const size_t FirstParameterIndex, const std::string & functionString, std::vector<std::string> & args, std::vector<std::string> & values);

	//extract
	virtual uint16_t _extractCovariate(const TBaseData & base){
		throw "No covariate defined for base class TRecalibrationEMCovariate!";
	};
	virtual uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		throw "No covariate defined for base class TRecalibrationEMCovariate!";
	};

public:
	TSequencingErrorCovariate(){};
	virtual ~TSequencingErrorCovariate(){};

	uint16_t numParameters();
	uint16_t numNonZeroFirstDerivatives();
	uint16_t numNonZeroSecondDerivatives();

	virtual std::string name(){ return SequencingErrorCovariateName_none; };

	//covariate function
	virtual void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable){};
	virtual void addFunction(const size_t FirstParameterIndex, const std::string & functionString){};
	std::string functionString();
	virtual bool checkParameterRange(TRecalibrationEMDataTable* dataTable){ return true; };
	virtual bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos){ return true; };
	TSequencingErrorCovariateFunction* getPointerToFunction(){ return _function.get(); };

	//calculate terms
	double getEtaTerm(const TBaseData & base){
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
// TSequencingErrorCovariate_quality
//-------------------------------------------
class TSequencingErrorCovariate_quality:public TSequencingErrorCovariate{
private:
	TRecalibrationEMQualityTransformationMap qualityToLogit;

	uint16_t _extractCovariate(const TBaseData & base){
		return base.qualityAsPhredInt;
	};
	uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		return data.qualityAsPhredInt;
	};

public:
	TSequencingErrorCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	TSequencingErrorCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString);

	std::string name(){ return SequencingErrorCovariateName_quality; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
	void fillFirstDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMFirstDerivatives & first, size_t & index);
	void fillSecondDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMSecondDerivatives & second, size_t & index);
};

//-------------------------------------------
// TSequencingErrorCovariate_position
//-------------------------------------------
class TSequencingErrorCovariate_position:public TSequencingErrorCovariate{
private:
	uint16_t _extractCovariate(const TBaseData & base){
		return base.distFrom5Prime;
	};
	uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		return data.positionFrom5Prime;
	};

public:
	TSequencingErrorCovariate_position(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	TSequencingErrorCovariate_position(const size_t FirstParameterIndex, const std::string & functionString);

	std::string name(){ return SequencingErrorCovariateName_position; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
};

//-------------------------------------------
// TSequencingErrorCovariate_context
//-------------------------------------------
class TSequencingErrorCovariate_context:public TSequencingErrorCovariate{
private:
	int numContext;

	uint16_t _extractCovariate(const TBaseData & base){
		return base.context;
	};
	uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		return data.context;
	};

public:
	TSequencingErrorCovariate_context(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	TSequencingErrorCovariate_context(const size_t FirstParameterIndex, const std::string & functionString);

	std::string name(){ return SequencingErrorCovariateName_context; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
};

//-------------------------------------------
// TSequencingErrorCovariate_fragmentLength
//-------------------------------------------

class TSequencingErrorCovariate_fragmentLength:public TSequencingErrorCovariate{
private:

	uint16_t _extractCovariate(const TBaseData & base){
		return base.fragmentLength;
	};
	uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		return data.fragmentLength;
	};

public:
	TSequencingErrorCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	TSequencingErrorCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString);

	std::string name(){ return SequencingErrorCovariateName_fragmentLength; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
	void fillFirstDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMFirstDerivatives & first, size_t & index);
	void fillSecondDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMSecondDerivatives & second, size_t & index);

};

//-------------------------------------------
// TSequencingErrorCovariate_mappingQuality
//-------------------------------------------

class TSequencingErrorCovariate_mappingQuality:public TSequencingErrorCovariate{
private:

	uint16_t _extractCovariate(const TBaseData & base){
		return base.mappingQuality;
	};
	uint16_t _extractCovariate(const TRecalibrationEMReadData & data){
		return data.mappingQuality;
	};

public:
	TSequencingErrorCovariate_mappingQuality(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	TSequencingErrorCovariate_mappingQuality(const size_t FirstParameterIndex, const std::string & functionString);

	std::string name(){ return SequencingErrorCovariateName_fragmentLength; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, TRecalibrationEMDataTable* dataTable);
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString);
	bool checkParameterRange(TRecalibrationEMDataTable* dataTable);
	bool checkParameterRange(std::vector<uint16_t> & usedQualities, uint16_t maxPos);
	void fillFirstDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMFirstDerivatives & first, size_t & index);
	void fillSecondDerivatives(const TRecalibrationEMReadData & data, TRecalibrationEMSecondDerivatives & second, size_t & index);

};

}; //end namespace recal

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_ */
