/*
 * TRecalibrationEMModules.h
 *
 *  Created on: Jan 24, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_

#include <memory>
#include "../BAM/TSequencedBase.h"
#include "GenotypeTypes.h"
#include "TSequencingErrorCovariateFunction.h"
#include "RecalEstimatorTools.h"

namespace GenotypeLikelihoods{

//define covariate names

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
	virtual uint16_t _extractCovariate(const BAM::TSequencedBase &) = 0;

public:
	static inline const std::string name = "none";
	TSequencingErrorCovariate(){};
	virtual ~TSequencingErrorCovariate(){};

	uint16_t numParameters();
	uint16_t numNonZeroFirstDerivatives();
	uint16_t numNonZeroSecondDerivatives();

	virtual std::string typeString() const = 0;

	//covariate function
	virtual void addFunction(const size_t, const std::string &, const RecalEstimatorTools::TRecalDataTable &) = 0;
	virtual void addFunction(const size_t, const std::string &) = 0;
	std::string functionString();

	virtual bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable &) = 0;
	virtual bool checkParameterRange(std::vector<uint16_t> &, uint16_t) = 0;
	virtual void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &) = 0;

	TSequencingErrorCovariateFunction* getPointerToFunction(){ return _function.get(); };

	//calculate terms
	double getEtaTerm(const BAM::TSequencedBase & base){
		return _function->getEtaTerm( _extractCovariate(base) );
	}

	void fillDerivatives(const BAM::TSequencedBase & base, TRecalibrationEMFirstDerivatives & first, TRecalibrationEMSecondDerivatives & second){
		_function->fillDerivatives(_extractCovariate(base), first, second);
	}

	double adjustParametersPostEstimation(){
		return _function->adjustParametersPostEstimation();
	}
};

//-------------------------------------------
// TSequencingErrorCovariate_quality
//-------------------------------------------
class TSequencingErrorCovariate_quality:public TSequencingErrorCovariate{
private:
	TRecalibrationEMQualityTransformationMap qualityToLogit;

	uint16_t _extractCovariate(const BAM::TSequencedBase & base) override {
		return base.originalQuality_phredInt.get();
	};

public:
	static inline const std::string name = "quality";
	TSequencingErrorCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable);
	TSequencingErrorCovariate_quality(const size_t FirstParameterIndex, const std::string & functionString);

	std::string typeString() const override { return name; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable) override;
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable) override;
	bool checkParameterRange(std::vector<uint16_t> & usedValues, uint16_t maxPos) override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable) override;
};

//-------------------------------------------
// TSequencingErrorCovariate_position
//-------------------------------------------
class TSequencingErrorCovariate_position:public TSequencingErrorCovariate{
private:
	uint16_t _extractCovariate(const BAM::TSequencedBase & base) override {
		return base.distFrom5Prime;
	};

public:
	static inline const std::string name =  "position";
	TSequencingErrorCovariate_position(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable);
	TSequencingErrorCovariate_position(const size_t FirstParameterIndex, const std::string & functionString);

	std::string typeString() const override  { return name; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable)  override;
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable) override;
	bool checkParameterRange(std::vector<uint16_t> & usedValues, uint16_t maxPos) override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable) override;
};

//-------------------------------------------
// TSequencingErrorCovariate_context
//-------------------------------------------
class TSequencingErrorCovariate_context:public TSequencingErrorCovariate{
private:
	int numContext;

	uint16_t _extractCovariate(const BAM::TSequencedBase & base) override {
		return genometools::index(base.context);
	};

public:
	static inline const std::string name =  "context";
	TSequencingErrorCovariate_context(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable);
	TSequencingErrorCovariate_context(const size_t FirstParameterIndex, const std::string & functionString);

	std::string typeString() const override { return name; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable) override;
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable) override;
	bool checkParameterRange(std::vector<uint16_t> & usedValues, uint16_t maxPos) override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable &) override {};
};

//-------------------------------------------
// TSequencingErrorCovariate_fragmentLength
//-------------------------------------------

class TSequencingErrorCovariate_fragmentLength:public TSequencingErrorCovariate{
private:
	uint16_t _extractCovariate(const BAM::TSequencedBase & base) override {
		return base.fragmentLength;
	};

public:
	static inline std::string name =  "fragmentLength";
	TSequencingErrorCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable);
	TSequencingErrorCovariate_fragmentLength(const size_t FirstParameterIndex, const std::string & functionString);

	std::string typeString() const override { return name; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable) override;
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable) override;
	bool checkParameterRange(std::vector<uint16_t> & usedValues, uint16_t maxPos) override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable) override;
};

//-------------------------------------------
// TSequencingErrorCovariate_mappingQuality
//-------------------------------------------

class TSequencingErrorCovariate_mappingQuality:public TSequencingErrorCovariate{
private:
	uint16_t _extractCovariate(const BAM::TSequencedBase & base) override {
		return base.mappingQuality;
	};

public:
	static inline const std::string name = "mappingQuality";
	TSequencingErrorCovariate_mappingQuality(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable);
	TSequencingErrorCovariate_mappingQuality(const size_t FirstParameterIndex, const std::string & functionString);

	std::string typeString() const override { return name; };
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString, const RecalEstimatorTools::TRecalDataTable & dataTable) override;
	void addFunction(const size_t FirstParameterIndex, const std::string & functionString) override;
	bool checkParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable) override;
	bool checkParameterRange(std::vector<uint16_t> & usedValues, uint16_t maxPos) override;
	void adjustParameterRange(const RecalEstimatorTools::TRecalDataTable & dataTable) override;
};

}; //end namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATE_H_ */
