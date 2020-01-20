/*
 * TRecalibration.h
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#ifndef TRecalIBRATION_H_
#define TRecalIBRATION_H_

#include "../bamtools/api/BamReader.h"
#include "../TSite.h"
#include <omp.h>
#include <TRecalibrationEMModel.h>
#include "../bamtools/api/SamHeader.h"
#include "../TQualityMap.h"

//---------------------------------------------------------------
//TRecalibration: default = no recalibration
//---------------------------------------------------------------
class TRecalibration{
protected:
	TQualityMap _qualityMap;
	std::string _type;

public:
	TRecalibration();

	virtual ~TRecalibration(){};

	virtual bool recalibrationChangesQualities(){
		return false;
	};

	std::string type(){ return _type; };

	void calcEmissionProbabilities(TSite & site);
	virtual double getErrorRate(TBase & base);
	virtual int getQuality(TBase & base);

	virtual bool _requiresEstimation(){ return false;};
};


//--------------------------------------------------------------------
// TRecalibrationEM
//--------------------------------------------------------------------
class TRecalibrationEM:public TRecalibration{
private:
	TLog* logfile;
	TReadGroups* readGroups;
	TRecalibrationEMModels* models;

	void _initializeRecalibrationParametersFromString(std::string & string);
	void _initializeRecalibrationParametersFromFile(std::string filename);

public:
	TRecalibrationEM(std::string string, TReadGroups* ReadGroups, TLog* Logfile);
	~TRecalibrationEM(){
		delete models;
	};

	bool recalibrationChangesQualities(){ return true; };
	void initializeRecalibrationParameters(std::string string);

	inline double getErrorRate(TBase & base){
		return models->getErrorRate(base);
	};
	inline int getQuality(TBase & base){
		double q = models->getErrorRate(base);
		return _qualityMap.errorToQuality(q);
	};
};

#endif /* TRecalIBRATION_H_ */
