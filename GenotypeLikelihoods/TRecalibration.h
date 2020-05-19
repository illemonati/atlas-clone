/*
 * TRecalibration.h
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#ifndef TRecalIBRATION_H_
#define TRecalIBRATION_H_

#include "../bamtools/api/BamReader.h"
#include "TSite.h"
#include <omp.h>
#include "../bamtools/api/SamHeader.h"

#include "TSequencingErrorModels.h"
#include "TQualityMap.h"

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
	virtual double getErrorRate(TBaseOld & base);
	virtual int getQuality(TBaseOld & base);

	virtual bool _requiresEstimation(){ return false;};
};


//--------------------------------------------------------------------
// TRecalibrationEM
//--------------------------------------------------------------------
class TRecalibrationEM:public TRecalibration{
private:
	TLog* logfile;
	GenotypeLikelihoods::TSequencingErrorModels models;
	TReadGroupMap* readGroupMap;

public:
	TRecalibrationEM(std::string string, TReadGroups* ReadGroups, TLog* Logfile);

	bool recalibrationChangesQualities(){ return true; };

	inline double getErrorRate(TBaseOld & base){
		return models.getErrorRate(base.data);
	};
	inline int getQuality(TBaseOld & base){
		double q = models.getErrorRate(base.data);
		return _qualityMap.errorToQuality(q);
	};
};

#endif /* TRecalIBRATION_H_ */
