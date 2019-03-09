/*
 * TRecalibration.h
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#ifndef TRecalIBRATION_H_
#define TRecalIBRATION_H_

#include "bamtools/api/BamReader.h"
#include "TSite.h"
#include <omp.h>
#include "TReadGroups.h"
#include "bamtools/api/SamHeader.h"
#include "TQualityMap.h"


//---------------------------------------------------------------
//TRecalibration: default = no recalibration
//---------------------------------------------------------------
class TRecalibration{
protected:
	TReadGroupMap& readGroupMapObject;
	TQualityMap qualityMap;

public:
	TRecalibration(TReadGroupMap& ReadGroupMap);

	virtual ~TRecalibration(){};

	virtual bool recalibrationChangesQualities(){
		return false;
	};

/*	char getQualityAsChar(const TBase & base, int & minOutQuality, int & maxOutQuality){
		int qual = getphredInt(base) + 33;
		if(qual > maxOutQuality) qual = maxOutQuality;
		if(qual < minOutQuality) qual = minOutQuality;
		return qual;
	};*/

	void calcEmissionProbabilities(TSite & site);
//	virtual double getErrorRate(const int & readGroupId, const int & quality, const int & pos, const int & posRev, const BaseContext & context);
	virtual double getErrorRate(TBase & base);
//	double getErrorRateFromBase(const TBase & base);
	virtual int getQuality(TBase & base);
//	virtual int getQualityFromBase(const TBase & base, TQualityMap & qualMap);
//	virtual int getphredInt(const TBase & base){
//		return base.phredInt;
//	};
//	virtual int getQuality(const TBase & base){
//		return base.quality;
//	};

	virtual bool requiresEstimation(){ return false;};
	int findReadGroupIndex(std::string & name, BamTools::SamReadGroupDictionary & readGroups);
};

#endif /* TRecalIBRATION_H_ */
