/*
 * TRecalibration.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#include "TRecalibration.h"

//---------------------------------------------------------------
//TRecalibration
//---------------------------------------------------------------

TRecalibration::TRecalibration(TReadGroupMap& ReadGroupMap):readGroupMapObject(ReadGroupMap){
	readGroupMapObject = ReadGroupMap;
}

int TRecalibration::findReadGroupIndex(std::string & name, BamTools::SamReadGroupDictionary & readGroups){
	int i = 0;
	for(std::vector<BamTools::SamReadGroup>::iterator it = readGroups.Begin(); it !=  readGroups.End(); ++it, ++i){
		if(name == it->ID) return i;
	}
	return -1;
}

/*
void TRecalibration::calcEmissionProbabilities(TSite & site){
	//first calculate for each base
	for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
		(*it)->fillEmissionProbabilitiesCore((**it).errorRate);
	}

	//then for the site
	site.calcEmissionProbabilities();
};
*/

double TRecalibration::getErrorRate(TBase & base){
	return base.errorRate;
}

/*
double TRecalibration::getErrorRateFromBase(const TBase & base){
	return base.errorRate;
//	return getErrorRate(base.readGroup, base.quality, base.posInRead, base.posInReadRev, base.context);
}
*/

int TRecalibration::getQuality(TBase & base){
	return qualityMap.errorToQuality(base.errorRate);
}

//int TRecalibration::getQualityFromBase(const TBase & base, TQualityMap & qualMap){
//	std::cout << "in base classe's getQualityFromBase" << std::endl;
//	return qualMap.errorToQuality(base.errorRate);
////	return getQuality(base.readGroup, base.quality, base.posInRead, base.posInReadRev, base.context);
//}
