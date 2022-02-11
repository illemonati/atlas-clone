/*
 * TPolymorhhicWindowIdentifier.h
 *
 *  Created on: Feb 14, 2020
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_
#define POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_

#include "TPopulationLikelihoods.h"

namespace PopulationTools{


class TPolymorhicWindowIdentifier{
private:
	coretools::TLog* logfile;

public:
	TPolymorhicWindowIdentifier(coretools::TParameters & Parameters, coretools::TLog* logfile);
	void identifyPolymorphicWindows(coretools::TParameters & Parameters, coretools::TRandomGenerator* randomGenerator);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_identifyPolymorphicWindows:public coretools::TTask{
public:
	TTask_identifyPolymorphicWindows(){ _explanation = "Identifying windows for which samples are polymoprhic"; };

	void run(coretools::TParameters & Parameters, coretools::TLog* Logfile){
		TPolymorhicWindowIdentifier identifier(Parameters, Logfile);
		identifier.identifyPolymorphicWindows(Parameters, _randomGenerator);
	}
};

}; //end namespace

#endif /* POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_ */
