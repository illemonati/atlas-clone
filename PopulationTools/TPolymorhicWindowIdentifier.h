/*
 * TPolymorhhicWindowIdentifier.h
 *
 *  Created on: Feb 14, 2020
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_
#define POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_

#include <string>
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TTask.h"


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
	TTask_identifyPolymorphicWindows(){ _explanation = "Identifying windows for which samples are polymorphic"; };

	void run(){
		using namespace coretools::instances;
		TPolymorhicWindowIdentifier identifier(parameters(), &logfile());
		identifier.identifyPolymorphicWindows(parameters(), &randomGenerator());
	}
};

}; //end namespace

#endif /* POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_ */
