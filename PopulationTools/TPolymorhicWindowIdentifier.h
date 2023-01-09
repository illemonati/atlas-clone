/*
 * TPolymorhhicWindowIdentifier.h
 *
 *  Created on: Feb 14, 2020
 *      Author: wegmannd
 */

#ifndef POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_
#define POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_

#include <string>
#include "coretools/Main/TTask.h"


namespace PopulationTools{


class TPolymorhicWindowIdentifier{
public:
	void identifyPolymorphicWindows();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_identifyPolymorphicWindows:public coretools::TTask{
public:
	TTask_identifyPolymorphicWindows(){ _explanation = "Identifying windows for which samples are polymorphic"; };

	void run(){
		TPolymorhicWindowIdentifier identifier;
		identifier.identifyPolymorphicWindows();
	}
};

}; //end namespace

#endif /* POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_ */
