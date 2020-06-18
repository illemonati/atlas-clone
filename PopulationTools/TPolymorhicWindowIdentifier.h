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
	TLog* logfile;
	TGlfConverter glfConverter;

public:
	TPolymorhicWindowIdentifier(TParameters & Parameters, TLog* logfile);
	void identifyPolymorphicWindows(TParameters & Parameters, TRandomGenerator* randomGenerator);
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_identifyPolymorphicWindows:public TTask{
public:
	TTask_identifyPolymorphicWindows(){ _explanation = "Identifying windows for which samples are polymoprhic"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TPolymorhicWindowIdentifier identifier(Parameters, Logfile);
		identifier.identifyPolymorphicWindows(Parameters, _randomGenerator);
	}
};

}; //end namespace

#endif /* POPULATIONTOOLS_TPOLYMORHICWINDOWIDENTIFIER_H_ */
