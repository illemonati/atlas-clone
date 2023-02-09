/*
 * TIlluminaIdentifier.h
 *
 *  Created on: Feb 9, 2023
 *      Author: raphael
 */

#ifndef GENOMETASKS_TILLUMINAIDENTIFIER_H_
#define GENOMETASKS_TILLUMINAIDENTIFIER_H_

#include <string>
#include <map>

#include "TGenome.h"
#include "coretools/Main/TTask.h"
#include "TBamFile.h"

namespace GenomeTasks{

class TIlluminaIdentifier:public TGenome_filtered{
private:
    std::map<std::string, std::string> rgPU_rgID;
	TOutputBamFile _out;
	void _handleAlignment();

public:
    TIlluminaIdentifier();
	void identify();
};



//-----------------------------------------
// Tasks
//-----------------------------------------

class TTask_identifyIllumina:public coretools::TTask{
public:
	TTask_identifyIllumina(){ _explanation = "Reassigning read groups based on the platform unit in their name"; };

	void run(){
		TIlluminaIdentifier identifyIllumina;
		identifyIllumina.identify();
	};
};

}; //end namespace GenomeTasks

#endif /* GENOMETASKS_TILLUMINAIDENTIFIER_H_ */