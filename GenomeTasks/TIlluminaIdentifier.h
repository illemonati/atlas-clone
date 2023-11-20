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

#include "coretools/Main/TTask.h"

#include "TBamFile.h"
#include "TGenome_OLD.h"
#include "TOutputBamFile.h"

namespace GenomeTasks{

class TIlluminaIdentifier : public old::TGenome_basic {
private:
	size_t _counter = 0;
	std::map<std::string, std::string> rgPU_rgID;
	BAM::TOutputBamFile _out;
	void _handleAlignment();

public:
    TIlluminaIdentifier();
	void run();
};

} //end namespace GenomeTasks

#endif /* GENOMETASKS_TILLUMINAIDENTIFIER_H_ */
