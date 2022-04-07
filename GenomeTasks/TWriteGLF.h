/*
 * TGLFWriter.h
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TWRITEGLF_H_
#define GENOMETASKS_TWRITEGLF_H_

#include <string>

#include "TGLF.h"
#include "TGenome.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TTask.h"

namespace GenomeTasks{

//-------------------------------------------
// TWriteGLF
//-------------------------------------------
class TWriteGLF:public TGenome_windows{
private:
	GLF::TGlfWriter _writer;
	bool _printAll;
	void _handleWindow();

public:
	TWriteGLF(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void writeGLF();
};

//-------------------------------------------
// Tasks
//-------------------------------------------
class TTask_writeGLF:public coretools::TTask{
public:
	TTask_writeGLF(){ _explanation = "Writing genotype likelihoods to a GLF file"; };

	void run(){
		using namespace coretools::instances;
		TWriteGLF writer(parameters(), &logfile(), &randomGenerator());
		writer.writeGLF();
	}
};

}; //end namespace

#endif /* GENOMETASKS_TWRITEGLF_H_ */

