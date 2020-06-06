/*
 * TGLFWriter.h
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TWRITEGLF_H_
#define GENOMETASKS_TWRITEGLF_H_

#include "TGenome.h"
#include "TGLF.h"
#include "TTask.h"

namespace GenomeTasks{

//-------------------------------------------
// TWriteGLF
//-------------------------------------------
class TWriteGLF:public TGenome_windows{
private:
	TGlfWriter _writer;
	bool _printAll;
	void _handleWindow();

public:
	TWriteGLF(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void writeGLF();
};

//-------------------------------------------
// Tasks
//-------------------------------------------
class TTask_writeGLF:public TTask{
public:
	TTask_writeGLF(){ _explanation = "Writing genotype likelihoods to a GLF file"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TWriteGLF writer(Parameters, Logfile, _randomGenerator);
		writer.writeGLF();
	}
};

}; //end namespace

#endif /* GENOMETASKS_TWRITEGLF_H_ */

