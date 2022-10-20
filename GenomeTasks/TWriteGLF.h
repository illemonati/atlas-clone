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
#include "coretools/Main/TTask.h"

namespace GenomeTasks{

//-------------------------------------------
// TWriteGLF
//-------------------------------------------
class TWriteGLF:public TGenome_windows{
private:
	GLF::TGlfWriter _writer;
	bool _printAll;
	void _handleWindow() override;
	void _handleAlignment() override {}

public:
	TWriteGLF();
	void writeGLF();
};

//-------------------------------------------
// Tasks
//-------------------------------------------
class TTask_writeGLF:public coretools::TTask{
public:
	TTask_writeGLF(){ _explanation = "Writing genotype likelihoods to a GLF file"; };

	void run(){
		TWriteGLF writer;
		writer.writeGLF();
	}
};

}; //end namespace

#endif /* GENOMETASKS_TWRITEGLF_H_ */

