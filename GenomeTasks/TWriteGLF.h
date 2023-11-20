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
#include "TGenome_OLD.h"
#include "coretools/Main/TTask.h"

namespace GenomeTasks{

//-------------------------------------------
// TWriteGLF
//-------------------------------------------
class TWriteGLF:public old::TGenome_windows{
private:
	GLF::TGlfWriter _writer;
	bool _printAll;
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _handleAlignment(BAM::TAlignment&) override {}

public:
	TWriteGLF();
	void run();
};

}; //end namespace

#endif /* GENOMETASKS_TWRITEGLF_H_ */

