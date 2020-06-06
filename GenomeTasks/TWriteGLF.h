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

namespace GenomeTasks{

//-------------------------------------------
// TWriter
//-------------------------------------------
class TWriteGLF:public TGenome_windows{
private:
	TGlfWriter _writer;
	bool _printAll;
	void _handleWindow();

public:
	TWriteGLF(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void writeGLF();
};


}; //end namespace

#endif /* GENOMETASKS_TWRITEGLF_H_ */

