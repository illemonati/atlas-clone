/*
 * TGLFWriter.h
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TWRITEGLF_H_
#define GENOMETASKS_TWRITEGLF_H_

#include "TGLF.h"
#include "TBamWindowTraverser.h"

namespace GenomeTasks{

//-------------------------------------------
// TWriteGLF
//-------------------------------------------
class TWriteGLF final : public TBamWindowTraverser{
private:
	GLF::TGlfWriter _writer;
	bool _printAll;
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;

public:
	TWriteGLF();
	void run();
};

}; //end namespace

#endif /* GENOMETASKS_TWRITEGLF_H_ */

