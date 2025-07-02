/*
 * TGLFWriter.h
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TWRITEGLF_H_
#define GENOMETASKS_TWRITEGLF_H_

#include "TBamWindowTraverser.h"
#include "genometools/GLF/TGLFWriter.h"

namespace GenomeTasks{

//-------------------------------------------
// TWriteGLF
//-------------------------------------------
class TWriteGLF final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	genometools::TGLFWriter _writer;
	bool _printAll;
	bool _curIsHapo = false;

	void _handleWindow(GenotypeLikelihoods::TWindow& Window) override;
	void _startChromosome(const genometools::TChromosome& Chr) override;
	void _endChromosome(const genometools::TChromosome&) override {}

public:
	TWriteGLF();
	void run();
};

}; //end namespace

#endif /* GENOMETASKS_TWRITEGLF_H_ */

