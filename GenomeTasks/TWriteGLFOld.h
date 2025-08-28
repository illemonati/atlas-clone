/*
 * TGLFWriter.h
 *
 *  Created on: Jun 6, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TWRITEGLFOLD_H_
#define GENOMETASKS_TWRITEGLFOLD_H_

#include "TBamWindowTraverser.h"
#include "genometools/GLF/TGLFWriter.h"

namespace GenomeTasks{

//-------------------------------------------
// TWriteGLFOld
//-------------------------------------------
class TWriteGLFOld final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	genometools::TGLFWriter _writer;
	bool _printAll;
	bool _curIsHapo = false;

	void _handleWindow(GenotypeLikelihoods::TWindow& Window) override;
	void _startChromosome(const genometools::TChromosome& Chr) override;
	void _endChromosome(const genometools::TChromosome&) override {}

public:
	TWriteGLFOld();
	void run();
};

}; //end namespace

#endif /* GENOMETASKS_TWRITEGLF_H_ */

