/*
 * TDepthWriter.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TDEPTHWRITER_H_
#define GENOMETASKS_TDEPTHWRITER_H_

#include <string>

#include "coretools/Files/TOutputFile.h"

#include "TBamWindowTraverser.h"

namespace GenomeTasks{

//----------------------------------------
// TDepthWriter
//----------------------------------------
class TDepthWriter final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	coretools::TOutputFile _out;
	coretools::TCountDistribution<> _distPerSite;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome &) override {}
	void _endChromosome(const genometools::TChromosome&) override {}
public:
	void run();
};

}; // end namespace



#endif /* GENOMETASKS_TDEPTHWRITER_H_ */
