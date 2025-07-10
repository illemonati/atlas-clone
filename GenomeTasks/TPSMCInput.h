/*
 * TPSMCInput.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPSMCINPUT_H_
#define GENOMETASKS_TPSMCINPUT_H_

#include <memory>

#include "TBamWindowTraverser.h"
#include "TThetaEstimator.h"
#include "coretools/Files/TLineWriter.h"

namespace GenomeTasks{

//----------------------------------------
// TPSMCInput
//----------------------------------------
class TPSMCInput final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	double _theta;
	double _confidence, _logConfidence, _logConfidenceHet;

	size_t _blockSize;
	coretools::TLineWriter _out;
	size_t _nCharOnLine;
	std::unique_ptr<GenotypeLikelihoods::TThetaEstimator> _thetaEstimator;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}
public:
	TPSMCInput();
	void run();
};

} // end namespace

#endif /* GENOMETASKS_TPSMCINPUT_H_ */
