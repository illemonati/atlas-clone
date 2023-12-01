/*
 * TPSMCInput.h
 *
 *  Created on: Jun 5, 2020
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TPSMCINPUT_H_
#define GENOMETASKS_TPSMCINPUT_H_

#include <iosfwd>
#include <memory>

#include "TBamWindowTraverser.h"
#include "TGenotypeData.h"
#include "TThetaEstimator.h"

namespace GenomeTasks{

//----------------------------------------
// TPSMCInput
//----------------------------------------
class TPSMCInput final : public TBamWindowTraverser{
private:
	double _theta;
	double _confidence, _logConfidence, _logConfidenceHet;

	size_t _blockSize;
	std::ofstream _out;
	size_t _nCharOnLine;
	std::unique_ptr<GenotypeLikelihoods::TThetaEstimator> _thetaEstimator;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _onChrChange(const genometools::TChromosome&) override {}
public:
	TPSMCInput();
	void run();
};

} // end namespace

#endif /* GENOMETASKS_TPSMCINPUT_H_ */
