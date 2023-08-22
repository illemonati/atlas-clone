/*
 * TestimateErrors.h
 *
 */

#ifndef GENOMETASKS_TESTIMATEERRORS_H_
#define GENOMETASKS_TESTIMATEERRORS_H_

#include "TErrorEstimator.h"
#include "TGenome.h"

namespace GenomeTasks {

//-----------------------------------------------------------
// TEstimateRecalibration_base
//-----------------------------------------------------------
class TEstimateErrors final : public TGenome_windows {
private:
	GenotypeLikelihoods::TErrorEstimator _errors;
	bool _onlyLL = false;

	void _handleWindow() override;
	void _handleAlignment() override {}

public:
	TEstimateErrors();
	void run();
};

} // namespace GenomeTasks

#endif /* GENOMETASKS_TESTIMATERECALIBRATION_H_ */
