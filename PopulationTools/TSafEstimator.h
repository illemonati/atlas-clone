#ifndef TSAFESTIMATOR_H_
#define TSAFESTIMATOR_H_

#include "TGlfMultiReader.h"

namespace PopulationTools {

class TSafEstimator {
	std::vector<double> _logProbs;
	GLF::TGlfMultiReader _glfReader;

	void _iterate(const TGenotypeLikelihoodsAllCombinationsVector& data, genometools::Base major);
public:
	TSafEstimator();
	void run();
};

} // namespace PopulationTools

#endif /* TMAJORMINOR_H_ */
