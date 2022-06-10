//
// Created by reynac on 2/22/22.
//

#ifndef ATLAS_TF2ESTIMATOR_H
#define ATLAS_TF2ESTIMATOR_H

#include <cstdint>
#include <string>
#include <vector>

#include "TParameters.h"
#include "TPopulation.h"
#include "TPopulationLikelihoods.h"
#include "TTask.h"

namespace PopulationTools {
//------------------------------------------------
// TF2Estimator
//------------------------------------------------

class TF2Estimator {
private:
	std::string _outname;

	// vcf reader
	genometools::TPopulationLikelihoodReaderLocus _reader;

	// specify population and samples to keep
	genometools::TPopulationSamples _samples;

	void _openVCF();
	void _matchSamples();
	void _readVCF();

public:
	void run();
	void calculateF2(const std::vector<uint64_t> &countsDiff);
	void writeF2(const std::vector<uint64_t> &countsDiff, const std::vector<double> &sampleF2,
	             const std::vector<double> &popF2);
};

//--------------------------------------
// Tasks
//--------------------------------------

class TTask_calculateF2 : public coretools::TTask {
public:
	TTask_calculateF2() {
		_explanation = "Calculate F2 between different samples, and within and between populations";
	};

	void run() override {
		TF2Estimator F2;
		F2.run();
	};
};
} // namespace PopulationTools

#endif // ATLAS_TF2ESTIMATOR_H
