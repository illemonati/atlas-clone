//
// Created by reynac on 2/22/22.
//

#ifndef ATLAS_TF2ESTIMATOR_H
#define ATLAS_TF2ESTIMATOR_H

#include <string>
#include <vector>

#include "genometools/VCF/TPopulationLikelihoods.h"

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
	void calculateF2(const std::vector<size_t> &countsDiff);
	void writeF2(const std::vector<size_t> &countsDiff, const std::vector<double> &sampleF2,
	             const std::vector<double> &popF2);
};

} // namespace PopulationTools

#endif // ATLAS_TF2ESTIMATOR_H
