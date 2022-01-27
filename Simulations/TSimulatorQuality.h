/*
 * TSimulatorQuality.h
 *
 *  Created on: Apr 27, 2021
 *      Author: phaentu
 */

#ifndef SIMULATIONS_TSIMULATORQUALITY_H_
#define SIMULATIONS_TSIMULATORQUALITY_H_

#include "../BAM/TSequencedBase.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TSequencingErrorModel.h"
#include "TSimulatorReadLength.h"
#include "algorithmsAndVectors.h"
#include "mathFunctions.h"

namespace Simulations {

using coretools::TRandomGenerator;
using genometools::PhredIntProbability;

// Pure abstract class
class TSimulatorQualityDist {
protected:
	TSimulatorQualityDist() = default;
public:
	virtual ~TSimulatorQualityDist() = default;
	virtual void sample(std::vector<PhredIntProbability> &phredInt) const noexcept final {
		for (auto &q : phredInt) q = sample();
	}
	virtual PhredIntProbability sample() const noexcept                                = 0;
	virtual void printDetails(coretools::TLog *logfile, const std::string &Name) const = 0;
};

//-------------------------------
// TSimulatorQualityDist
// Used for quality and mapping quality
//-------------------------------
// Class of a fixed value
class TSimulatorQualityDistFixed : public TSimulatorQualityDist {
private:
	PhredIntProbability _max{30};
public:
	TSimulatorQualityDistFixed(std::string &s);

	PhredIntProbability sample() const noexcept override { return _max; };
	void printDetails(coretools::TLog *logfile, const std::string &Name) const override;
};

//------------------------------------------------
// TSimulatorQualityDistBinned
// Class of a binned distribution
//------------------------------------------------
class TSimulatorQualityDistBinned : public TSimulatorQualityDist {
private:
	TRandomGenerator *_randomGenerator;
	std::vector<PhredIntProbability> _qualBins;
public:
	TSimulatorQualityDistBinned(std::string &s, TRandomGenerator *RandomGenerator);

	PhredIntProbability sample() const noexcept override;
	void printDetails(coretools::TLog *logfile, const std::string &Name) const override;
};

//------------------------------------------------
// TSimulatorQualityDistFreq
// Class of a binned distribution with frequencies
//------------------------------------------------
class TSimulatorQualityDistFreq : public TSimulatorQualityDist {
private:
	TRandomGenerator *_randomGenerator;
	std::vector<PhredIntProbability> _qualBins;
	std::vector<coretools::Probability> _frequencies;
	std::vector<coretools::Probability> _cumulativeFrequencies;

public:
	TSimulatorQualityDistFreq(std::string &s, TRandomGenerator *RandomGenerator);

	PhredIntProbability sample() const noexcept override;
	void printDetails(coretools::TLog *logfile, const std::string &Name) const override;
};

//------------------------------------------------
// TSimulatorQualityDistNormal
// Class of a normal distribution
//------------------------------------------------
class TSimulatorQualityDistNormal : public TSimulatorQualityDist {
private:
	TRandomGenerator *_randomGenerator;
	// densities
	std::vector<double> _densities;
	std::vector<double> _cumulDensities;
	double _mean = -1.;
	double _sd   = -1.;
	PhredIntProbability _min{0};
	PhredIntProbability _max{30};
public:
	TSimulatorQualityDistNormal(std::string &s, TRandomGenerator *RandomGenerator);
	TSimulatorQualityDistNormal(double mean, double sd, int min, int max, TRandomGenerator *RandomGenerator);

	void parseFunctionString(std::string &s);
	void fillDensities();
	double mean() const noexcept { return _mean; };
	double sd() const noexcept { return _sd; };

	PhredIntProbability sample() const noexcept override;
	void printDetails(coretools::TLog *logfile, const std::string &Name) const override;
};

} // namespace Simulations

#endif /* SIMULATIONS_TSIMULATORQUALITY_H_ */
