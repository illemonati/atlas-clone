/*
 * TSimulatorQuality.h
 *
 *  Created on: Apr 27, 2021
 *      Author: phaentu
 */

#ifndef SIMULATIONS_TSIMULATORQUALITY_H_
#define SIMULATIONS_TSIMULATORQUALITY_H_

#include <string>
#include <vector>

#include "PhredProbabilityTypes.h"
#include "probability.h"

namespace Simulations {

// Pure abstract class
class TSimulatorQualityDist {
protected:
	TSimulatorQualityDist() = default;
public:
	virtual ~TSimulatorQualityDist() = default;
	virtual void sample(std::vector<genometools::PhredIntProbability> &phredInt) const noexcept final {
		for (auto &q : phredInt) q = sample();
	}
	virtual genometools::PhredIntProbability sample() const noexcept                                = 0;
	virtual void printDetails(const std::string &Name) const = 0;
};

//-------------------------------
// TSimulatorQualityDist
// Used for quality and mapping quality
//-------------------------------
// Class of a fixed value
class TSimulatorQualityDistFixed : public TSimulatorQualityDist {
private:
	genometools::PhredIntProbability _max{30};
public:
	TSimulatorQualityDistFixed(std::string &s);

	genometools::PhredIntProbability sample() const noexcept override { return _max; };
	void printDetails(const std::string &Name) const override;
};

//------------------------------------------------
// TSimulatorQualityDistBinned
// Class of a binned distribution
//------------------------------------------------
class TSimulatorQualityDistBinned : public TSimulatorQualityDist {
private:
	std::vector<genometools::PhredIntProbability> _qualBins;
public:
	TSimulatorQualityDistBinned(std::string &s);

	genometools::PhredIntProbability sample() const noexcept override;
	void printDetails(const std::string &Name) const override;
};

//------------------------------------------------
// TSimulatorQualityDistFreq
// Class of a binned distribution with frequencies
//------------------------------------------------
class TSimulatorQualityDistFreq : public TSimulatorQualityDist {
private:
	std::vector<genometools::PhredIntProbability> _qualBins;
	std::vector<coretools::Probability> _frequencies;
	std::vector<coretools::Probability> _cumulativeFrequencies;

public:
	TSimulatorQualityDistFreq(std::string &s);

	genometools::PhredIntProbability sample() const noexcept override;
	void printDetails(const std::string &Name) const override;
};

//------------------------------------------------
// TSimulatorQualityDistNormal
// Class of a normal distribution
//------------------------------------------------
class TSimulatorQualityDistNormal : public TSimulatorQualityDist {
private:
	// densities
	std::vector<double> _densities;
	std::vector<double> _cumulDensities;
	double _mean = -1.;
	double _sd   = -1.;
	genometools::PhredIntProbability _min{0};
	genometools::PhredIntProbability _max{30};
public:
	TSimulatorQualityDistNormal(std::string &s);
	TSimulatorQualityDistNormal(double mean, double sd, int min, int max);

	void parseFunctionString(std::string &s);
	void fillDensities();
	double mean() const noexcept { return _mean; };
	double sd() const noexcept { return _sd; };

	genometools::PhredIntProbability sample() const noexcept override;
	void printDetails(const std::string &Name) const override;
};

} // namespace Simulations

#endif /* SIMULATIONS_TSIMULATORQUALITY_H_ */
