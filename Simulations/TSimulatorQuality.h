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

//-------------------------------
// TSimulatorQualityDist
// Used for quality and mapping quality
//-------------------------------
// Class of a fixed value
class TSimulatorQualityDist {
protected:
	TRandomGenerator *_randomGenerator;
	PhredIntProbability _min{0};
	PhredIntProbability _max{30};
	double _mean = -1.;
	double _sd   = -1.;

	virtual std::string _details() const;
public:
	TSimulatorQualityDist(TRandomGenerator *RandomGenerator) : _randomGenerator(RandomGenerator){};
	TSimulatorQualityDist(std::string &s, TRandomGenerator *RandomGenerator);
	virtual ~TSimulatorQualityDist(){};
	PhredIntProbability min() const noexcept { return _min; };
	PhredIntProbability max() const noexcept { return _max; };
	double mean() const noexcept { return _mean; };
	double sd() const noexcept { return _sd; };

	virtual PhredIntProbability sample() const noexcept { return _max; };
	virtual void sample(std::vector<PhredIntProbability> &phredInt) const noexcept;
	void printDetails(coretools::TLog *logfile, const std::string &Name) const;
};

//------------------------------------------------
// TSimulatorQualityDistBinned
// Class of a binned distribution
//------------------------------------------------
class TSimulatorQualityDistBinned : public TSimulatorQualityDist {
private:
	std::vector<PhredIntProbability> _qualBins;
	std::string _details() const override;
public:
	TSimulatorQualityDistBinned(std::string &s, TRandomGenerator *RandomGenerator);
	PhredIntProbability sample() const noexcept override;
};

//------------------------------------------------
// TSimulatorQualityDistFreq
// Class of a binned distribution with frequencies
//------------------------------------------------
class TSimulatorQualityDistFreq : public TSimulatorQualityDist {
private:
	std::vector<PhredIntProbability> _qualBins;
	std::vector<coretools::Probability> _frequencies;
	std::vector<coretools::Probability> _cumulativeFrequencies;

	std::string _details() const override;
public:
	TSimulatorQualityDistFreq(std::string &s, TRandomGenerator *RandomGenerator);
	PhredIntProbability sample() const noexcept override;
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

	std::string _details() const override;
public:
	TSimulatorQualityDistNormal(std::string &s, TRandomGenerator *RandomGenerator);
	TSimulatorQualityDistNormal(double mean, double sd, int min, int max, TRandomGenerator *RandomGenerator);
	void parseFunctionString(std::string &s);
	void fillDensities();
	PhredIntProbability sample() const noexcept override;
};

} // namespace Simulations

#endif /* SIMULATIONS_TSIMULATORQUALITY_H_ */
