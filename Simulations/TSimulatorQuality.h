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
#include "TRandomGenerator.h"
#include "TQualityMap.h"
#include "TParameters.h"
#include "TSimulatorReadLength.h"
#include "TSequencingErrorModel.h"

namespace Simulations{

//-------------------------------
// TSimulatorQualityDist
// Used for quality and mapping quality
//-------------------------------
//Class of a fixed value
class TSimulatorQualityDist{
protected:
	TRandomGenerator* _randomGenerator;
	uint8_t _min;
	uint8_t _max;
	uint8_t _maxPlusOne;
	double _mean, _sd;

	virtual std::string _details() const;

public:
	TSimulatorQualityDist(TRandomGenerator* RandomGenerator);
	TSimulatorQualityDist(std::string & s, TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorQualityDist(){};
	uint8_t min(){ return _min; };
	uint8_t max(){ return _max; };
	double mean(){return _mean; };
	double sd(){return _sd; };

	virtual PhredIntProbability sample(){ return _max; };
	void sample(std::vector<PhredIntProbability> & phredInt) const;
	void printDetails(TLog* logfile, const std::string & Name) const;
};

//------------------------------------------------
// TSimulatorQualityDistBinned
// Class of a binned distribution
//------------------------------------------------
class TSimulatorQualityDistBinned:public TSimulatorQualityDist{
private:
	std::vector<uint16_t> _qualBins;

	std::string _details() const override;

public:
	TSimulatorQualityDistBinned(std::string & s, TRandomGenerator* RandomGenerator);
	PhredIntProbability sample() const override;
	void sample(std::vector<uint8_t> & qualities) const override;
};

//------------------------------------------------
// TSimulatorQualityDistFreq
// Class of a binned distribution with frequencies
//------------------------------------------------
class TSimulatorQualityDistFreq:public TSimulatorQualityDist{
private:
	std::vector<uint8_t> _qualBins;
	std::vector<double> _frequencies;
	std::vector<double> _cumulativeFrequencies;

	std::string _details() const override;

public:
	TSimulatorQualityDistFreq(std::string & s, TRandomGenerator* RandomGenerator);
	PhredIntProbability sample() const override;
};

//------------------------------------------------
// TSimulatorQualityDistNormal
// Class of a normal distribution
//------------------------------------------------
class TSimulatorQualityDistNormal:public TSimulatorQualityDist{
private:
	//densities
	int _size;
	std::vector<double> _densities;
	std::vector<double> _cumulDensities;

	std::string _details() const override;

public:
	TSimulatorQualityDistNormal(std::string & s, TRandomGenerator* RandomGenerator);
	TSimulatorQualityDistNormal(double & mean, double & sd, int min, int max, TRandomGenerator* RandomGenerator);;
	void parseFunctionString(std::string & s);
	void fillDensities();
	PhredIntProbability sample() const override;
};

}; //end namespace


#endif /* SIMULATIONS_TSIMULATORQUALITY_H_ */
