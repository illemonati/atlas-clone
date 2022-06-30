/*
 * TSimulatorSoftClip.h
 *
 *  Created on: Jun 30, 2022
 *      Author: raphael
 */

#ifndef SIMULATIONS_TSIMULATORSOFTCLIP_H_
#define SIMULATIONS_TSIMULATORSOFTCLIP_H_

#include <string>
#include <vector>

#include "probability.h"

namespace Simulations {

// Pure abstract class
class TSimulatorSoftClipDist {
protected:
	TSimulatorSoftClipDist() = default;
public:
	virtual ~TSimulatorSoftClipDist() = default;
	virtual void sample(std::vector<uint16_t> &softClipLengths) const noexcept final {
		for (auto &q : softClipLengths) q = sample();
	}
	virtual uint16_t sample() const noexcept                                = 0;
};

//-------------------------------
// TSimulatorSoftClipDist
// Used for lengths of softclipped reads
//-------------------------------
// Class of a fixed value
class TSimulatorSoftClipDistFixed : public TSimulatorSoftClipDist {
private:
	uint16_t _max;
public:
	TSimulatorSoftClipDistFixed(std::string &s);

	uint16_t sample() const noexcept override { return _max; };
};

//------------------------------------------------
// TSimulatorSoftClipDistBinned
// Class of a binned distribution
//------------------------------------------------
class TSimulatorSoftClipDistBinned : public TSimulatorSoftClipDist {
private:
	std::vector<uint16_t> _sCBins;
public:
	TSimulatorSoftClipDistBinned(std::string &s);

	uint16_t sample() const noexcept override;
};

//------------------------------------------------
// TSimulatorSoftClipDistFreq
// Class of a binned distribution with frequencies
//------------------------------------------------
class TSimulatorSoftClipDistFreq : public TSimulatorSoftClipDist {
private:
	std::vector<uint16_t> _sCBins;
	std::vector<coretools::Probability> _frequencies;
	std::vector<coretools::Probability> _cumulativeFrequencies;

public:
	TSimulatorSoftClipDistFreq(std::string &s);

	uint16_t sample() const noexcept override;
};

//------------------------------------------------
// TSimulatorSoftClipDistPoisson
// Class of a Poisson distribution
//------------------------------------------------
class TSimulatorSoftClipDistPois : public TSimulatorSoftClipDist {
private:
	double _lambda;

public:
	TSimulatorSoftClipDistPois(std::string &s);

	uint16_t sample() const noexcept override;
};


} // namespace Simulations

#endif /* SIMULATIONS_TSIMULATORQUALITY_H_ */
