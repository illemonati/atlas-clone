/*
 * TSimulatorReadLength.h
 *
 *  Created on: Oct 6, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREADLENGTH_H_
#define TSIMULATORREADLENGTH_H_

#include <cstdint>
#include <string>
#include <vector>

namespace Simulations {


//---------------------------------------------------------
// ReadLength
//---------------------------------------------------------
struct TReadLength {
	uint16_t read;
	uint16_t fragment;
	TReadLength(uint16_t Read, uint16_t Fragment) : read(Read), fragment(Fragment) {}
	uint16_t diff() { return fragment - read; } //danger if read > fragment!
};

//---------------------------------------------------------
// TReadLengthDistribution
//---------------------------------------------------------
class TReadLengthDistribution {
protected:
	uint32_t _meanLength = -1;
	std::vector<double> _positionProbs; // normalized (1 - cumulDensity)

public:
	TReadLengthDistribution(std::string &s);
	TReadLengthDistribution() = default;
	virtual ~TReadLengthDistribution() = default;

	double operator[](const uint32_t position) const { return _positionProbs[position]; };

	virtual TReadLength sample() const noexcept { return TReadLength(_meanLength, _meanLength); }
	virtual uint32_t max() const noexcept { return _meanLength; }
	virtual double mean() const noexcept { return _meanLength; }
	virtual double probAcceptance() const noexcept { return 1.0; }
	virtual void printDetails();
};

class TSimulatorReadLengthGamma : public TReadLengthDistribution {
protected:
	double _meanLength = -1.;
	double _alpha = -1.;
	double _beta = -1.;
	uint32_t _min = 0;
	uint32_t _maxPlusOne = 1;

	std::vector<double> _gammaDensity;
	std::vector<double> _gammaCumulDensity;

	void parseFunctionString(std::string &s, double &param1, double &param2);
	void initiate();

public:
	TSimulatorReadLengthGamma(std::string &s);
	TSimulatorReadLengthGamma() = default;
	TReadLength sample() const noexcept override;
	uint32_t max() const noexcept override { return _maxPlusOne - 1; };
	double mean() const noexcept override { return _meanLength; };
	void printDetails() override;
	double probAcceptance() const noexcept override{ return 1.0 - _gammaCumulDensity.front(); };
};

class TSimulatorReadLengthGammaMode : public TSimulatorReadLengthGamma {
protected:
	double _mode, _var;

public:
	TSimulatorReadLengthGammaMode(std::string &s);
	void printDetails() override;
};

}; // namespace Simulations

#endif /* TSIMULATORREADLENGTH_H_ */
