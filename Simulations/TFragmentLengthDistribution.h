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

#include "TCategoricalDistribution.h"
#include "TLog.h"

namespace Simulations {


//---------------------------------------------------------
// ReadLength
//---------------------------------------------------------
struct TReadAndFragmentLength {
	uint16_t read;
	uint16_t fragment;
	TReadAndFragmentLength(uint16_t Read, uint16_t Fragment) : read(Read), fragment(Fragment) {}
	uint16_t diff() { return fragment - read; } //danger if read > fragment!
};

//---------------------------------------------------------
// TReadLengthDistribution
//---------------------------------------------------------
class TFragmentLengthDistribution {
protected:
	uint16_t _maxReadLength{}; //number of cycles on an illumina machine
	coretools::probdist::TCategoricalDistribution<uint16_t> _fragmentLengthDistr;

public:
	TFragmentLengthDistribution(const std::string & FragmentLengthDistribution, const uint16_t MaxReadLength)
	: _maxReadLength(MaxReadLength),
	  _fragmentLengthDistr(FragmentLengthDistribution){};

	~TFragmentLengthDistribution() = default;

	[[nodiscard]] uint32_t max() const {
		return _fragmentLengthDistr.max();
	}

	[[nodiscard]] double mean() const {
		return _fragmentLengthDistr.mean();
	}

	[[nodiscard]] std::string functionString() const {
		return _fragmentLengthDistr.functionString();
	};

	void printDetails() const {
		//TODO: add max length info
		coretools::instances::logfile().list(_fragmentLengthDistr.verbalizedString());
	};

	[[nodiscard]] TReadAndFragmentLength sample() const {
		uint16_t fragmentLength = _fragmentLengthDistr.sample();
		if(fragmentLength > _maxReadLength){
			return TReadAndFragmentLength(_maxReadLength, fragmentLength);
		} else {
			return TReadAndFragmentLength(fragmentLength, fragmentLength);
		}
	}
};

}; // namespace Simulations

#endif /* TSIMULATORREADLENGTH_H_ */
