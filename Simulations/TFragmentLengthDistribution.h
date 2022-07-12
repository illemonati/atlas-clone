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
	std::unique_ptr<coretools::probdist::TCategoricalDistribution<uint16_t>> _fragmentLengthDistr;
	uint16_t _maxReadLength{}; //number of cycles on an illumina machine

public:
	TFragmentLengthDistribution(const uint16_t MaxReadLength, std::string &s){
		set(MaxReadLength, s);
	}
	TFragmentLengthDistribution() = default;

	virtual ~TFragmentLengthDistribution() = default;

	void set(const uint16_t MaxReadLength, std::string &s){
		_maxReadLength = MaxReadLength;
		coretools::probdist::createDiscreteDistribution(_fragmentLengthDistr, s);
	}

	operator bool() const {
		return (bool) _fragmentLengthDistr;
	}

	[[nodiscard]] uint32_t max() const noexcept {
		return _fragmentLengthDistr->max();
	}

	[[nodiscard]] double mean() const noexcept {
		return _fragmentLengthDistr->mean();
	}

	[[nodiscard]] std::string functionString(){
		return _fragmentLengthDistr->functionString();
	};

	void printDetails(){
		//TODO: add max length info
		coretools::instances::logfile().list(_fragmentLengthDistr->verbalizedString());
	};

	[[nodiscard]] TReadAndFragmentLength sample() const noexcept {
		uint16_t fragmentLength = _fragmentLengthDistr->sample();
		if(fragmentLength > _maxReadLength){
			return TReadAndFragmentLength(_maxReadLength, fragmentLength);
		} else {
			return TReadAndFragmentLength(fragmentLength, fragmentLength);
		}
	}
};

}; // namespace Simulations

#endif /* TSIMULATORREADLENGTH_H_ */
