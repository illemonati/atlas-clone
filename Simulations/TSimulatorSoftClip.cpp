/*
 * TSimulatorSoftClip.cpp
 *
 *  Created on: Apr 27, 2021
 *      Author: phaentu
 */

#include "TSimulatorSoftClip.h"

#include <stddef.h>
#include <stdint.h>
#include <algorithm>
#include <memory>

#include "TLog.h"
#include "TRandomGenerator.h"
#include "algorithms.h"
#include "mathFunctions.h"
#include "stringFunctions.h"
#include "strongTypes.h"
#include "weakTypes.h"

namespace Simulations {
using coretools::instances::randomGenerator;
using coretools::str::toString;
using coretools::instances::logfile;

//----------------------------------
// TSimulatorSoftClipDist
//----------------------------------
TSimulatorSoftClipDistFixed::TSimulatorSoftClipDistFixed(std::string &s) {
	const auto pos1 = s.find("(");
	if (pos1 == std::string::npos)
		_max = coretools::str::convertStringCheck<uint8_t>(s);
	else if (pos1 == 0) {
		const auto pos2 = s.find(')');
		if (pos2 != s.size() - 1)
			throw "Failed to understand fixed distribution '" + s + "'!";
		_max = coretools::str::convertStringCheck<uint16_t>(s.substr(1, pos2 - 1));
	} else {
		throw "Failed to understand fixed distribution '" + s + "'!";
	}
}


//------------------------------------------------
// TSimulatorSoftClipDistBinned
//------------------------------------------------
TSimulatorSoftClipDistBinned::TSimulatorSoftClipDistBinned(std::string &s){
	const auto pos1 = s.find("(");
	if (pos1 == 0) {
		s.erase(0, 1);
		const auto pos2 = s.find(')');
		if (pos2 == std::string::npos || pos2 != s.size() - 1)
			throw "Failed to understand binned distribution '" + s + "'! Use binned(softClipLength_1,softClipLength_2,..,softClipLength_n).";
		s.erase(pos2, 1);
		coretools::str::fillContainerFromString(s, _sCBins, ',', true, true, false);
	} else {
		throw "Failed to understand binned distribution '" + s + "'! Use binned(softClipLength_1,softClipLength_2,..,softClipLength_n).";
	}
}

uint16_t TSimulatorSoftClipDistBinned::sample() const noexcept {
	return _sCBins[randomGenerator().sample(_sCBins.size())];
}


//------------------------------------------------
// TSimulatorSoftClipDistFreq
//------------------------------------------------
TSimulatorSoftClipDistFreq::TSimulatorSoftClipDistFreq(std::string &s){
	const auto pos1 = s.find("(");
	if (pos1 == 0) {
		s.erase(0, 1);
		const auto pos2 = s.find(')');
		if (pos2 == std::string::npos || pos2 != s.size() - 1) {
			throw "Failed to understand frequency distribution '" + s +
				"'! Use frequency(softClipLength_1:frequency_1,softClipLength_2:frequency_2, ... ,softClipLength_n:frequency_n).";
		}
		s.erase(pos2, 1);
		std::vector<std::string> tmp;
		coretools::str::fillContainerFromString(s, tmp, ',');

		// now parse each bin
		_sCBins.resize(tmp.size());
		_frequencies.resize(tmp.size());

		for (size_t i = 0; i < tmp.size(); ++i) {
			const auto pos3 = tmp[i].find(':');
			if (pos3 == std::string::npos) {
				throw "Failed to understand frequency distribution '" + s +
					"'! Use frequency(softClipLength_1:frequency_1,softClipLength_2:frequency_2, ... ,softClipLength_n:frequency_n).";
			}
			_sCBins[i]    = std::stoi(tmp[i].substr(0, pos3));
			_frequencies[i] = tmp[i].substr(pos3 + 1);
		}

		// fill cumulative
		coretools::fillCumulative(_frequencies, _cumulativeFrequencies);

	} else {
		throw "Failed to understand frequency distribution '" + s +
			"'! Use frequency(softClipLength_1:frequency_1,softClipLength_2:frequency_2, ... ,softClipLength_n:frequency_n).";
	}
}

uint16_t TSimulatorSoftClipDistFreq::sample() const noexcept {
	return _sCBins[randomGenerator().pickOne(_cumulativeFrequencies)];
}

//------------------------------------------------
// TSimulatorSoftClipDistPois
//------------------------------------------------
TSimulatorSoftClipDistPois::TSimulatorSoftClipDistPois(std::string &s) {
	const auto pos1 = s.find("(");
	if (pos1 == 0) {
		const auto pos2 = s.find(')');
		if (pos2 != s.size() - 1)
			throw "Failed to understand Poisson distribution '" + s + "'!";
		_lambda = coretools::str::convertStringCheck<double>(s.substr(1, pos2 - 1));
		} else {
			throw "Failed to understand Poisson distribution '" + s + "'!";
		}
}

uint16_t TSimulatorSoftClipDistPois::sample() const noexcept {
	return randomGenerator().getPoissonRandom(_lambda);
}

} // namespace Simulations
