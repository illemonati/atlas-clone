/*
 * TRecalibrationEMAuxiliaryTools.h
 *
 *  Created on: Mar 11, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMAUXILIARYTOOLS_H_
#define TRECALIBRATIONEMAUXILIARYTOOLS_H_

#include <stddef.h>
#include <stdint.h>
#include <algorithm>
#include <vector>

#include "PhredProbabilityTypes.h"
#include "mathFunctions.h"
#include "probability.h"
#include "strongTypes.h"

namespace GenotypeLikelihoods{

struct TNoTransformation {
	constexpr double operator()(uint val) const noexcept {return static_cast<double>(val);}
	constexpr bool checkRange(uint) const noexcept {return true;}
};

class TQualityTransformation {
private:
	std::array<double, 256> _map;
public:	
	TQualityTransformation() {
		for(size_t v = 0; v <_map.size(); ++v) {
			_map[v] = logit(coretools::Probability(genometools::PhredIntProbability(v)));
		}
	}
	constexpr double operator()(uint16_t val) const noexcept {return _map[val];}
	constexpr bool checkRange(uint val) const noexcept {return val <= genometools::PhredIntProbability::max().get();}
};

//--------------------------------------------------------------------
// TRecalibrationEMFirstDerivative
//--------------------------------------------------------------------
struct TRecalibrationEMFirstDerivative{
	uint16_t index;
	double derivative;
};

struct TRecalibrationEMSecondDerivative{
	uint16_t index1;
	uint16_t index2;
	double derivative;
};

typedef std::vector<TRecalibrationEMFirstDerivative>::iterator TRecalibrationEMFirstDerivativesIterator;

class TRecalibrationEMFirstDerivatives{
private:
	std::vector<TRecalibrationEMFirstDerivative> _derivatives;
	TRecalibrationEMFirstDerivativesIterator _cur;

public:
	TRecalibrationEMFirstDerivatives(){};
	TRecalibrationEMFirstDerivatives(size_t Size);

	void resize(size_t Size);
	size_t size() const;
	void  restart();
	void add(uint16_t parameterIndex, double derivative);
	TRecalibrationEMFirstDerivativesIterator begin();
	TRecalibrationEMFirstDerivativesIterator end();
};

//--------------------------------------------------------------------
// TRecalibrationEMSecondDerivatives
//--------------------------------------------------------------------
typedef std::vector<TRecalibrationEMSecondDerivative>::iterator TRecalibrationEMSecondDerivativesIterator;

class TRecalibrationEMSecondDerivatives{
private:
	std::vector<TRecalibrationEMSecondDerivative> _derivatives;
	TRecalibrationEMSecondDerivativesIterator _cur;

public:
	TRecalibrationEMSecondDerivatives(){};
	TRecalibrationEMSecondDerivatives(size_t Size);

	void resize(size_t Size);
	size_t size() const;
	void  restart();
	void add(uint16_t parameterIndex1, uint16_t parameterIndex2, double derivative);
	TRecalibrationEMSecondDerivativesIterator begin();
	TRecalibrationEMSecondDerivativesIterator end();
};


}; //end namespace


#endif /* TRECALIBRATIONEMAUXILIARYTOOLS_H_ */
