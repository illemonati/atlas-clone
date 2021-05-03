/*
 * TRecalibrationEMAuxiliaryTools.h
 *
 *  Created on: Mar 11, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMAUXILIARYTOOLS_H_
#define TRECALIBRATIONEMAUXILIARYTOOLS_H_

#include <TSite.h>
#include <TGenotypeData.h>
#include <cstddef>
#include "TQualityMap.h"
#include "TGenotypeMap.h"
#include "stringFunctions.h"
#include "TLog.h"
#include "TReadGroups.h"
#include <algorithm>
#include "TFile.h"
#include <set>

namespace GenotypeLikelihoods{



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
	void add(const uint16_t & parameterIndex, const double & derivative);
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
	void add(const uint16_t & parameterIndex1, const uint16_t & parameterIndex2, const double & derivative);
	TRecalibrationEMSecondDerivativesIterator begin();
	TRecalibrationEMSecondDerivativesIterator end();
};


}; //end namespace


#endif /* TRECALIBRATIONEMAUXILIARYTOOLS_H_ */
