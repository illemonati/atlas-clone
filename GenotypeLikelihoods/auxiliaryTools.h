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
//#include <cstddef>
#include "TQualityMap.h"
#include "TGenotypeMap.h"
#include "stringFunctions.h"
#include "TLog.h"
#include "TReadGroups.h"
#include <algorithm>
#include "TFile.h"
#include <set>

namespace GenotypeLikelihoods{

//--------------------------------------------------------------
// TRecalibrationEMTransformationMap
//--------------------------------------------------------------
class TRecalibrationEMTransformationMap{
protected:
	uint16_t size;
	uint16_t max;
	std::vector<double> map;

public:
	TRecalibrationEMTransformationMap(){
		clear();
	};

	TRecalibrationEMTransformationMap(const uint16_t & Max){
		initialize(Max);
	};

	~TRecalibrationEMTransformationMap(){
		clear();
	};

	void clear(){
		max = 0;
		size = 0;
		map.clear();
	};

	void initialize(const uint16_t & Max){
		max = Max;
		size = Max + 1;
		map.resize(Max + 1);
		std::fill(map.begin(), map.end(), 0.0);
	};

	void set(const uint16_t & x, const double & value){
		map[x] = value;
	};

	bool checkRange(const uint16_t & val) const{
		if(val <= max) return true;
		else return false;
	};

	double operator[](const int & x){
		return map[x];
	};

	double get(const int & x) const{
		return map[x];
	};
};

class TRecalibrationEMQualityTransformationMap:public TRecalibrationEMTransformationMap{
public:
	TRecalibrationEMQualityTransformationMap(){
		initialize(BAM::PhredIntErrorRate::max().get());

		//now set
		for(BAM::PhredIntErrorRate q = BAM::PhredIntErrorRate::min(); q <= BAM::PhredIntErrorRate::max(); ++q){
			//convert phred int quality to error
			coretools::Probability eps(q);

			//ensure range
			if(eps < 0.000'000'000'1) eps = 0.000'000'000'1;
			else if(eps > 0.999'999'999'9) eps = 0.999'999'999'9;

			//then transform into logit space
			map[q.get()] = logit(eps);
		}
	}
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
