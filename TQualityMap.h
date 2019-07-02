/*
 * TQualityMap.h
 *
 *  Created on: Dec 7, 2018
 *      Author: phaentu
 */

#ifndef TQUALITYMAP_H_
#define TQUALITYMAP_H_

#include <math.h>
#include <cstdint>



//---------------------------------------------------------------
//TQualityMap
//---------------------------------------------------------------
class TQualityMap{
public:
	//IMPORTANT NOMENCLATURE
	//error is error rate between 0 and 1
	//phred is phred-scaled error as phred = -10 * log10(error)
	//phredInt is (int) phred
	//quality is phredInt + 33
	double* phredIntToErrorMap;
	double* qualityToErrorMap;
	double* phredIntToLogErrorMap;
	int* illuminaQualityBins;
	double min;
	int minPhredInt, sizeQual;

    TQualityMap();
    ~TQualityMap();

    double phredIntToError(int phredInt);
    double qualityToError(int qual);
    int phredIntToQuality(int phredInt);
    double& operator[](int phred);
    double& operator[](uint8_t & phred);

    inline double phredToError(double phred){
        return pow(10.0, -phred/10.0);
    }

    inline int errorToPhredInt(const double & errorRate){
        return round(errorToPhred(errorRate));
    }

    inline double errorToPhred(const double & errorRate){
        if(errorRate < min)
            return minPhredInt;
        else
            return -10.0 * log10(errorRate);
    }

    inline int errorToQuality(const double & errorRate){
        return round(-10.0 * log10(errorRate)) + 33;
    }

};







#endif /* TQUALITYMAP_H_ */
