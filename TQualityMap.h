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

	TQualityMap(){
		//only up to phred = 255, else always return 255
		minPhredInt = 255 - 33;
		sizeQual = minPhredInt + 34;
		phredIntToErrorMap = new double[minPhredInt+1];
		phredIntToLogErrorMap = new double[minPhredInt+1];
		qualityToErrorMap = new double[sizeQual];


		//initialize quality <= 0
		for(int i=0; i<33; ++i)
			qualityToErrorMap[i] = 1.0;
		phredIntToErrorMap[0] = 1.0;

		//and now others
		double tmp = - log(10) / 10.0;;
		for(int i=0; i<(minPhredInt+1); ++i){
			phredIntToErrorMap[i] = phredToError(i);
			phredIntToLogErrorMap[i] = i * tmp;
			qualityToErrorMap[i+33] = phredIntToErrorMap[i];
		}
		min = phredToError(minPhredInt);

		//Create map of illumina quality bins
		illuminaQualityBins = new int[sizeQual];
		for(int i=0; i<35; ++i)
			illuminaQualityBins[i] = 33;
		for(int i=35; i<43; ++i)
			illuminaQualityBins[i] = 39;
		for(int i=43; i<53; ++i)
			illuminaQualityBins[i] = 48;
		for(int i=53; i<58; ++i)
			illuminaQualityBins[i] = 55;
		for(int i=58; i<63; ++i)
			illuminaQualityBins[i] = 60;
		for(int i=63; i<68; ++i)
			illuminaQualityBins[i] = 66;
		for(int i=68; i<72; ++i)
			illuminaQualityBins[i] = 70;
		for(int i=72; i<sizeQual; ++i)
			illuminaQualityBins[i] = 73;
	};

	~TQualityMap(){
		delete[] phredIntToErrorMap;
		delete[] qualityToErrorMap;
		delete[] illuminaQualityBins;
		delete[] phredIntToLogErrorMap;
	};

	double phredIntToError(int phredInt){
		if(phredInt >= minPhredInt)
			return min;
		else return phredIntToErrorMap[phredInt];
	};

	inline double phredToError(double phred){
		return pow(10.0, -phred/10.0);
	};

	double qualityToError(int qual){
		if(qual < 33)
			return 1.0;
		return phredIntToError(qual - 33);
	};

	int phredIntToQuality(int phredInt){
		return phredInt + 33;
	};

	inline int errorToPhredInt(const double & errorRate){
		return round(errorToPhred(errorRate));
	};

	inline double errorToPhred(const double & errorRate){
		if(errorRate < min)
			return minPhredInt;
		else
			return -10.0 * log10(errorRate);
	};

	inline int errorToQuality(const double & errorRate){
		return round(-10.0 * log10(errorRate)) + 33;
	};

	double& operator[](int phred){
		//Note: no check on range!
		return phredIntToErrorMap[phred];
	};

	double& operator[](uint8_t & phred){
		//Note: no check on range!
		return phredIntToErrorMap[phred];
	};
};







#endif /* TQUALITYMAP_H_ */
