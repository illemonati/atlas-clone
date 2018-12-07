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
	int* illuminaQualityBins;
	double min;
	int minPhred, sizeQual;

	TQualityMap(){
		//only up to phred = 255, else always return 256
		minPhred = 256;
		sizeQual = minPhred + 33;
		phredIntToErrorMap = new double[minPhred];
		qualityToErrorMap = new double[sizeQual];

		//initialize quality <= 0
		for(int i=0; i<33; ++i)
			qualityToErrorMap[i] = 1.0;
		phredIntToErrorMap[0] = 1.0;

		//and now others
		for(int i=0; i<256; ++i){
			phredIntToErrorMap[i] = phredToError(i);
			qualityToErrorMap[i+33] = phredIntToErrorMap[i];
		}
		min = phredToError(minPhred);

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
	};

	double phredIntToError(int phredInt){
		if(phredInt >= minPhred)
			return min;
		else return phredIntToErrorMap[phredInt];
	};

	inline double phredToError(double phred){
		return pow(10.0, -phred/10.0);
	};

	double qualityToError(int qual){
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
			return minPhred;
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
