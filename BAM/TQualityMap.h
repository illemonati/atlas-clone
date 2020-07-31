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

#include "TParameters.h"

//TODO: make things private

namespace BAM{

//---------------------------------------------------------------
//TQualityMap
//---------------------------------------------------------------
class TQualityMap{
private:
	bool _binQualities;

public:
	//IMPORTANT NOMENCLATURE
	//error is error rate between 0 and 1
	//phred is phred-scaled error as phred = -10 * log10(error)
	//phredInt is (int) phred
	//quality is phredInt + 33
	double* phredIntToErrorMap;
	double* qualityToErrorMap;
	double* phredIntToLogErrorMap;
	uint8_t* illuminaQualityBins;
	double min;
	uint8_t minPhredInt;
	uint16_t sizeQual;
	bool  qualityLimitSet;
	uint8_t minQuality, maxQuality;

	TQualityMap();
	~TQualityMap();

	void setQualityLimits(uint8_t MinQualityForPrinting, uint8_t MaxQualityForPrinting);
	void  setBinQualityScores(bool BinQualityScores);
	double phredIntToError(uint8_t phredInt) const;

	inline double phredToError(double phred) const{
		return pow(10.0, -phred/10.0);
	};

	inline uint8_t errorToPhredInt(const double & errorRate) const{
		uint16_t phred = round(errorToPhred(errorRate));
		if(phred > minPhredInt)
			return minPhredInt;
		else
			return phred;
	};

	inline double errorToPhred(const double & errorRate) const{
		return -10.0 * log10(errorRate);
	};

	double qualityToError(int qual) const;
	uint8_t phredIntToIlluminaPhredInt(const uint8_t phredInt) const;
	char qualityToIlluminaQuality(char quality) const;
	uint8_t phredIntToQuality(uint8_t phredInt) const;
	uint8_t errorToQuality(const double & errorRate) const;
	uint8_t qualityToPhredInt(uint8_t quality) const;
	void adjustQualitiesForWriting(std::string & qualities) const;
};

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
class TQualityFilter{
private:
	uint8_t _minPhredInt, _maxPhredInt;
	char _minQuality, _maxQuality;

	void _default(){ _set(1,93); };

public:
	TQualityFilter(){ _default(); };

	TQualityFilter(TParameters & params, TLog* logfile);
	bool set(TParameters & params, TLog* logfile);
	void report(TLog* logfile);
	void _set(const uint8_t MinPhredInt, const uint8_t MaxPhredInt);
	uint8_t minPhredInt() const{ return _minPhredInt; };
	uint8_t maxPhredInt() const{ return _maxPhredInt; };
	bool pass(uint8_t phredInt) const;
	bool pass(char quality) const;
};

}; //end namespace



#endif /* TQUALITYMAP_H_ */
