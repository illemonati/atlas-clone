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
#include "TSequencedBase.h"

//TODO: make things private

namespace BAM{

//---------------------------------------------------------------
//TQualityMap
//---------------------------------------------------------------
/*
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

	void setQualityLimits(const uint8_t MinQualityForPrinting, const uint8_t MaxQualityForPrinting);
	void  setBinQualityScores(const bool BinQualityScores);
	double phredIntToError(const uint8_t phredInt) const;

	inline double phredToError(const double & phred) const{
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

	double qualityToError(const int & qual) const;
	uint8_t phredIntToIlluminaPhredInt(const uint8_t & phredInt) const;
	char qualityToIlluminaQuality(const char & quality) const;
	uint8_t phredIntToQuality(const uint8_t & phredInt) const;
	uint8_t phredIntToPhredInt(const uint8_t & phredInt) const;
	uint8_t errorToQuality(const double & errorRate) const;
	uint8_t qualityToPhredInt(const uint8_t & quality) const;
	void adjustQualitiesForWriting(std::string & qualities) const;
};
*/

//-------------------------------------
// TBaseFilter
//-------------------------------------
class TBaseFilter{
protected:
	bool _filter;

	constexpr operator bool() const{
		return _filter;
	};

public:
	explicit constexpr TBaseFilter() : _filter(false) {};
	~TBaseFilter() = default;

	virtual bool pass(const TSequencedBase & base) const = 0;
};

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
class TQualityFilter : public TBaseFilter{
private:
	PhredIntErrorRate _minPhredInt, _maxPhredInt;
	void _default();

public:
	explicit constexpr TQualityFilter(){
		_default();
	};

	TQualityFilter(TParameters & params, TLog* logfile){
		set(params, logfile);
	};

	~TQualityFilter() = default;

	void set(TParameters & params, TLog* logfile);

	bool pass(const TSequencedBase & base) const;
};

//-------------------------------------
// TContextFilter
//-------------------------------------
class TContextFilter : public TBaseFilter{
private:
	std::array<bool, static_cast<uint8_t>(BaseContextEnum::cNN) + 1> _keptContexts;

public:
	explicit constexpr TContextFilter(){
		_keptContexts.fill(true);
	};
	~TContextFilter() = default;

	void set(TParameters & params, TLog* logfile);

	bool pass(const TSequencedBase & base) const;
};


}; //end namespace BAM



#endif /* TQUALITYMAP_H_ */
