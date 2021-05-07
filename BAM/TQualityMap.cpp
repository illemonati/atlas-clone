/*
 * TQualityMap.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: phaentu
 */

#include "TQualityMap.h"

namespace BAM{

//TODO: change to strong types!

//---------------------------------------------------------------
//TQualityMap
//---------------------------------------------------------------
TQualityMap::TQualityMap(){
	//only up to phred = 255, else always return 255
	minPhredInt = 255;
	sizeQual = minPhredInt + 35;
	phredIntToErrorMap = new double[minPhredInt+1];
	phredIntToLogErrorMap = new double[minPhredInt+1];
	qualityToErrorMap = new double[sizeQual];
	qualityLimitSet = false;
	minQuality = 0;
	maxQuality = 255;

	//initialize quality <= 0
	for(uint8_t i=0; i<33; ++i)
		qualityToErrorMap[i] = 1.0;
	phredIntToErrorMap[0] = 1.0;

	//and now others
	double tmp = - log(10) / 10.0;
	for(uint32_t i=0; i<=minPhredInt; ++i){
		phredIntToErrorMap[i] = phredToError(i);
		phredIntToLogErrorMap[i] = i * tmp;
		qualityToErrorMap[i+33] = phredIntToErrorMap[i];
	}
	min = phredToError(minPhredInt);

	//Create map of illumina quality bins
	_binQualities = false; //whether or not to bin qualities when printing
	illuminaQualityBins = new uint8_t[sizeQual];
	for(uint16_t i=0; i<35; ++i)
		illuminaQualityBins[i] = 33;
	for(uint16_t i=35; i<43; ++i)
		illuminaQualityBins[i] = 39;
	for(uint16_t i=43; i<53; ++i)
		illuminaQualityBins[i] = 48;
	for(uint16_t i=53; i<58; ++i)
		illuminaQualityBins[i] = 55;
	for(uint16_t i=58; i<63; ++i)
		illuminaQualityBins[i] = 60;
	for(uint16_t i=63; i<68; ++i)
		illuminaQualityBins[i] = 66;
	for(uint16_t i=68; i<72; ++i)
		illuminaQualityBins[i] = 70;
	for(uint16_t i=72; i<sizeQual; ++i)
		illuminaQualityBins[i] = 73;
};

TQualityMap::~TQualityMap(){
	delete[] phredIntToErrorMap;
	delete[] qualityToErrorMap;
	delete[] illuminaQualityBins;
	delete[] phredIntToLogErrorMap;
};

void TQualityMap::setQualityLimits(uint8_t MinQualityForPrinting, uint8_t MaxQualityForPrinting){
	qualityLimitSet = true;
	minQuality = MinQualityForPrinting;
	maxQuality = MaxQualityForPrinting;
};

void TQualityMap::setBinQualityScores(bool BinQualityScores){
	_binQualities = BinQualityScores;
};

//to error
double TQualityMap::phredIntToError(uint8_t phredInt) const{
	if(phredInt >= minPhredInt)
		return min;
	else return phredIntToErrorMap[phredInt];
};

double TQualityMap::qualityToError(const int & qual) const{
	if(qual < 33)
		return 1.0;
	return phredIntToError(qual - 33);
};

//to illumina error
uint8_t TQualityMap::phredIntToIlluminaPhredInt(const uint8_t & phredInt) const{
	return illuminaQualityBins[phredInt + 33] - 33;
};

char TQualityMap::qualityToIlluminaQuality(const char & quality) const{
	return illuminaQualityBins[(uint8_t) quality];
};

//to quality
uint8_t TQualityMap::phredIntToQuality(const uint8_t & phredInt) const{
	uint16_t q = phredInt + 33;
	if(q > maxQuality)
		return maxQuality;
	else if(q < minQuality)
		return minQuality;
	else
		return q;
};

uint8_t TQualityMap::errorToQuality(const double & errorRate) const{
	uint16_t q = round(errorToPhred(errorRate)) + 33;
	if(q > maxQuality)
		return maxQuality;
	else if(q < minQuality)
		return minQuality;
	else
		return q;
};

//to phred and phred int
uint8_t TQualityMap::qualityToPhredInt(const uint8_t & quality) const{
	return quality - 33;
};

uint8_t TQualityMap::phredIntToPhredInt(const uint8_t & phredInt) const{
	return phredInt;
};

//function used when writing qualities
void TQualityMap::adjustQualitiesForWriting(std::string & qualities) const{
	if(_binQualities){
		for(auto& q : qualities){
			q = qualityToIlluminaQuality(q);
		}
	}
	if(qualityLimitSet){
		for(auto& q : qualities){
			if(q > maxQuality){
				q = maxQuality;
			} else if(q < minQuality){
				q = minQuality;
			}
		}
	}
};

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
TQualityFilter::TQualityFilter(TParameters & params, TLog* logfile){
	set(params, logfile);
};

bool TQualityFilter::set(TParameters & params, TLog* logfile){
	if(params.parameterExists("minQual") || params.parameterExists("maxQual")){
		int MinPhredInt = params.getParameterWithDefault<int>("minQual", 1);
		int MaxPhredInt = params.getParameterWithDefault<int>("maxQual", 93);

		if(MinPhredInt < 0 || MinPhredInt > 255) throw "minQual " + toString(MinPhredInt) + " is outside accepted range [0, 255]!";
		if(MaxPhredInt < 0 || MaxPhredInt > 255) throw "maxQual " + toString(MaxPhredInt) + " is outside accepted range [0, 255]!";
		if(MaxPhredInt < MinPhredInt) throw "maxQual must be >= minQual!";

		_set(MinPhredInt, MaxPhredInt);
		report(logfile);
		return true;
	} else {
		_default();
		return false;
	}
};

void TQualityFilter::report(TLog* logfile){
	logfile->list("Will filter out bases with quality outside the range [" + toString(_minPhredInt) + ", " + toString(_maxPhredInt) + "] (parameters 'minQual', 'maxQual')");
};

void TQualityFilter::_set(const uint8_t MinPhredInt, const uint8_t MaxPhredInt){
	_minPhredInt = MinPhredInt;
	_maxPhredInt = MaxPhredInt;
	_minQuality = _minPhredInt + 33;
	_maxQuality = _maxPhredInt + 33;
};

bool TQualityFilter::pass(uint8_t phredInt) const{
	if(phredInt < _minPhredInt || phredInt > _maxPhredInt)
		return false;
	return true;
};

bool TQualityFilter::pass(char quality) const{
	if(quality < _minQuality || quality > _maxQuality)
		return false;
	return true;
};

} // end namespace


