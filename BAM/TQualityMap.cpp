/*
 * TQualityMap.cpp
 *
 *  Created on: Jul 31, 2020
 *      Author: phaentu
 */

#include "TQualityMap.h"

namespace BAM{

//TODO: change to strong types!
/*

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

*/

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
void TQualityFilter::_default(){
	//default values according to SAM specifications
	_filter = false;
	_range.set(genometools::PhredIntProbability(1), true, genometools::PhredIntProbability(93), true);
};

void TQualityFilter::set(coretools::TParameters & params, coretools::TLog* logfile){
	if(params.parameterExists("filterBaseQual")){
		params.fillParameter("filterBaseQual", _range);
		logfile->list("Will filter out bases with quality outside the range " + _range.rangeString() + " (parameter 'filterBaseQual')");
		_filter = true;
	} else {
		_default();
		logfile->list("Will keep bases regardless of quality (use 'filterBaseQual' to filter)");
	}
};

//-------------------------------------
// TContextFilter
//-------------------------------------
void TContextFilter::set(coretools::TParameters & params, coretools::TLog* logfile){
	_filter = false;
	if(params.parameterExists("ignoreContexts")){
		std::vector<std::string> contexts;
		params.fillParameterIntoContainer("ignoreContexts", contexts, ',');

		if(contexts.size() > 0){
			for(auto& c : contexts){
				if(c.size() != 2){
					throw "Context " + c + " does not consist of two bases! (parameter 'ignoreContexts')";
				}

				genometools::Base first(c[0]);
				genometools::Base second(c[1]);

				if((char) first != c[0] || (char) second != c[1]){
					throw "Unable to understand context '" + c + "'!  (parameter 'ignoreContexts')";
				}

				//save context
				genometools::BaseContext cc(first, second);
				_keptContexts[static_cast<uint8_t>(cc.get())] = false;
			}

			std::vector<std::string> rep;
			for(auto i = 0; i <= static_cast<uint8_t>(genometools::BaseContextEnum::cNN); ++i){
				if(!_keptContexts[i]){
					rep.push_back( (std::string) genometools::BaseContext(static_cast<genometools::BaseContextEnum>(i)));
				}
			}
			logfile->list("Will ignore the following contexts: " + coretools::str::concatenateString(rep, ", ")  + ". (parameter 'ignoreContexts')");
			_filter = true;
		}
	}

	if(!_filter){
		logfile->list("Will keep bases regardless of base context. (use 'ignoreContexts' to filter)");
	}
};

bool TContextFilter::pass(const TSequencedBase & base) const{
	return _keptContexts[static_cast<uint8_t>(base.context.get())];
};

} // end namespace


