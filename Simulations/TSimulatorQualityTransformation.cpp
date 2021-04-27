/*
 * TSimulatorQualityTransformation.cpp
 *
 *  Created on: Nov 28, 2017
 *      Author: phaentu
 */


#include "TSimulatorQualityTransformation.h"

namespace Simulations{

//-----------------------------------------------
//TSimulatorQualityTransformation
//-----------------------------------------------
TSimulatorQualityTransformation::TSimulatorQualityTransformation(TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator){
	qualityDist = QualityDist;
	randomGenerator = RandomGenerator;
	p = 0;
};


void TSimulatorQualityTransformation::simulateQualitiesAndErrors(Base* bases, int* qualities, const int & len, const bool isReverseStrand){
	//simulate qualities
	qualityDist->sample(qualities, len);

	//add errors
	for(p=0; p<len; ++p){
		if(randomGenerator->getRand() < qualityMap.phredIntToErrorMap[qualities[p]])
			bases[p] = static_cast<Base>((bases[p] + randomGenerator->sample(3) + 1) % 4);
	}
};

void TSimulatorQualityTransformation::printDetails(TLog* logfile){
	logfile->list("Will write original qualities to BAM file (no transformation).");
};


//-------------------------------------
//TSimulatorQualityTransformationRecal
//-------------------------------------
TSimulatorQualityTransformationRecal::TSimulatorQualityTransformationRecal(std::string string, int maxReadLength, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator)
	:TSimulatorQualityTransformation(QualityDist, RandomGenerator){

	//code copied from TRecalbrationEM
	//read model tag
	size_t pos = string.find_first_of('[');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing '['!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";
	std::string modelTag = string.substr(0,pos);
	string.erase(0, pos+1);

	//read parameters: quality, position and context separted by semicolon (;)
	pos = string.find_first_of(']');
	if(pos == std::string::npos)
		throw "Failed to understand recal string: missing ']'!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";
	std::vector<std::string> tmpVec, vec;
	fillVectorFromString(string.substr(0, pos), tmpVec, ';');
	if(tmpVec.size() != 3)
		throw "Failed to understand recal string: wrong number of parameter sets (" + toString(tmpVec.size()) + " instead of 3)!\nEither provide a valid file name or a string of format 'modelTag[quality parameters; position parameters; context parameters]'.";

	//prec-alculate transformation table
	fillTransformationTable(modelTag, tmpVec, maxReadLength);
};

TSimulatorQualityTransformationRecal::TSimulatorQualityTransformationRecal(const std::string & modelTag, std::vector<std::string> & values, int maxReadLength, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator)
	:TSimulatorQualityTransformation(QualityDist, RandomGenerator){

	//pre-calculate transformation table
	fillTransformationTable(modelTag, values, maxReadLength);
};

void TSimulatorQualityTransformationRecal::fillTransformationTable(const std::string & modelTag, std::vector<std::string> & values, int maxReadLength){
	//set size
	_maxReadLengthPlusOne = maxReadLength + 1;
	maxQualPlusOne = qualityDist->max() + 1;
	numContext = genoMap.numContextsNotN;

	//allocate table
	transformedQuality = new int**[maxQualPlusOne];
	for(int q=0; q<maxQualPlusOne; ++q){
		transformedQuality[q] = new int*[_maxReadLengthPlusOne];
		for(int p=0; p<_maxReadLengthPlusOne; ++p)
			transformedQuality[q][p] = new int[numContext];
	}

	//fill table using a recal model
	//TODO!!!!
	//model = createTRecalibrationEMModel(modelTag, values, false, NULL);
	//model->fillTransformationTableForSimulation(transformedQuality, maxReadLengthPlusOne, maxQualPlusOne, qualityDist->min());
};

void TSimulatorQualityTransformationRecal::clearTransformationTable(){
	for(int q=0; q<maxQualPlusOne; ++q){
		for(p=0; p<_maxReadLengthPlusOne; ++p)
			delete[] transformedQuality[q][p];
		delete[] transformedQuality[q];
	}
	delete[] transformedQuality;
	//delete model;
};

void TSimulatorQualityTransformationRecal::printDetails(TLog* logfile){
	//TODO!!!!
	//logfile->list("Will transform qualities using the recal model " + model->getModelString());
};

void TSimulatorQualityTransformationRecal::simulateQualitiesAndErrors(Base* bases, int* qualities, const int & len, const bool isReverseStrand){
	//call base function to simulate
	TSimulatorQualityTransformation::simulateQualitiesAndErrors(bases, qualities, len, isReverseStrand);

	//transform quality
	Base previousBase = N;
	if(isReverseStrand){
		for(p=len - 1; p>=0; --p){
			qualities[p] = transformedQuality[qualities[p]][len - p - 1][genoMap.contextMapFlipped[previousBase][bases[p]]];
			previousBase = bases[p];
		}
	} else {
		for(p=0; p<len; ++p){
			qualities[p] = transformedQuality[qualities[p]][p][genoMap.contextMap[previousBase][bases[p]]];
			previousBase = bases[p];
		}
	}
};

}; //end namespace
