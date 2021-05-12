/*
 * TSimulatorQualityTransformation.h
 *
 *  Created on: Nov 28, 2017
 *      Author: phaentu
 */

#ifndef TSIMULATORQUALITYTRANSFORMATION_H_
#define TSIMULATORQUALITYTRANSFORMATION_H_

#include "../BAM/TSequencedBase.h"
#include "TLog.h"
#include "TRandomGenerator.h"
#include "TQualityMap.h"
#include "TParameters.h"
#include "TSimulatorReadLength.h"
#include "TSequencingErrorModel.h"

namespace Simulations{

//-----------------------------------------------
//TSimulatorQualityTransformation
//-----------------------------------------------
class TSimulatorQualityTransformation{
protected:
	TSimulatorQualityDist* qualityDist;
	TRandomGenerator* randomGenerator;
	BAM::TQualityMap qualityMap;

	//tmp vars
	int p;

public:
	TSimulatorQualityTransformation(TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorQualityTransformation(){};
	virtual void simulateQualitiesAndErrors(Base* bases, int* qualities, const int & len, const bool isReverseStrand);
	virtual void printDetails(TLog* logfile);
};


//-------------------------------------
//TSimulatorQualityTransformationRecal
//-------------------------------------
class TSimulatorQualityTransformationRecal:public TSimulatorQualityTransformation{
private:
	//TRecalibrationEMModel model;
	int _maxReadLengthPlusOne;
	int maxQualPlusOne;
	int numContext;
	int*** transformedQuality; //index are [qual][pos][context]
	GenotypeLikelihoods::TGenotypeMap genoMap;

	//private functions
	void fillTransformationTable(const std::string & modelTag, std::vector<std::string> & values, int maxReadLength);
	void clearTransformationTable();
	void transformQualities(Base* bases, int* qualities, const int & len, bool & isReverseStrand);

public:
	TSimulatorQualityTransformationRecal(std::string string, int maxReadLength, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator);
	TSimulatorQualityTransformationRecal(const std::string & modelTag, std::vector<std::string> & values, int maxReadLength, TSimulatorQualityDist* QualityDist, TRandomGenerator* RandomGenerator);
	~TSimulatorQualityTransformationRecal(){
		clearTransformationTable();
	};

	int getTransformedQuality(int qual, int pos, int context){
		return(transformedQuality[qual][pos][context]);
	};

	void simulateQualitiesAndErrors(Base* bases, int* qualities, const int & len, const bool isReverseStrand);
	void printDetails(TLog* logfile);
};

//---------------------------------------------------------
//TSimulatorQualityTransformParameters
//---------------------------------------------------------
struct TSimulatorQualityTransformParameters{
	std::string type; //either none, recal or BQSR
	std::string parameters_firstMate;
	std::string parameters_secondMate;

	TSimulatorQualityTransformParameters(){};
	TSimulatorQualityTransformParameters(std::string Type, std::string Parameters_firstMate, std::string Parameters_secondMate){
		type = Type;
		parameters_firstMate = Parameters_firstMate;
		parameters_secondMate = Parameters_secondMate;
	};

};

}; //end namespace

#endif /* TSIMULATORQUALITYTRANSFORMATION_H_ */
