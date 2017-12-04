/*
 * TSimulatorQuality.h
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREAD_H_
#define TSIMULATORREAD_H_

#include "TGenotypeMap.h"
#include "TSimulatorReadLength.h"
#include "TPostMortemDamage.h"
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorQualityTransformation.h"


//-------------------------------
//TSimulatorRead
//-------------------------------
class TSimulatorRead{
protected:
	TRandomGenerator* randomGenerator;
	int maxPrintQual;
	std::string _name;

	TSimulatorReadLength* readLengthDist;
	bool readLengthInitialized;
	TSimulatorQualityDist* qualityDist;
	bool qualityDistInitialized;
	TSimulatorQualityTransformation* qualityTransform;
	bool qualityTransformInitialized;
	TPMD pmdObject;
	bool hasPMD;
	bool isInitialized;

	TGenotypeMap genoMap;
	TQualityMap qualityMap;

	//alignment
	BamTools::BamAlignment bamAlignment;
	int fragmentLength;
	Base* bases;
	int* qualities;


	//tmp variables
	int p;
	Base previousBase;
	int tmp_qual;

	//general functions
	void simulateQualitiesAndErrors(Base* _bases, int* _qualities, int & len);
	void applyPMD(Base* _bases, int & len, int & fragmentLength);

public:
	TSimulatorRead(std::string readGroupName, int MaxPrintQual, TRandomGenerator* RandomGenerator);
	~TSimulatorRead(){
		if(readLengthInitialized){
			delete readLengthDist;
			delete[] bases;
			delete[] qualities;
		}
		if(qualityDistInitialized)
			delete qualityDist;
		if(qualityTransformInitialized)
			delete qualityTransform;
	};

	bool checkInitialization();
	void setReadLengthDistribution(std::string s);
	void setQualityDistribution(std::string s);
	void setQualityTransformation(const std::string & type, const std::string & arg, TLog* logfile);
	void setPMD(const std::string & pmdStringCT, const std::string & pmdStringGA);

	std::string name(){ return _name; }
	double meanReadLength(){
		if(!readLengthInitialized) throw "Read length distribution not initialized!";
		return readLengthDist->mean();
	};
	double maxReadLength(){
		if(!readLengthInitialized) throw "Read length distribution not initialized!";
		return readLengthDist->max();
	};

	void setRefId(int refId){bamAlignment.RefID = refId; };
	void simulate(Base* haplotype, const long & pos, TSimulatorBamFile & bamFile);

	void printDetails(TLog* logfile);
};



#endif /* TSIMULATORREAD_H_ */
