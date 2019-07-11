/*
 * TSimulatorQuality.h
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREAD_H_
#define TSIMULATORREAD_H_

#include "TParameters.h"
#include "TGenotypeMap.h"
#include "TSimulatorReadLength.h"
#include "TPostMortemDamage.h"
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorQualityTransformation.h"


//-------------------------------
//TSimulatorSingleEndRead
//-------------------------------
class TSimulatorSingleEndRead{
protected:
	TRandomGenerator* randomGenerator;
	int maxPrintPhredInt;
	std::string _name;
	std::string readNamePrefix;
	int readXPos, readYPos;
	bool isInitialized;


	//read length
	TSimulatorReadLength* readLengthDist;
	bool readLengthInitialized;

	//qualities
	TSimulatorQualityDist* qualityDist;
	bool qualityDistInitialized;
	std::string qualDistType;
	TSimulatorQualityTransformation* qualityTransform;
	bool qualityTransformInitialized;

	//PMD
	TPMD pmdObject;
	bool hasPMD;

	//contamination
	bool isContaminated;
	double contaminationRate;
	TSimulatorReference* contaminationSource;

	TGenotypeMap genoMap;
	TQualityMap qualityMap;

	//alignment
	BamTools::BamAlignment bamAlignment;
	Base* bases;
	int* phredIntQualities;

	//general functions
	void simulateQualitiesAndErrors(Base* _bases, int* _qualities, int & len);
	void applyPMD(Base* _bases, BamTools::BamAlignment & alignment, int & fragmentLength);
	std::string getNextReadName();
	void fillAlignmentDetails(BamTools::BamAlignment & alignment, const Base* theBases, const int* thePhredIntQualities);

public:

	std::string type;

	TSimulatorSingleEndRead(std::string readGroupName, int readGroupNumber, int MaxPrintQual, TRandomGenerator* RandomGenerator);
	virtual ~TSimulatorSingleEndRead();

	bool checkInitialization();
	void setReadLengthDistribution(std::string s, TLog* logfile);
	void setQualityDistribution(std::string s);
	virtual void setQualityTransformation(TSimulatorQualityTransformParameters & parameters, TLog* logfile);
	void setPMD(const std::string & pmdStringCT, const std::string & pmdStringGA);
	void setContamination(double rate, TSimulatorReference* source);

	std::string name(){ return _name; };
	double meanReadLength(){
		if(!readLengthInitialized) throw "Read length distribution not initialized!";
		return readLengthDist->mean();
	};
	double maxReadLength(){
		if(!readLengthInitialized) throw "Read length distribution not initialized!";
		return readLengthDist->max();
	};

	void setRefId(int refId){bamAlignment.RefID = refId; };
	virtual void simulate(Base* haplotype, const long & pos, TSimulatorBamFile & bamFile);

	void printDetails(TLog* logfile);
	virtual void writeUnwrittenAlignments(const long & pos, TSimulatorBamFile & bamFile){};
};

//-------------------------------
//TSimulatorPairedEndReads
//-------------------------------
class TSimulatorPairedEndReads:public TSimulatorSingleEndRead{
private:
	std::vector<BamTools::BamAlignment*> bamAlignmentSecondMates;
	std::vector<BamTools::BamAlignment*> bamAlignmentSecondMates_idle;

	TSimulatorQualityTransformation* qualityTransform_secondMate;

	void initializeSecondMateAlignment(BamTools::BamAlignment & alignment);

public:
	TSimulatorPairedEndReads(std::string readGroupName, int readGroupNumber, int MaxPrintQual, TRandomGenerator* RandomGenerator);
	~TSimulatorPairedEndReads();

	void setQualityTransformation(TSimulatorQualityTransformParameters & parameters, TLog* logfile);
	void simulate(Base* haplotype, const long & pos, TSimulatorBamFile & bamFile);
	void writeUnwrittenAlignments(const long & pos, TSimulatorBamFile & bamFile);
};


#endif /* TSIMULATORREAD_H_ */
