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

namespace Simulations{

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
	GenotypeLikelihoods::TPMDDoubleStrand pmdObject;
	bool hasPMD;

	//contamination
	bool isContaminated;
	double contaminationRate;
	TSimulatorReference* contaminationSource;

	TGenotypeMap genoMap;
	TQualityMap qualityMap;

	//alignment
	BAM::TCigar _cigar;
	BAM::TSamFlags _flags;
	BAM::TAlignment _alignment;
	Base* bases;
	int* phredIntQualities;

	//general functions
	void _simulateQualitiesAndErrors(Base* _bases, int* _qualities, int & len);
	void _applyPMD(Base* _bases, const uint16_t readLength, const uint16_t fragmentLength, const bool isReverseStrand);
	std::string _getNextReadName();
	void _fillAlignmentDetails(BAM::TAlignment & alignment, const uint16_t Length, const Base* theBases, const int* thePhredIntQualities);

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

	virtual void simulate(Base* haplotype, const uint32_t refID, const long & pos, TSimulatorBamFile & bamFile);

	void printDetails(TLog* logfile);
	virtual void writeUnwrittenAlignments(const long & pos, TSimulatorBamFile & bamFile){};
};

//-------------------------------
//TSimulatorPairedEndReads
//-------------------------------
class TSimulatorPairedEndReads:public TSimulatorSingleEndRead{
private:
	BAM::TSamFlags _mateFlags;
	//BAM::TAlignment _mate;

	std::vector<BAM::TAlignment> bamAlignmentSecondMates;
	std::vector<BAM::TAlignment> bamAlignmentSecondMates_idle;

	TSimulatorQualityTransformation* qualityTransform_secondMate;

	void initializeSecondMateAlignment(BamTools::BamAlignment & alignment);

public:
	TSimulatorPairedEndReads(std::string readGroupName, int readGroupNumber, int MaxPrintQual, TRandomGenerator* RandomGenerator);
	~TSimulatorPairedEndReads();

	void setQualityTransformation(TSimulatorQualityTransformParameters & parameters, TLog* logfile);
	void simulate(Base* haplotype, const uint32_t refID, const long & pos, TSimulatorBamFile & bamFile);
	void writeUnwrittenAlignments(const long & pos, TSimulatorBamFile & bamFile);
};

}; //end namespace

#endif /* TSIMULATORREAD_H_ */
