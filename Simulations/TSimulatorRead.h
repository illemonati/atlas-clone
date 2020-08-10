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
	TRandomGenerator* _randomGenerator;
	TGenotypeMap& _genoMap;

	std::string _type;
	int _maxPrintPhredInt;
	std::string _name;
	uint32_t _readGroupID;
	std::string _readNamePrefix;
	int _readXPos, _readYPos;
	bool _isInitialized;

	//read length
	TReadLengthDistribution* _readLengthDist;
	bool _readLengthInitialized;

	//qualities
	TSimulatorQualityDist* _qualityDist;
	bool _qualityDistInitialized;
	std::string _qualDistType;
	TSimulatorQualityTransformation* _qualityTransform;
	bool _qualityTransformInitialized;

	//PMD
	GenotypeLikelihoods::TPMDDoubleStrand _pmdObject;
	bool _hasPMD;

	//contamination
	bool _isContaminated;
	double _contaminationRate;
	TSimulatorReference* contaminationSource;

	//alignment
	BAM::TCigar _cigar;
	BAM::TSamFlags _flags;
	BAM::TAlignment _alignment;
	Base* bases;
	int* phredIntQualities;

	//general functions
	void _simulateQualitiesAndErrors(Base* _bases, int* _qualities, int & len);
	void _applyPMD(Base* _bases, const TReadLength & readLength, const bool isReverseStrand);
	std::string _getNextReadName();
	void _simulateBasesQualities(BAM::TAlignment & alignment, Base* haplotype, const uint64_t pos, const TReadLength & readLength, const bool readIsContaminated, const bool isReverse, TSimulatorQualityTransformation* qualityTransform);

public:
	TSimulatorSingleEndRead(std::string readGroupName, int readGroupNumber, int MaxPrintQual, TRandomGenerator* RandomGenerator, TGenotypeMap & GenoMap);
	virtual ~TSimulatorSingleEndRead();

	bool checkInitialization();
	void setReadLengthDistribution(std::string s, TLog* logfile);
	void setQualityDistribution(std::string s);
	virtual void setQualityTransformation(TSimulatorQualityTransformParameters & parameters, TLog* logfile);
	void setPMD(const std::string & pmdStringCT, const std::string & pmdStringGA);
	void setContamination(double rate, TSimulatorReference* source);

	std::string name() const{ return _name; };
	std::string type() const{ return _type; };
	double meanReadLength(){
		if(!_readLengthInitialized) throw "Read length distribution not initialized!";
		return _readLengthDist->mean();
	};
	double maxReadLength(){
		if(!_readLengthInitialized) throw "Read length distribution not initialized!";
		return _readLengthDist->max();
	};

	virtual void simulate(Base* haplotype, const uint32_t & refID, const uint32_t & pos, TSimulatorBamFile & bamFile);

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

	std::vector<BAM::TAlignment*> bamAlignmentSecondMates;
	std::vector<BAM::TAlignment*> bamAlignmentSecondMates_idle;

	TSimulatorQualityTransformation* qualityTransform_secondMate;

public:
	TSimulatorPairedEndReads(std::string readGroupName, int readGroupNumber, int MaxPrintQual, TRandomGenerator* RandomGenerator, TGenotypeMap & GenoMap);
	~TSimulatorPairedEndReads();

	void setQualityTransformation(TSimulatorQualityTransformParameters & parameters, TLog* logfile);
	void simulate(Base* haplotype, const uint32_t refID, const uint32_t & pos, TSimulatorBamFile & bamFile);
	void writeUnwrittenAlignments(const long & pos, TSimulatorBamFile & bamFile);
};

}; //end namespace

#endif /* TSIMULATORREAD_H_ */
