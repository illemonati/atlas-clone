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
#include "TSimulatorQuality.h"
#include "TSimulatorQualityTransformation.h"
#include "TReadGroups.h"

namespace Simulations{

//-------------------------------
//TSimulatorSingleEndRead
//-------------------------------
class TSimulatorSingleEndRead{
protected:
	TRandomGenerator* _randomGenerator;
	GenotypeLikelihoods::TGenotypeMap& _genoMap;
	BAM::TQualityMap& _qualMap;

	const BAM::TReadGroup& _readGroup;
	std::string _type;
	std::string _readNamePrefix;
	int _readXPos, _readYPos;
	bool _isInitialized;

	//read length
	std::unique_ptr<TReadLengthDistribution> _readLengthDist;

	//qualities
	std::unique_ptr<TSimulatorQualityDist> _qualityDist;
	std::unique_ptr<TSimulatorQualityDist> _mappingQualityDist;

	TSimulatorQualityTransformation* _qualityTransform;
	bool _qualityTransformInitialized;

	//PMD
	GenotypeLikelihoods::TPMDType const* _pmdObject;
	bool _hasPMD;

	//contamination
	bool _isContaminated;
	double _contaminationRate;
	TSimulatorReference* contaminationSource;

	//alignment
	BAM::TSamFlags _flags;
	BAM::TCigar _cigar;
	BAM::TAlignment _alignment;

	//function initialize
	void _initializeQualityDistribution(std::string s, std::unique_ptr<TSimulatorQualityDist> & pointer);

	//general functions
	void _simulateQualitiesAndErrors(Base* _bases, int* _qualities, int & len);
	void _applyPMD(std::vector<Base>& bases, const TReadLength & readLength, const bool isReverseStrand);
	std::string _getNextReadName();
	void _simulateBasesQualities(BAM::TAlignment & alignment, Base* haplotype, const uint64_t pos, const TReadLength & readLength, const bool readIsContaminated, TSimulatorQualityTransformation* qualityTransform);

public:
	TSimulatorSingleEndRead(const BAM::TReadGroup& ReadGroup, TRandomGenerator* RandomGenerator, GenotypeLikelihoods::TGenotypeMap & GenoMap);
	virtual ~TSimulatorSingleEndRead();

	bool checkInitialization();
	void setReadLengthDistribution(std::string s, TLog* logfile);
	void setQualityDistribution(std::string s);
	void setMappingQualityDistribution(std::string s);
	void setPMD(GenotypeLikelihoods::TPMDType const* PmdObject);
	virtual void setQualityTransformation(TSimulatorQualityTransformParameters & parameters, TLog* logfile);
	void setContamination(double rate, TSimulatorReference* source);

	std::string name() const{ return _readGroup.name_ID; };
	std::string type() const{ return _type; };
	double meanReadLength(){
		if(!_readLengthDist) throw "Read length distribution not initialized!";
		return _readLengthDist->mean();
	};
	double maxReadLength(){
		if(!_readLengthDist) throw "Read length distribution not initialized!";
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
	TSimulatorPairedEndReads(const BAM::TReadGroup&, TRandomGenerator* RandomGenerator, GenotypeLikelihoods::TGenotypeMap & GenoMap);
	~TSimulatorPairedEndReads();

	void setQualityTransformation(TSimulatorQualityTransformParameters & parameters, TLog* logfile);
	void simulate(Base* haplotype, const uint32_t refID, const uint32_t & pos, TSimulatorBamFile & bamFile);
	void writeUnwrittenAlignments(const long & pos, TSimulatorBamFile & bamFile);
};

}; //end namespace

#endif /* TSIMULATORREAD_H_ */
