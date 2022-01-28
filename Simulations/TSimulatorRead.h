/*
 * TSimulatorQuality.h
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREAD_H_
#define TSIMULATORREAD_H_

#include "TParameters.h"
#include "TPostMortemDamage.h"
#include "TReadGroups.h"
#include "TSequencingErrorModels.h"
#include "TSimulatorAuxiliaryTools.h"
#include "TSimulatorQuality.h"
#include "TSimulatorReadLength.h"
#include <memory>

namespace Simulations {

using genometools::Base;

//-------------------------------
// TSimulatorSingleEndRead
//-------------------------------
class TSimulatorSingleEndRead {
protected:
	TRandomGenerator *_randomGenerator;

	const BAM::TReadGroup &_readGroup;
	std::string _readNamePrefix;
	int _readXPos = 1;
	int _readYPos = 1;

	// read length
	std::unique_ptr<TReadLengthDistribution> _readLengthDist;

	// qualities
	std::unique_ptr<TSimulatorQualityDist> _qualityDist;
	std::unique_ptr<TSimulatorQualityDist> _mappingQualityDist;

	// PMD
	GenotypeLikelihoods::TPMDType const *_pmd = nullptr;

	// Recal
	GenotypeLikelihoods::TSequencingErrorModelsOneReadGroup const *_recal = nullptr;

	// contamination
	double _contaminationRate = 0.;
	TSimulatorReference *_contaminationSource = nullptr;

	// alignment
	BAM::TSamFlags _flags;
	BAM::TCigar _cigar;
	BAM::TAlignment _alignment;

	// function initialize
	std::unique_ptr<TSimulatorQualityDist> _initializeQualityDistribution(std::string s);

	// general functions
	void _simulateQualitiesAndErrors(Base *_bases, int *_qualities, int &len);
	void _applyPMD(std::vector<Base> &bases, const TReadLength &readLength, const bool isReverseStrand);
	std::string _getNextReadName();
	void _simulateBasesQualities(BAM::TAlignment &alignment, const std::vector<Base>& haplotype, const uint64_t pos,
				     const TReadLength &readLength, bool readIsContaminated);//, TSimulatorQualityTransformation *qualityTransform);

public:
	TSimulatorSingleEndRead(const BAM::TReadGroup &ReadGroup, coretools::TRandomGenerator *RandomGenerator);
	virtual ~TSimulatorSingleEndRead() = default;

	bool checkInitialization();
	void setReadLengthDistribution(std::string s, TLog *logfile);
	void setQualityDistribution(std::string s);
	void setMappingQualityDistribution(std::string s);
	void setPMD(GenotypeLikelihoods::TPMDType const *Pmd);
	virtual void setQualityTransformation(GenotypeLikelihoods::TSequencingErrorModelsOneReadGroup const *Recal);
	void setContamination(double rate, TSimulatorReference *source);

	std::string name() const { return _readGroup.name_ID; };
	virtual std::string type() const {return "single-end";}
	double meanReadLength() {
		if (!_readLengthDist) throw "Read length distribution not initialized!";
		return _readLengthDist->mean();
	};
	double maxReadLength() {
		if (!_readLengthDist) throw "Read length distribution not initialized!";
		return _readLengthDist->max();
	};

	virtual void simulate(const std::vector<Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile);

	void printDetails(TLog *logfile);
	virtual void writeUnwrittenAlignments(long, TSimulatorBamFile &){};
};

//-------------------------------
// TSimulatorPairedEndReads
//-------------------------------
class TSimulatorPairedEndReads : public TSimulatorSingleEndRead {
private:
	BAM::TSamFlags _mateFlags;
	// BAM::TAlignment _mate;

	std::vector<BAM::TAlignment> bamAlignmentSecondMates;

public:
	TSimulatorPairedEndReads(const BAM::TReadGroup &, coretools::TRandomGenerator *RandomGenerator);

	void simulate(const std::vector<Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile) override;
	void writeUnwrittenAlignments(long pos, TSimulatorBamFile &bamFile) override;
	std::string type() const override {return "paired-end";} 
};

}; // namespace Simulations

#endif /* TSIMULATORREAD_H_ */
