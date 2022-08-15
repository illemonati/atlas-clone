/*
 * TSimulatorQuality.h
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREAD_H_
#define TSIMULATORREAD_H_

#include <stdint.h>
#include <array>
#include <memory>
#include <string>
#include <vector>
#include "../BAM/TReadGroupInfo.h"
#include "GenotypeTypes.h"
#include "TAlignment.h"
#include "TCigar.h"
#include "TReadGroups.h"
#include "TSamFlags.h"
#include "PhredProbabilityTypes.h"
#include "TCategoricalDistribution.h"
#include "TFragmentLengthDistribution.h"
#include "TReadGroupInfo.h"

namespace GenotypeLikelihoods { class TPMDType; }
namespace GenotypeLikelihoods { namespace SequencingError { class TModel; } }
namespace Simulations { class TSimulatorBamFile; }
namespace Simulations { class TSimulatorReference; }

namespace Simulations {

using genometools::Base;
using genometools::PhredIntProbability;
using coretools::probdist::TCategoricalDistribution;


//-------------------------------
// TSimulatorSingleEndRead
//-------------------------------
class TSimulatorSingleEndRead {
protected:
	const BAM::TReadGroup &_readGroup;
	const TSimulatorReadGroupInfoEntry & _readGroupInfo;
	std::string _readNamePrefix;

	// required distributions
	TFragmentLengthDistribution _fragmentLengthDist;
	TCategoricalDistribution<PhredIntProbability> _qualityDist;
	TCategoricalDistribution<PhredIntProbability> _mappingQualityDist;

	// Additional info
	coretools::StrictlyPositive<uint16_t> _numCycles;
	int _readXPos = 1;
	int _readYPos = 1;
	std::unique_ptr<TCategoricalDistribution<uint16_t>> _softClipDist5;
	std::unique_ptr<TCategoricalDistribution<uint16_t>> _softClipDist3;
	GenotypeLikelihoods::TPMDType const *_pmd = nullptr;
	std::array<GenotypeLikelihoods::SequencingError::TModel const *, 2> _recal;

	// contamination
	double _contaminationRate = 0.;
	TSimulatorReference *_contaminationSource = nullptr;

	// alignment
	BAM::TSamFlags _flags;
	BAM::TCigar _cigar;
	BAM::TAlignment _alignment;

	// general functions
	void _simulateQualitiesAndErrors(Base *_bases, int *_qualities, int &len);
	std::string _getNextReadName();
	void _addSoftclippedBases(std::vector<Base> & bases, const std::unique_ptr<TCategoricalDistribution<uint16_t>> & softClippedDist);
	void _simulateBasesQualities(BAM::TAlignment &alignment, const std::vector<Base>& haplotype, const uint64_t pos,
				     const TReadAndFragmentLength &readLength, bool readIsContaminated);//, TSimulatorQualityTransformation *qualityTransform);

public:
	TSimulatorSingleEndRead(const BAM::TReadGroup & ReadGroup, const TSimulatorReadGroupInfoEntry & RGInfo);
	virtual ~TSimulatorSingleEndRead() = default;

	void setPMD(GenotypeLikelihoods::TPMDType const *Pmd);
	void setRecal(GenotypeLikelihoods::SequencingError::TModel const *Recal1, GenotypeLikelihoods::SequencingError::TModel const *Recal2);
	void setContamination(double rate, TSimulatorReference *source);

	std::string name() const { return _readGroup.name_ID; };
	virtual std::string type() const {return "single-end";}
	double meanReadLength() {
		return _fragmentLengthDist.mean();
	};
	double maxReadLength() {
		return _fragmentLengthDist.max();
	};

	virtual void simulate(const std::vector<Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile);

	void printDetails(double frequency);
	virtual void writeUnwrittenAlignments(long, TSimulatorBamFile &){};
};

//-------------------------------
// TSimulatorPairedEndReads
//-------------------------------
class TSimulatorPairedEndReads : public TSimulatorSingleEndRead {
private:
	BAM::TSamFlags _mateFlags;
	// BAM::TAlignment _mate;
	uint16_t _numCyclesSecond;

	std::vector<BAM::TAlignment> bamAlignmentSecondMates;

public:
	TSimulatorPairedEndReads(const BAM::TReadGroup &, const uint16_t NumCyclesFirst, const uint16_t NumCyclesSecond);

	void simulate(const std::vector<genometools::Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile) override;
	void writeUnwrittenAlignments(long pos, TSimulatorBamFile &bamFile) override;
	std::string type() const override {return "paired-end";} 
};

}; // namespace Simulations

#endif /* TSIMULATORREAD_H_ */
