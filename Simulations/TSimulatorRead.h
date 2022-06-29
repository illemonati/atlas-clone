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

#include "GenotypeTypes.h"
#include "TAlignment.h"
#include "TCigar.h"
#include "TReadGroups.h"
#include "TSamFlags.h"
#include "TSimulatorQuality.h"
#include "TSimulatorReadLength.h"
#include "TSimulatorSoftClip.h"

namespace GenotypeLikelihoods { class TPMDType; }
namespace GenotypeLikelihoods { namespace SequencingError { class TModel; } }
namespace Simulations { class TSimulatorBamFile; }
namespace Simulations { class TSimulatorReference; }

namespace Simulations {


//-------------------------------
// TSimulatorSingleEndRead
//-------------------------------
class TSimulatorSingleEndRead {
protected:
	const BAM::TReadGroup &_readGroup;
	std::string _readNamePrefix;
	int _readXPos = 1;
	int _readYPos = 1;

	// read length
	std::unique_ptr<TReadLengthDistribution> _readLengthDist;

	// qualities
	std::unique_ptr<TSimulatorQualityDist> _qualityDist;
	std::unique_ptr<TSimulatorQualityDist> _mappingQualityDist;
	std::unique_ptr<TSimulatorSoftClipDist> _softClipDist;

	GenotypeLikelihoods::TPMDType const *_pmd = nullptr;
	std::array<GenotypeLikelihoods::SequencingError::TModel const *, 2> _recal;

	// contamination
	double _contaminationRate = 0.;
	TSimulatorReference *_contaminationSource = nullptr;

	// alignment
	BAM::TSamFlags _flags;
	BAM::TCigar _cigar;
	BAM::TAlignment _alignment;

	// function initialize
	std::unique_ptr<TSimulatorQualityDist> _initializeQualityDistribution(std::string s);
	std::unique_ptr<TSimulatorSoftClipDist> _initializeSoftClipDistribution(std::string s);

	// general functions
	void _simulateQualitiesAndErrors(genometools::Base *_bases, int *_qualities, int &len);
	std::string _getNextReadName();
	void _simulateBasesQualities(BAM::TAlignment &alignment, const std::vector<genometools::Base>& haplotype, const uint64_t pos,
				     const TReadLength &readLength, bool readIsContaminated);//, TSimulatorQualityTransformation *qualityTransform);

public:
	TSimulatorSingleEndRead(const BAM::TReadGroup &ReadGroup);
	virtual ~TSimulatorSingleEndRead() = default;

	bool checkInitialization();
	void setReadLengthDistribution(std::string s);
	void setQualityDistribution(std::string s);
	void setMappingQualityDistribution(std::string s);
	void setSoftClipDistribution(std::string s);
	void setPMD(GenotypeLikelihoods::TPMDType const *Pmd);
	void setRecal(GenotypeLikelihoods::SequencingError::TModel const *Recal1, GenotypeLikelihoods::SequencingError::TModel const *Recal2);
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

	virtual void simulate(const std::vector<genometools::Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile);

	void printDetails();
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
	TSimulatorPairedEndReads(const BAM::TReadGroup &);

	void simulate(const std::vector<genometools::Base>& haplotype, uint32_t refID, uint32_t pos, TSimulatorBamFile &bamFile) override;
	void writeUnwrittenAlignments(long pos, TSimulatorBamFile &bamFile) override;
	std::string type() const override {return "paired-end";} 
};

}; // namespace Simulations

#endif /* TSIMULATORREAD_H_ */
