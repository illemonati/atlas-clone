/*
 * TReadSimulator.h
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#ifndef TSIMULATORREAD_H_
#define TSIMULATORREAD_H_

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "TOutputBamFile.h"
#include "PMD/TModel.h"
#include "SequencingError/TModels.h"
#include "TAlignment.h"
#include "TCigar.h"
#include "TReadGroupInfo.h"
#include "TSamFlags.h"
#include "coretools/Main/TLog.h"
#include "coretools/Math/TCategoricalDistribution.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"

namespace Simulations { class TSimulatorReference; }

namespace Simulations {


//-------------------------------------
// TSimulatorRead
// a pure abstract base-class
//-------------------------------------
class TReadSimulator{
protected:
	const BAM::TReadGroup *_readGroup;
	//const TReadGroupInfoEntry & _readGroupInfo;
	std::string _readNamePrefix;

	// required distributions
	coretools::probdist::TCategoricalDistribution<size_t> _fragmentLengthDistr;
	coretools::probdist::TCategoricalDistribution<genometools::PhredIntProbability> _qualityDist;
	coretools::probdist::TCategoricalDistribution<genometools::PhredIntProbability> _mappingQualityDist;

	// Additional info
	int _readXPos = 1;
	int _readYPos = 1;
	std::unique_ptr<coretools::probdist::TCategoricalDistribution<size_t>> _softClipDist5;
	std::unique_ptr<coretools::probdist::TCategoricalDistribution<size_t>> _softClipDist3;
	const GenotypeLikelihoods::PMD::TModel *_pmd;
	const GenotypeLikelihoods::SequencingError::RGModels _recal;
	coretools::Probability _duplicationRate, _duplicationRateAmongSimulated;

	// contamination
	double _contaminationRate = 0.;
	TSimulatorReference *_contaminationSource = nullptr;

	// alignment
	BAM::TSamFlags _flags;
	BAM::TAlignment _alignment;

	//initialization functions
	template <typename Distr>
	void _initDistribution(Distr & Dist, const BAM::RGInfo::TReadGroupInfoEntry & RGInfo, const BAM::RGInfo::InfoType & Info){
		coretools::instances::logfile().list(coretools::str::capitalizeFirst(BAM::RGInfo::infos[Info].description), ": ", RGInfo.getString(Info));
		Dist.set(RGInfo.getString(Info));
	};

	// general functions
	double _calcMeanReadLength(size_t maxLen) const;
	std::string _getNextReadName();
	void _simulateAlignmentDetails(const genometools::TGenomePosition & Position);
	bool _simulateContamination();
	void _addSoftclippedBases(std::vector<genometools::Base> & bases, const std::unique_ptr<coretools::probdist::TCategoricalDistribution<size_t>> & softClippedDist, BAM::TCigar & Cigar);
	void _simulateBasesQualities(BAM::TAlignment &alignment, const std::vector<genometools::Base> &haplotype, size_t fragmentLength,
								 size_t readLength, bool readIsContaminated);

	virtual void _simulate(const genometools::TGenomePosition & Position, const std::vector<genometools::Base> & Haplotype) = 0;
	virtual void _writeSimulatedAlignments(BAM::TOutputBamFile & BamFile) = 0;

public:
	TReadSimulator(const BAM::TReadGroup & ReadGroup, const BAM::RGInfo::TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal);
	virtual ~TReadSimulator() = default;

	//setters
	void setPMD(GenotypeLikelihoods::PMD::TModel const *Pmd);
	void setContamination(double rate, TSimulatorReference *source);

	//simulate
	void simulate(const genometools::TGenomePosition & Position, const std::vector<genometools::Base> & Haplotype, BAM::TOutputBamFile &BamFile);

	//getters
	std::string name() const { return _readGroup->name_ID; };
	[[nodiscard]] virtual double meanReadLength() const = 0;
	double maxFragmentLength() {
		return _fragmentLengthDistr.max();
	};
};

//-------------------------------
// TSimulatorSingleEndRead
//-------------------------------
class TReadSimulatorSingleEnd final : public TReadSimulator {
private:
	coretools::StrictlyPositive<size_t> _numCycles;

	void _simulate(const genometools::TGenomePosition & Position, const std::vector<genometools::Base> & Haplotype) override;
	void _writeSimulatedAlignments(BAM::TOutputBamFile & BamFile) override;

public:
	TReadSimulatorSingleEnd(const BAM::TReadGroup & ReadGroup, const BAM::RGInfo::TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal);
	~TReadSimulatorSingleEnd() = default;

	[[nodiscard]] double meanReadLength() const override;
};

//-------------------------------
// TSimulatorPairedEndReads
//-------------------------------
class TReadSimulatorPairedEnd final : public TReadSimulator {
private:
	BAM::TAlignment _secondMate;
	BAM::TSamFlags _mateFlags;
	std::array<coretools::StrictlyPositive<size_t>, 2> _numCycles;

	void _simulate(const genometools::TGenomePosition & Position, const std::vector<genometools::Base> & Haplotype) override;
	void _writeSimulatedAlignments(BAM::TOutputBamFile & BamFile) override;

public:
	TReadSimulatorPairedEnd(const BAM::TReadGroup & ReadGroup, const BAM::RGInfo::TReadGroupInfoEntry & RGInfo, const GenotypeLikelihoods::PMD::TModel & Pmd, const GenotypeLikelihoods::SequencingError::RGModels& Recal);
	~TReadSimulatorPairedEnd() = default;

	[[nodiscard]] double meanReadLength() const override;
};

}; // namespace Simulations

#endif /* TSIMULATORREAD_H_ */
