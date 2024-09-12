/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include "TEpsilon.h"
#include "genometools/Genotypes/Containers.h"
#include "TReadGroupInfo.h"
#include "TRho.h"

namespace BAM {
struct TSequencedBase;
class TAlignment;
}

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------
// TModel
// pure abstract base class
//--------------------------------------------------------------------
struct TModel {
	virtual ~TModel()                                                                                 = default;
	// Per Base
	virtual coretools::PhredInt phredInt(const BAM::TSequencedBase &data) const noexcept = 0;
	virtual genometools::TBaseLikelihoods P_dij(const BAM::TSequencedBase &data) const noexcept                    = 0;
	// Per Alignment
	virtual void simulate(BAM::TAlignment &aln) const noexcept                                        = 0;
	virtual void recalibrate(BAM::TAlignment &aln) const noexcept                                     = 0;
	// Model Info
	virtual bool recalibrates() const noexcept                                                        = 0;
	virtual BAM::RGInfo::TInfo info() const                                                           = 0;
	virtual TRho *rho() noexcept                                                                      = 0;
	virtual TEpsilon *epsilon() noexcept                                                              = 0;
};

//------------------------------------------------
// TModelNoRecal
//------------------------------------------------
struct TNoRecal final : public TModel {
	// Per Base
	coretools::PhredInt phredInt(const BAM::TSequencedBase &data) const noexcept override;
	genometools::TBaseLikelihoods P_dij(const BAM::TSequencedBase &data) const noexcept override;

	// Per Alignment
	void simulate(BAM::TAlignment &aln) const noexcept override;
	void recalibrate(BAM::TAlignment &aln) const noexcept override;

	// Model Info
	bool recalibrates() const noexcept override { return false; };
	BAM::RGInfo::TInfo info() const override { return BAM::RGInfo::TInfo{}; };
	TRho *rho() noexcept override {return nullptr;}
	TEpsilon *epsilon() noexcept override {return nullptr;}
};

//------------------------------------------------
// TModelRecal
//------------------------------------------------
class TWithRecal final : public TModel {
private:
	TRho _rho;
	TEpsilon _epsilon;
public:
	TWithRecal(std::string_view EpsilonDef, std::string_view RhoDef = "")
		: _rho(RhoDef), _epsilon(EpsilonDef) {}
	TWithRecal(const BAM::RGInfo::TInfo &info) : _rho(info["rho"]), _epsilon(info) {}

	// Per Base
	coretools::PhredInt phredInt(const BAM::TSequencedBase &data) const noexcept override;
	genometools::TBaseLikelihoods P_dij(const BAM::TSequencedBase &data) const noexcept override;

	void simulate(BAM::TAlignment &aln) const noexcept override;
	void recalibrate(BAM::TAlignment &aln) const noexcept override;

	// Model Info
	bool recalibrates() const noexcept override { return true; };
	BAM::RGInfo::TInfo info() const override;
	TRho* rho() noexcept override {return &_rho;}
	TEpsilon* epsilon() noexcept override {return &_epsilon;}
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
