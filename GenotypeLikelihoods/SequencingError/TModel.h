/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include "TEpsilon.h"
#include "coretools/Main/TError.h"
#include "genometools/Genotypes/Containers.h"
#include "TReadGroupInfo.h"
#include "TRho.h"
#include <stdexcept>

namespace BAM {
struct TSequencedData;
class TAlignment;
}

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------
// TModel
// pure abstract base class
//--------------------------------------------------------------------
struct TModel {
	virtual ~TModel()                                                                               = default;
	// Per Base
	virtual coretools::PhredInt phredInt(const BAM::TSequencedData &data) const noexcept            = 0;
	virtual genometools::TBaseLikelihoods P_dij(const BAM::TSequencedData &data) const noexcept     = 0;
	// Per Alignment
	virtual void simulate(BAM::TAlignment &aln) const noexcept                                      = 0;
	virtual void recalibrate(BAM::TAlignment &aln) const noexcept                                   = 0;
	// Model Info
	virtual bool recalibrates() const noexcept                                                      = 0;
	virtual BAM::RGInfo::TInfo info() const                                                         = 0;
	virtual void addtoRho(const BAM::TSequencedData &data, coretools::Probability P_g_I_d,
						  const genometools::TBaseProbabilities &P_bbar_I_d) noexcept               = 0;
	virtual void addToEpsilon(const BAM::TSequencedData &data, const TGenotypeFloats &P_g_I_ds,
							  const TGenotypeFloats &P_bbar_I_gds, bool updateJF, bool isInvariant) = 0;
	virtual TRho *rho()                                                                             = 0;
	virtual TEpsilon *epsilon()                                                                     = 0;
};

//------------------------------------------------
// TModelNoRecal
//------------------------------------------------
struct TNoRecal final : public TModel {
	// Per Base
	coretools::PhredInt phredInt(const BAM::TSequencedData &data) const noexcept override;
	genometools::TBaseLikelihoods P_dij(const BAM::TSequencedData &data) const noexcept override;

	// Per Alignment
	void simulate(BAM::TAlignment &aln) const noexcept override;
	void recalibrate(BAM::TAlignment &aln) const noexcept override;

	// Model Info
	bool recalibrates() const noexcept override { return false; };
	BAM::RGInfo::TInfo info() const override { return "-"; };

	void addtoRho(const BAM::TSequencedData &, coretools::Probability,
				  const genometools::TBaseProbabilities &) noexcept override {}
	void addToEpsilon(const BAM::TSequencedData &, const TGenotypeFloats &, const TGenotypeFloats &, bool,
					  bool) noexcept override {}
	TRho *rho() override { throw coretools::TDevError("TNoRecal does not have rho!"); }
	TEpsilon *epsilon() override {throw coretools::TDevError("TNoRecal does not have epsilon!");}
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
	coretools::PhredInt phredInt(const BAM::TSequencedData &data) const noexcept override;
	genometools::TBaseLikelihoods P_dij(const BAM::TSequencedData &data) const noexcept override;

	void simulate(BAM::TAlignment &aln) const noexcept override;
	void recalibrate(BAM::TAlignment &aln) const noexcept override;

	// Model Info
	bool recalibrates() const noexcept override { return true; };
	BAM::RGInfo::TInfo info() const override;

	void addtoRho(const BAM::TSequencedData &data, coretools::Probability P_g_I_d,
	              const genometools::TBaseProbabilities &P_bbar_I_d) noexcept override {
		_rho.add(data, P_g_I_d, P_bbar_I_d);
	}
	void addToEpsilon(const BAM::TSequencedData &data, const TGenotypeFloats &P_g_I_ds,
					  const TGenotypeFloats &P_bbar_I_gds, bool updateJF, bool isInvariant) noexcept override {
		_epsilon.add(data, P_g_I_ds, P_bbar_I_gds, updateJF, isInvariant);
	}
	TRho *rho() noexcept override { return &_rho; }
	TEpsilon* epsilon() noexcept override {return &_epsilon;}
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
