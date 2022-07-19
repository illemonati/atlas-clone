/*
 * TRecalibrationEMModel.h
 *
 *  Created on: Mar 7, 2019
 *      Author: phaentu
 */

#ifndef TRECALIBRATIONEMMODEL_H_
#define TRECALIBRATIONEMMODEL_H_

#include <limits>
#include <stdint.h>
#include <string>

#include "GenotypeTypes.h"
#include "PhredProbabilityTypes.h"
#include "TGenotypeData.h"
#include "TSequencedBase.h"
#include "TStrongArray.h"
#include "TEpsilon.h"
#include "probability.h"

namespace BAM {
class TSequencedBase;
}
namespace GenotypeLikelihoods {
namespace RecalEstimatorTools {
class TRecalDataTable;
}
} // namespace GenotypeLikelihoods

namespace GenotypeLikelihoods {
namespace SequencingError {

//--------------------------------------------------------------------
// TRho
//--------------------------------------------------------------------
class TRho {
private:
	coretools::TStrongArray<TBaseProbabilities, genometools::Base> _rho{
		{coretools::TStrongArray<double, genometools::Base>{{0., 1. / 3, 1. / 3, 1. / 3}},
		 coretools::TStrongArray<double, genometools::Base>{{1. / 3, 0., 1. / 3, 1. / 3}},
		 coretools::TStrongArray<double, genometools::Base>{{1. / 3, 1. / 3, 0., 1. / 3}},
		 coretools::TStrongArray<double, genometools::Base>{{1. / 3, 1. / 3, 1. / 3, 0.}}}}; //[from][to]

	coretools::TStrongArray<coretools::TStrongArray<double, genometools::Base>, genometools::Base> _rhoSum{};
public:
	TRho() = default;
	TRho(const std::string &Def);
	TRho(const TRho &other) = default;
	TRho &operator=(const TRho &other) = default;


	const TBaseProbabilities& operator[](genometools::Base from) const noexcept {
		return _rho[from];
	}
	std::string getDefinition() const noexcept;

	// functions used to estimate
	void add(genometools::Base base, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept;
	void estimate() noexcept;
};


//--------------------------------------------------------------------
// TModel
// pure abstract base class
//--------------------------------------------------------------------
class TModel {
public:
	virtual ~TModel()                                                                                    = default;
	virtual bool estimatable() const noexcept                                                            = 0;
	virtual bool recalibrates() const noexcept                                                           = 0;
	virtual coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept          = 0;
	virtual genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept = 0;
	virtual TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept          = 0;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept                                      = 0;
	virtual std::string getCovariateDefinition() const noexcept                                          = 0;
	virtual std::string getRhoDefinition() const noexcept                                                = 0;
};

//------------------------------------------------
// TModelNoRecal
//------------------------------------------------
class TModelNoRecal final : public TModel {
public:
	virtual bool estimatable() const noexcept override { return false; };
	virtual bool recalibrates() const noexcept override { return false; };

	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept override;
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept override;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept override;

	virtual std::string getCovariateDefinition() const noexcept override { return "-"; };
	virtual std::string getRhoDefinition() const noexcept override { return "-"; };
};

//------------------------------------------------
// TModelRecal
//------------------------------------------------
class TModelRecal final : public TModel {
private:
	TRho _rho;
	TEpsilon _epsilon;
public:
	TModelRecal(const std::string& RhoDef, const std::string &EpsilonDef);

	bool estimatable() const noexcept override { return true; };
	bool recalibrates() const noexcept override { return true; };
	std::string getCovariateDefinition() const noexcept override {return _epsilon.getDefinition();};
	std::string getRhoDefinition() const noexcept override { return _rho.getDefinition(); };

	// get error rates
	coretools::Probability getErrorRate(const BAM::TSequencedBase &base) const noexcept override;
	genometools::PhredIntProbability getPhredInt(const BAM::TSequencedBase &base) const noexcept override;
	TBaseLikelihoods getBaseLikelihoods(const BAM::TSequencedBase &base) const noexcept override;
	virtual void simulate(BAM::TSequencedBase &base) const noexcept override;

	// functions to estimate
	void checkOrInit(const RecalEstimatorTools::TRecalDataTable &DataTable) { _epsilon.checkOrInit(DataTable); };

	// functions to estimate rho
	void addToRho(const BAM::TSequencedBase &data, coretools::Probability P_g_I_d, const TBaseProbabilities &P_bbar_I_d) noexcept {_rho.add(data.base, P_g_I_d, P_bbar_I_d);}; 
	void estimateRho() noexcept {_rho.estimate();};

	// functions to estimate betas
	double resetQ() noexcept {return _epsilon.resetQ();}
	void addToEpsilon(const BAM::TSequencedBase &base, coretools::Probability P_g_I_d, coretools::Probability P_bbar_I_gd, bool update) {
		_epsilon.addToEpsilon(base, P_g_I_d, P_bbar_I_gd, update);
	}
	void solveJxF() {_epsilon.solveJxF();}
	void proposeNewParameters(double lambda) {_epsilon.proposeNewParameters(lambda);}
	bool acceptProposedParametersBasedOnQ() {return _epsilon.acceptProposedParametersBasedOnQ();}
	void adjustParametersPostEstimation() {_epsilon.adjustParametersPostEstimation();}
	double maxF() const noexcept {return _epsilon.maxF();}
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
