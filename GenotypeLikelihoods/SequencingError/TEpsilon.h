/*
 * TEpsilon.h
 *
 *  Created on: Jul 19, 2022
 *      Author: Andreas
 */

#ifndef TEPSILON_H_
#define TEPSILON_H_

#include <armadillo>

#include "SequencingError/TFunctions.h"
#include "coretools/Types/probability.h"
#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/Genotype.h"

namespace RecalEstimatorTools {class TRecalDataTable;}

namespace GenotypeLikelihoods {
namespace SequencingError {

using TGenotypeFloats = coretools::TStrongArray<float, genometools::Genotype, 10>;

class TEpsilon {
	TFunctions _functions;

	double _Q      = 0.;
	double _oldQ   = 0.;
	double _maxF   = 0.;
	double _maxJxF = 0.;

	arma::mat _Jacobian;
	arma::vec _F;
	arma::mat _JxF;
	size_t _numSitesAdded = 0;
	bool _accepted        = false;

	coretools::Probability _calcErrorRate(const BAM::TSequencedData &data, std::vector<T1stDerivative> &der1,
												   std::vector<T2ndDerivative> &der2) const noexcept;

	struct Genotypes {
		static constexpr std::array<genometools::Genotype, 4> homo{
			genometools::Genotype::AA, genometools::Genotype::CC, genometools::Genotype::GG, genometools::Genotype::TT};
		static constexpr std::array<genometools::Genotype, 6> het{genometools::Genotype::AC, genometools::Genotype::AG,
																  genometools::Genotype::AT, genometools::Genotype::CG,
																  genometools::Genotype::CT, genometools::Genotype::GT};
	};

	void _addToQ(const BAM::TSequencedData &data, const TGenotypeFloats &P_g_I_ds, const TGenotypeFloats &P_bbar_I_gds,
				 bool isInvariant);

	void _addToQJF(const BAM::TSequencedData &data, const TGenotypeFloats &P_g_I_ds,
				   const TGenotypeFloats &P_bbar_I_gds, bool isInvariant);

public:
	TEpsilon(std::string_view Def);
	TEpsilon(const BAM::RGInfo::TInfo & info);

	void init(const RecalEstimatorTools::TRecalDataTable &DataTable, size_t MinData);

	coretools::Probability calcErrorRate(const BAM::TSequencedData &data) const noexcept; 
	double deltaQ() const noexcept {return _Q - _oldQ;};
	double Q() const noexcept {return _Q;};
	double maxF() const noexcept {return _maxF;};
	double maxChange() const noexcept {return _maxJxF;}
	bool accepted() const noexcept {return _accepted; }

	void add(const BAM::TSequencedData &data, const TGenotypeFloats &P_g_I_ds, const TGenotypeFloats &P_bbar_I_gds,
			 bool updateJF, bool isInvariant);
	void solveJxF();
	void propose(double lambda);
	bool acceptOrReject();
	void adjust();

	void log() const;
	BAM::RGInfo::TInfo info() const;
};
} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* TRECALIBRATIONEMMODEL_H_ */
