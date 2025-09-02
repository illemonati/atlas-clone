/*
 * TDistanceCalculator.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#include "TDistanceEstimator.h"
#include "coretools/Files/TLineWriter.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/concatenateString.h"
#include "coretools/algorithms.h"

namespace PopulationTools{

using genometools::TBaseData;
using genometools::Genotype;
using genometools::Base;
using coretools::Probability;
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::user_assert;

//----------------------------------------------------
//TDistanceClass
//----------------------------------------------------
TDistance::TDistance(){
	using DP = DistancePhi;
	using coretools::index;
	//squared difference between genotypes
	_distanceWeight[DP::aa_aa] = 0.0;
	_distanceWeight[DP::ab_aa] = 1.0;
	_distanceWeight[DP::aa_ab] = 1.0;
	_distanceWeight[DP::aa_bb] = 4.0;
	_distanceWeight[DP::ab_ab] = 0.0;
	_distanceWeight[DP::ab_ac] = 1.0;
	_distanceWeight[DP::aa_bc] = 4.0;
	_distanceWeight[DP::ab_cc] = 4.0;
	_distanceWeight[DP::ab_cd] = 4.0;
}

double TDistance::calculateDistance(const TDistanceData & phi){
	using DP = DistancePhi;
	return    phi[DP::aa_aa] * _distanceWeight[DP::aa_aa]
			+ phi[DP::ab_aa] * _distanceWeight[DP::ab_aa]
			+ phi[DP::aa_ab] * _distanceWeight[DP::aa_ab]
			+ phi[DP::aa_bb] * _distanceWeight[DP::aa_bb]
			+ phi[DP::ab_ab] * _distanceWeight[DP::ab_ab]
			+ phi[DP::ab_ac] * _distanceWeight[DP::ab_ac]
			+ phi[DP::aa_bc] * _distanceWeight[DP::aa_bc]
			+ phi[DP::ab_cc] * _distanceWeight[DP::ab_cc]
			+ phi[DP::ab_cd] * _distanceWeight[DP::ab_cd];
}

TDistanceProbMismatch::TDistanceProbMismatch(){
	using DP = DistancePhi;
	//probability that a random allele from each individual is different
	_distanceWeight[DP::aa_aa] = 0.0;
	_distanceWeight[DP::ab_aa] = 0.5;
	_distanceWeight[DP::aa_ab] = 0.5;
	_distanceWeight[DP::aa_bb] = 1.0;
	_distanceWeight[DP::ab_ab] = 0.5;
	_distanceWeight[DP::ab_ac] = 0.75;
	_distanceWeight[DP::aa_bc] = 1.0;
	_distanceWeight[DP::ab_cc] = 1.0;
	_distanceWeight[DP::ab_cd] = 1.0;
}

double TDistanceEuclidian::calculateDistance(const TDistanceData & phi){
	return sqrt(TDistance::calculateDistance(phi));
}

TDistanceUser::TDistanceUser(std::vector<double> vec){
	using coretools::index;
	for(DistancePhi d = DistancePhi::min; d < DistancePhi::max; ++d)
		_distanceWeight[d] = vec[index(d)];
}

//----------------------------------------------------
//TDistanceEstimate
//----------------------------------------------------

coretools::TStrongArray<
	coretools::TStrongArray<coretools::TStrongArray<bool, genometools::Base>, genometools::Genotype>,
	genometools::Genotype>
TEMforDistanceEstimation::_genoToBaseMap() {
	static const coretools::TStrongArray<
		coretools::TStrongArray<coretools::TStrongArray<bool, genometools::Base>, genometools::Genotype>,
		genometools::Genotype>
		map = []() {
			coretools::TStrongArray<
				coretools::TStrongArray<coretools::TStrongArray<bool, genometools::Base>, genometools::Genotype>,
				genometools::Genotype>
				map;

		    for (Genotype g1 = Genotype::min; g1 < Genotype::max; ++g1) {
			    for (Genotype g2 = Genotype::min; g2 < Genotype::max; ++g2) {
				    for (Base b = Base::min; b < Base::max; ++b) { map[g1][g2][b] = false; }
				    map[g1][g2][first(g1)]  = true;
				    map[g1][g2][second(g1)] = true;
				    map[g1][g2][first(g2)]  = true;
				    map[g1][g2][second(g2)] = true;
			    }
		    }
		    return map;
	    }();
	return map;
}

TEMforDistanceEstimation::TEMforDistanceEstimation(){
	//prepare storage
	_phi.fill(0.0);
	_LL = 0.0;
	_old_LL = 0.0;
	_distance = -1.0;

	//read EM parameters
	logfile().startIndent("Parameters of EM algorithm:");
	_maxIterations = parameters().get<int>("iterations", 100);
	logfile().list("Will perform at max ", _maxIterations, " EM iterations. (parameter 'iterations')");
	_minDeltaLL = parameters().get<double>("minDeltaLL", 1e-6);
	logfile().list("Will stop EM when deltaLL < ", _minDeltaLL, ". (parameter 'minDeltaLL')");
	logfile().endIndent();

	//set how to calculate distances
//	distanceWeight = new double[9];
	if(parameters().exists("distWeights")){
		logfile().list("Using user-provided distance weights.");
		const auto vec = parameters().get<std::vector<double>>("distWeights");
		user_assert(vec.size() == 9, "Wrong number of distance weights! Required are nine values for aa/aa, "
									 "aa/ab, ab/aa, aa/bb, ab/ab, ab/ac, aa/bc, ab/cc, ab/cd");

		_distanceObject = std::make_unique<TDistanceUser>(vec);

	} else {
		std::string distType = parameters().get<std::string>("distType", "squaredDiff");
		logfile().list("Using distance type '" + distType + "'.");
		if(distType == "probMismatch"){
			_distanceObject = std::make_unique<TDistanceProbMismatch>();
		} else if(distType == "squaredDiff"){
			_distanceObject = std::make_unique<TDistance>();
		} else if(distType == "euclidian"){
			_distanceObject = std::make_unique<TDistanceEuclidian>();
		} else
			throw coretools::TUserError("Unknown distance type '", distType, "'! Use either squaredDiff, euclidian, or probMismatch.");
	}
	logfile().conclude("Using distance weights " + coretools::str::concatenateString(_distanceObject->weights(), ", ") + ".");

}

void TEMforDistanceEstimation::_guessPi(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2){
	//check sizes are equal
	using coretools::index;
	user_assert(genoQual1.size() == genoQual2.size(), "Provided genotype quality vectors are of different size in TEMforDistanceEstimation::guessPi!");

	//just estimate pi as average posterior probability
	_pi.fill(0.0);

	//now loop over sites
	auto it1 = genoQual1.begin();
	auto it2 = genoQual2.begin();
	for(; it1 != genoQual1.end(); ++it1, ++it2){
		double sum1 = 0.0;
		double sum2 = 0.0;
		for(Genotype g = Genotype::min; g < Genotype::max; ++g){
			sum1 += (Probability) (*it1)[g];
			sum2 += (Probability) (*it2)[g];
		}
		for(Genotype g = Genotype::min; g < Genotype::max; ++g){
			double tmp = (Probability) (*it1)[g] / sum1;
			_pi[first(g)] += tmp;
			_pi[second(g)] += tmp;
			tmp = (Probability) (*it2)[g] / sum2;
			_pi[first(g)] += tmp;
			_pi[second(g)] += tmp;
		}
	}

	//normalize
	normalize(_pi);
}

void TEMforDistanceEstimation::_guessPhi(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2){
	using coretools::index;
	//check sizes are equal
	user_assert(genoQual1.size() == genoQual2.size(), "Provided genotype quality vectors are of different size in TEMforDistanceEstimation::guessPhi!");

	//set to zero
	_phi.fill(0.0);

	//now loop over sites and add posterior probs
	auto it1 = genoQual1.begin();
	auto it2 = genoQual2.begin();
	for(; it1 != genoQual1.end(); ++it1, ++it2){
		double sum1 = 0.0; double sum2 = 0.0;
		for(Genotype g = Genotype::min; g < Genotype::max; ++g){
			sum1 += (Probability) (*it1)[g];
			sum2 += (Probability) (*it2)[g];
		}
		for(Genotype g1 = Genotype::min; g1 < Genotype::max; ++g1){
			double tmp = ((Probability) (*it1)[g1] / sum1);
			for(Genotype g2 = Genotype::min; g2 < Genotype::max; ++g2){
				_phi[distancePhi(g1, g2)] += tmp * ((Probability) (*it2)[g2] / sum2);
			}
		}
	}

	//normalize
	normalize(_phi);
}

void TEMforDistanceEstimation::_fill_K(TBaseData & thesePi){
	using genometools::Base;
	using DP = DistancePhi;
	//normalizing constant for each phi class
	//case of one base
	_K[DP::aa_aa] = 1.0;

	//cases of two bases
	_K[DP::ab_ab] = thesePi[Base::A] * thesePi[Base::C]
		 + thesePi[Base::A] * thesePi[Base::G]
		 + thesePi[Base::A] * thesePi[Base::T]
		 + thesePi[Base::C] * thesePi[Base::G]
		 + thesePi[Base::C] * thesePi[Base::T]
		 + thesePi[Base::G] * thesePi[Base::T];
	_K[DP::aa_ab] = 2.0 * _K[DP::ab_ab]; //account for AC vs CA
	_K[DP::ab_aa] = _K[DP::aa_ab];
	_K[DP::aa_bb] = _K[DP::aa_ab];

	//cases of three bases
	_K[DP::aa_bc] = thesePi[Base::C] * thesePi[Base::G] * thesePi[Base::T]
		 + thesePi[Base::A] * thesePi[Base::G] * thesePi[Base::T]
		 + thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::T]
		 + thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::G];
	_K[DP::aa_bc] = 3.0 * _K[DP::aa_bc]; //account for ways to distribute

	_K[DP::ab_cc] = _K[DP::aa_bc];
	_K[DP::ab_ac] = 2.0 * _K[DP::aa_bc]; //twice as many cases than other cases with three bases!

	//case of four bases: each of the 6 cases is equally likely
	//Note: product of pi's cancels out when calculating P_g_given_phi_pi
	_K[DP::ab_cd] = 6.0;
}

void TEMforDistanceEstimation::_fill_P_g_given_phi_pi(const TDistanceData & thesePhi, TBaseData & thesePi){
	using genometools::Base;
	using coretools::index;
	using DP = DistancePhi;
	using GT = genometools::Genotype;
	//0: case aa/aa (K[0]=1)
	_probGeno[index(GT::AA)][index(GT::AA)] = thesePhi[DP::aa_aa] * thesePi[Base::A];
	_probGeno[index(GT::CC)][index(GT::CC)] = thesePhi[DP::aa_aa] * thesePi[Base::C];
	_probGeno[index(GT::GG)][index(GT::GG)] = thesePhi[DP::aa_aa] * thesePi[Base::G];
	_probGeno[index(GT::TT)][index(GT::TT)] = thesePhi[DP::aa_aa] * thesePi[Base::T];

	//1: cases aa/ab
	double tmp = thesePhi[DP::aa_ab] / _K[DP::aa_aa];
	double tmp2 = tmp * thesePi[Base::A];
	_probGeno[index(GT::AA)][index(GT::AC)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::AA)][index(GT::AG)] = tmp2 * thesePi[Base::G];
	_probGeno[index(GT::AA)][index(GT::AT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::C];
	_probGeno[index(GT::CC)][index(GT::AC)] = tmp2 * thesePi[Base::A];
	_probGeno[index(GT::CC)][index(GT::CG)] = tmp2 * thesePi[Base::G];
	_probGeno[index(GT::CC)][index(GT::CT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::G];
	_probGeno[index(GT::GG)][index(GT::AG)] = tmp2 * thesePi[Base::A];
	_probGeno[index(GT::GG)][index(GT::CG)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::GG)][index(GT::GT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::T];
	_probGeno[index(GT::TT)][index(GT::AT)] = tmp2 * thesePi[Base::A];
	_probGeno[index(GT::TT)][index(GT::CT)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::TT)][index(GT::GT)] = tmp2 * thesePi[Base::G];

	//2: case ab/aa
	tmp = thesePhi[DP::ab_aa] / _K[DP::ab_aa];
	tmp2 = tmp * thesePi[Base::A];
	_probGeno[index(GT::AC)][index(GT::AA)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::AG)][index(GT::AA)] = tmp2 * thesePi[Base::G];
	_probGeno[index(GT::AT)][index(GT::AA)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::C];
	_probGeno[index(GT::AC)][index(GT::CC)] = tmp2 * thesePi[Base::A];
	_probGeno[index(GT::CG)][index(GT::CC)] = tmp2 * thesePi[Base::G];
	_probGeno[index(GT::CT)][index(GT::CC)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::G];
	_probGeno[index(GT::AG)][index(GT::GG)] = tmp2 * thesePi[Base::A];
	_probGeno[index(GT::CG)][index(GT::GG)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::GT)][index(GT::GG)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::T];
	_probGeno[index(GT::AT)][index(GT::TT)] = tmp2 * thesePi[Base::A];
	_probGeno[index(GT::CT)][index(GT::TT)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::GT)][index(GT::TT)] = tmp2 * thesePi[Base::G];

	//3: case aa/bb
	tmp = thesePhi[DP::aa_bb] / _K[DP::aa_bb];
	tmp2 = tmp * thesePi[Base::A];
	_probGeno[index(GT::AA)][index(GT::CC)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::AA)][index(GT::GG)] = tmp2 * thesePi[Base::G];
	_probGeno[index(GT::AA)][index(GT::TT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::C];
	_probGeno[index(GT::CC)][index(GT::AA)] = tmp2 * thesePi[Base::A];
	_probGeno[index(GT::CC)][index(GT::GG)] = tmp2 * thesePi[Base::G];
	_probGeno[index(GT::CC)][index(GT::TT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::G];
	_probGeno[index(GT::GG)][index(GT::AA)] = tmp2 * thesePi[Base::A];
	_probGeno[index(GT::GG)][index(GT::CC)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::GG)][index(GT::TT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::T];
	_probGeno[index(GT::TT)][index(GT::AA)] = tmp2 * thesePi[Base::A];
	_probGeno[index(GT::TT)][index(GT::CC)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::TT)][index(GT::GG)] = tmp2 * thesePi[Base::G];

	//4: case ab/ab
	tmp = thesePhi[DP::ab_ab] / _K[DP::ab_ab];
	tmp2 = tmp * thesePi[Base::A];
	_probGeno[index(GT::AC)][index(GT::AC)] = tmp2 * thesePi[Base::C];
	_probGeno[index(GT::AG)][index(GT::AG)] = tmp2 * thesePi[Base::G];
	_probGeno[index(GT::AT)][index(GT::AT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::C];
	_probGeno[index(GT::CG)][index(GT::CG)] = tmp2 * thesePi[Base::G];
	_probGeno[index(GT::CT)][index(GT::CT)] = tmp2 * thesePi[Base::T];
	_probGeno[index(GT::GT)][index(GT::GT)] = tmp * thesePi[Base::G] * thesePi[Base::T];

	//5: case ab/ac
	tmp = thesePhi[DP::ab_ac] / _K[DP::ab_ac];
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::G];
	_probGeno[index(GT::AC)][index(GT::AG)] = tmp2;
	_probGeno[index(GT::AC)][index(GT::CG)] = tmp2;
	_probGeno[index(GT::AG)][index(GT::AC)] = tmp2;
	_probGeno[index(GT::AG)][index(GT::CG)] = tmp2;
	_probGeno[index(GT::CG)][index(GT::AC)] = tmp2;
	_probGeno[index(GT::CG)][index(GT::AG)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::T];
	_probGeno[index(GT::AC)][index(GT::AT)] = tmp2;
	_probGeno[index(GT::AC)][index(GT::CT)] = tmp2;
	_probGeno[index(GT::AT)][index(GT::AC)] = tmp2;
	_probGeno[index(GT::AT)][index(GT::CT)] = tmp2;
	_probGeno[index(GT::CT)][index(GT::AC)] = tmp2;
	_probGeno[index(GT::CT)][index(GT::AT)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::G] * thesePi[Base::T];
	_probGeno[index(GT::AG)][index(GT::GT)] = tmp2;
	_probGeno[index(GT::AG)][index(GT::AT)] = tmp2;
	_probGeno[index(GT::AT)][index(GT::AG)] = tmp2;
	_probGeno[index(GT::AT)][index(GT::GT)] = tmp2;
	_probGeno[index(GT::GT)][index(GT::AG)] = tmp2;
	_probGeno[index(GT::GT)][index(GT::AT)] = tmp2;
	tmp2 = tmp * thesePi[Base::C] * thesePi[Base::G] * thesePi[Base::T];
	_probGeno[index(GT::CG)][index(GT::CT)] = tmp2;
	_probGeno[index(GT::CG)][index(GT::GT)] = tmp2;
	_probGeno[index(GT::CT)][index(GT::CG)] = tmp2;
	_probGeno[index(GT::CT)][index(GT::GT)] = tmp2;
	_probGeno[index(GT::GT)][index(GT::CG)] = tmp2;
	_probGeno[index(GT::GT)][index(GT::CT)] = tmp2;

	//6: case aa/bc
	tmp = thesePhi[DP::aa_bc] / _K[DP::aa_bc];
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::G];
	_probGeno[index(GT::AA)][index(GT::CG)] = tmp2;
	_probGeno[index(GT::CC)][index(GT::AG)] = tmp2;
	_probGeno[index(GT::GG)][index(GT::AC)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::T];
	_probGeno[index(GT::AA)][index(GT::CT)] = tmp2;
	_probGeno[index(GT::CC)][index(GT::AT)] = tmp2;
	_probGeno[index(GT::TT)][index(GT::AC)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::G] * thesePi[Base::T];
	_probGeno[index(GT::AA)][index(GT::GT)] = tmp2;
	_probGeno[index(GT::GG)][index(GT::AT)] = tmp2;
	_probGeno[index(GT::TT)][index(GT::AG)] = tmp2;
	tmp2 = tmp * thesePi[Base::C] * thesePi[Base::G] * thesePi[Base::T];
	_probGeno[index(GT::CC)][index(GT::GT)] = tmp2;
	_probGeno[index(GT::GG)][index(GT::CT)] = tmp2;
	_probGeno[index(GT::TT)][index(GT::CG)] = tmp2;

	//7: case ab/cc
	tmp = thesePhi[DP::ab_cc] / _K[DP::ab_cc];
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::G];
	_probGeno[index(GT::AC)][index(GT::GG)] = tmp2;
	_probGeno[index(GT::AG)][index(GT::CC)] = tmp2;
	_probGeno[index(GT::CG)][index(GT::AA)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::T];
	_probGeno[index(GT::AC)][index(GT::TT)] = tmp2;
	_probGeno[index(GT::AT)][index(GT::CC)] = tmp2;
	_probGeno[index(GT::CT)][index(GT::AA)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::G] * thesePi[Base::T];
	_probGeno[index(GT::AG)][index(GT::TT)] = tmp2;
	_probGeno[index(GT::AT)][index(GT::GG)] = tmp2;
	_probGeno[index(GT::GT)][index(GT::AA)] = tmp2;
	tmp2 = tmp * thesePi[Base::C] * thesePi[Base::G] * thesePi[Base::T];
	_probGeno[index(GT::CG)][index(GT::TT)] = tmp2;
	_probGeno[index(GT::CT)][index(GT::GG)] = tmp2;
	_probGeno[index(GT::GT)][index(GT::CC)] = tmp2;

	//8: case ab/cd
	tmp = thesePhi[DP::ab_cd] / _K[DP::ab_cd];
	_probGeno[index(GT::AC)][index(GT::GT)] = tmp;
	_probGeno[index(GT::AG)][index(GT::CT)] = tmp;
	_probGeno[index(GT::AT)][index(GT::CG)] = tmp;
	_probGeno[index(GT::CG)][index(GT::AT)] = tmp;
	_probGeno[index(GT::CT)][index(GT::AG)] = tmp;
	_probGeno[index(GT::GT)][index(GT::AC)] = tmp;
}

bool TEMforDistanceEstimation::estimatePhiWithEM(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2){
	//prepare estimates
	using coretools::index;
	using genometools::Genotype;
	logfile().listFlush("Estimating initial base frequencies pi ...");
	_guessPi(genoQual1, genoQual2);
	logfile().done();
	logfile().conclude("Initial pi are A=", _pi[Base::A], ", C=", _pi[Base::C], ", G=", _pi[Base::G], " and T=", _pi[Base::T], ".");
	logfile().listFlush("Estimating initial genotype classes phi ...");
	_guessPhi(genoQual1, genoQual2);
	logfile().done();
	logfile().conclude("Initial phi are " + coretools::str::concatenateString(_phi, ", ") + ".");

	//variables
	double old_LL, LL = 0.0;
	double LL_diff;

	//now run EM
	logfile().startIndent("Estimating phi using an EM algorithm:");
	for(int iter=0; iter<_maxIterations; ++iter){
		logfile().listFlush("Running EM iteration ", iter+1, " ...");
		//save old LL
		old_LL = LL;
		LL = 0.0;

		//calculate P(g|phi, pi)
		_fill_K(_pi);
		_fill_P_g_given_phi_pi(_phi, _pi);

		//set P_G to zero
		for(int g1 = 0; g1<10; ++g1){
			for(int g2 = 0; g2<10; ++g2){
				_P_G[g1][g2] = 0.0;
			}
		}

		//loop across loci to calculate P_G
		auto it1 = genoQual1.begin();
		auto it2 = genoQual2.begin();
		for(; it1 != genoQual1.end(); ++it1, ++it2){
			//calculate P_G per site
			double sum = 0.0;
			//for(int g1 = 0; g1<10; ++g1){
			for (Genotype g1 = Genotype::min; g1 < Genotype::max; ++g1) {
				const auto ig1 = index(g1);
				for (Genotype g2 = Genotype::min; g2 < Genotype::max; ++g2) {
				const auto ig2 = index(g2);
					_P_G_one_site[ig1][ig2] = (coretools::Probability) (*it1)[g1] * (coretools::Probability) (*it2)[g2] * _probGeno[ig1][ig2];
					sum += _P_G_one_site[ig1][ig2];
				}
			}

			//now add to P_G
			for(int g1 = 0; g1<10; ++g1){
				for(int g2 = 0; g2<10; ++g2){
					_P_G[g1][g2] += _P_G_one_site[g1][g2] / sum;
				}
			}
			LL += log(sum);
		}

		//update phi
		_phi.fill(0.0);

		double sum = 0.0;
		for(Genotype g1 = Genotype::min; g1 < Genotype::max; ++g1){
			for(Genotype g2 = Genotype::min; g2 < Genotype::max; ++g2){
				_phi[distancePhi(g1,g2)] += _P_G[index(g1)][index(g2)];
				sum += _P_G[index(g1)][index(g2)];
			}
		}
		normalize(_phi ,sum);

		//update pi
		_pi.fill(0.0);
		for(Genotype g1 = Genotype::min; g1 < Genotype::max; ++g1){
			for(Genotype g2 = Genotype::min; g2 < Genotype::max; ++g2){
				for(Base b = Base::min; b < Base::max; ++b){
					if(_genoToBaseMap()[g1][g2][b]){
						_pi[b] += _P_G[index(g1)][index(g2)];
					}
				}
			}
		}
		normalize(_pi);

		//check if EM converged
		logfile().done();

		if(iter > 0 ){
			LL_diff = LL - old_LL;
			logfile().conclude("LL = ", LL, " (deltaLL = ", LL_diff, ").");
			if(LL_diff < _minDeltaLL){
				logfile().conclude("EM converged, delatLL = ", LL_diff, " < ", _minDeltaLL);
				_distance = _distanceObject->calculateDistance(_phi);
				logfile().conclude("Resulting distance is ", _distance);
				logfile().endIndent();
				return true;
			}
		} else
			logfile().conclude("LL = ", LL, ".");
	}
	logfile().warning("EM reached maximum number of iterations (", _maxIterations, ") without converging!");
	_distance = _distanceObject->calculateDistance(_phi);
	logfile().conclude("Resulting distance is ", _distance);
	logfile().endIndent();
	return false;
}

//----------------------------------------------------
//TDistanceEstimator
//----------------------------------------------------
TDistanceEstimator::TDistanceEstimator(){
	maxterations = 0;
	_epsilon = 0.0;

	//outputname
	_outputName = coretools::instances::parameters().get<std::string>("out", "ATLAS");
	coretools::instances::logfile().list("Writing output files with prefix '" + _outputName + "'. (parameter 'out')");

	user_assert(_GLFs.size() >= 2, "At least two GLF files have to be provided to estimate distances!");
}

//------------------------------------------------------------------
void TDistanceEstimator::run(){
	//open EM object
	TEMforDistanceEstimation EM_object;

	//in windows or whole genome?
	long windowLen = coretools::instances::parameters().get("window", -1L);
	if(windowLen < 0)
		_estimateDistanceGenomeWide(EM_object);
	else
		_estimateDistanceInWindows(EM_object, windowLen);
}

//--------------------------------------------
// Estimation Genome Wide
//--------------------------------------------
void TDistanceEstimator::_estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object){
	logfile().list("Will estimate genetic distances genome wide.");

	//open output file
	std::string filename = _outputName + "_distanceEstimates.txt.gz";
	coretools::TLineWriter out(filename);

	//write header to output file
	out.writeln("individual1\tindividual2\tnumSitesWithData\tfreqA\tfreqC\tfreqG\tfreqT\tfreqaa/aa\tfreqaa/ab\tfreqab/aa\tfreqaa/bb\tfreqab/ab\tfreqab/ac\tfreqaa/bc\tfreqab/cc\tfreqab/cd\tgeneticDist");

	//prepare storage for distance matrix
	std::vector<double> distMatrix;
	auto numGLFs = _GLFs.size();
	distMatrix.resize(numGLFs*numGLFs, 0.);

	//loop over all pairs
	for(size_t g1=0; g1<(numGLFs-1); ++g1){
		for(size_t g2 = g1+1; g2 < numGLFs; ++g2){
			logfile().startIndent("Estimating distance between individuals ", g1+1, " (" + _GLFs.sampleName(g1), ") and ", g2+1, " (", _GLFs.sampleName(g2), "):");

			//write names to file
			out.write(_GLFs.sampleName(g1), "\t", _GLFs.sampleName(g2));

			//run estimation
			_estimateDistanceGenomeWide(EM_object, _GLFs[g1], _GLFs[g2], out);

			//write to matrix
			distMatrix[g1*numGLFs + g2] = EM_object.distance();
			distMatrix[g1*numGLFs + g2] = EM_object.distance();
			logfile().endIndent();
		}
	}

	out.close();

	//open matrix file
	filename = _outputName + "_distanceMatrix.txt";
	coretools::TLineWriter distMatrixFile(filename);

	//write header to matrix file
	distMatrixFile.write("/");
	for(size_t g=0; g<numGLFs; ++g)
		distMatrixFile.write("\t",_GLFs.sampleName(g));
	distMatrixFile.endln();

	//write rows
	for(size_t g1 = 0; g1 < numGLFs; ++g1){
		distMatrixFile.write(_GLFs.sampleName(g1));
		for(size_t g2 = 0; g2 < numGLFs; ++g2)
			distMatrixFile.write("\t",distMatrix[g1*numGLFs + g2]);
		distMatrixFile.endln();
	}

	//close file
	distMatrixFile.close();
}

void TDistanceEstimator::_moveToNextCommonChr(genometools::TGLFReader & g1, genometools::TGLFReader & g2){
	while((!g1.empty() && !g2.empty()) && (g1.refID() != g2.refID() || g1.curChr().isHaploid() || g2.curChr().isHaploid())) {
		//advance the one laging behind
		if(g1.refID() < g2.refID()){
			if(!g1.jumpToNextChr()) return;
		} else {
			if(!g2.jumpToNextChr()) return;
		}
	}
}

void TDistanceEstimator::_advance(genometools::TGLFReader & g1, genometools::TGLFReader & g2){
	//advance
	if(g2.position() == g1.position()){
		//advance both
		g1.popFront();
		g2.popFront();
	} else if(g2.position() < g1.position()){
		g2.popFront();
	} else {
		g1.popFront();
	}

	//make sure we are on same chromosome
	return(_moveToNextCommonChr(g1, g2));
}

void TDistanceEstimator::_readCommonSites(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2, genometools::TGLFReader & g1, genometools::TGLFReader & g2){
	//parse GLFs. Only keep sites where both individuals have data!

	//rewind GLFs
	g1.rewind();
	g2.rewind();

	//if not both are good at least one file reach end. So we are done!
	for(; !g1.empty() && !g2.empty(); _advance(g1, g2)) {
		if(g2.position() == g1.position()){
			//add data
			genoQual1.push_back(g1.front().likelihoods);
			genoQual2.push_back(g2.front().likelihoods);
		}
	}
}

void TDistanceEstimator::_estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object, genometools::TGLFReader & g1, genometools::TGLFReader & g2, coretools::TLineWriter & out){
	//initialize storage for two windows
	logfile().listFlush("Reading common sites ...");
	GenotypeQualityVector genoQual1, genoQual2;
	_readCommonSites(genoQual1, genoQual2, g1, g2);
	logfile().done();
	logfile().conclude("Read data for ", genoQual1.size(), " sites.");

	//now estimate
	if(genoQual1.size() > 0){
		logfile().startIndent("Estimating genetic distance:");
		EM_object.estimatePhiWithEM(genoQual1, genoQual2);
		_writeDistanceEstimates(out, genoQual1.size(), EM_object);
		logfile().endIndent();
	} else {
		logfile().conclude("Not enough data to estimate distance.");
		_writeDistanceEstimatesNoData(out);
	}

	//clean up memory
	logfile().listFlush("Cleaning up memory ...");
	logfile().done();
}

//--------------------------------------------
// Estimation in windows
//--------------------------------------------
void TDistanceEstimator::_estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, uint32_t windowLen){
	logfile().list("Will estimate genetic distance in windows of length ", windowLen, ".");
	user_assert(windowLen >= 100, "Window size must be at least 100bp!");

	//loop over all pairs
	for(size_t g1=0; g1<(_GLFs.size() - 1); ++g1){
		for(size_t g2 = g1+1; g2 < _GLFs.size(); ++g2){
			logfile().startIndent("Estimating distance between individuals ", g1+1, " (", _GLFs.fileName(g1), ") and ", g2+1, " (", _GLFs.fileName(g2), "):");

			//output file
			std::string filename = _outputName + "_" + _GLFs.fileName(g1) + "_" + _GLFs.fileName(g2) + "_distanceEstimates.txt.gz";
			logfile().list("Will write estimates to file '" + filename + "'.");

			//rewind GLFs
			_GLFs[g1].rewind();
			_GLFs[g2].rewind();

			//now run estimation
			_estimateDistanceInWindows(EM_object, filename, _GLFs[g1], _GLFs[g2], windowLen);

			logfile().endIndent();
		}
	}
}

void TDistanceEstimator::_estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, std::string filename, genometools::TGLFReader & g1, genometools::TGLFReader & g2, uint32_t windowLen){
	//initialize variables
	bool isGood1 = true;
	bool isGood2 = true;

	//initialize storage for two windows
	//TODO: share across pairs? Do all pairs at once?
	//Todo: store arrays rather than pointers
	GenotypeQualityVector genoQual1, genoQual2;
	genoQual1.resize(windowLen);
	genoQual2.resize(windowLen);

	//open output file
	coretools::TLineWriter out(filename.c_str());

	//write header to output file
	out.writeln("chr\twindowStart\twindowEnd\tnumSitesWithData\tfreqA\tfreqC\tfreqG\tfreqT\tfreq00_00\tfreq00_01\tfreq01_00\tfreq00_11\tfreq01_01\tfreq01_02\tfreq00_12\tfreq01_22\tfreq01_23\tgeneticDist");

	//prepare variables
	genometools::TGenomeWindow window;

	int numSitesWithData = 100;

	//parse GLFs in windows
	logfile().startIndent("Will estimate distance in windows of size ", windowLen, ":");
	while(!g1.empty() && !g2.empty()){
		//move to new chromosome
		window.move(g1.refID(), 0, windowLen);
		const genometools::TChromosome& curChr = g1.curChr();

		logfile().startNumbering("Chromosome ", curChr.name(), ":");

		//parse all windows of chromosome
		while(window < curChr.to() && !g1.empty() && !g2.empty()){
			logfile().number("Window [", window, ")");
			logfile().addIndent();

			//read data
			isGood1 = g1.readNextWindow(genoQual1, window);
			if(isGood1 || g1.empty()){
				isGood2 = g2.readNextWindow(genoQual2, window);
				if(isGood2 || g2.empty()){
					//estimate distance
					EM_object.estimatePhiWithEM(genoQual1, genoQual2);

					//write to file
					_writeDistanceEstimates(out, curChr, window, numSitesWithData, EM_object);


				} else _writeDistanceEstimatesNoData(out, curChr, window);
			} else _writeDistanceEstimatesNoData(out, curChr, window);

			//move window
			window += windowLen;
			logfile().endIndent();
		}
		logfile().endNumbering();
	}

	logfile().endIndent();
}

//--------------------------------------------
// Writing estimates
//--------------------------------------------
void TDistanceEstimator::_writeDistanceEstimates(coretools::TLineWriter & out, const genometools::TChromosome& Chr, genometools::TGenomeWindow& Window, uint32_t numsitesWithData, TEMforDistanceEstimation & EM_object){
	out.write(Chr.name(),"\t",Window.from().position() + 1,"\t",Window.to().position()); //internal position is zero-based
	_writeDistanceEstimates(out, numsitesWithData, EM_object);
}

void TDistanceEstimator::_writeDistanceEstimates(coretools::TLineWriter & out, int numsitesWithData, TEMforDistanceEstimation & EM_object){
	using coretools::index;
	out.write("\t",numsitesWithData);
	//write pi
	for(Base b = Base::min; b < Base::max; ++b){
		out.write("\t",EM_object.pi()[b]);
	}
	//write phi
	for(DistancePhi p = DistancePhi::min; p < DistancePhi::max; ++p){
		out.write("\t",EM_object.phi()[p]);
	}
	//write distance
	out.write("\t",EM_object.distance());
	out.endln();
}

void TDistanceEstimator::_writeDistanceEstimatesNoData(coretools::TLineWriter & out, const genometools::TChromosome& Chr, genometools::TGenomeWindow& Window){
	out.write(Chr.name(),"\t",Window.from().position() + 1,"\t",Window.to().position(),"\t"); //internal position is zero-based
	_writeDistanceEstimatesNoData(out);
}

void TDistanceEstimator::_writeDistanceEstimatesNoData(coretools::TLineWriter & out){
	out.write("\t0");
	for(int i=0; i<14; ++i)
		out.write("\t-");
	out.endln();
}

} //end namespace
