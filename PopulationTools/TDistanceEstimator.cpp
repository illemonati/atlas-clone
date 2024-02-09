/*
 * TDistanceCalculator.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#include "TDistanceEstimator.h"
#include "coretools/Files/gzstream.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/concatenateString.h"
#include "coretools/algorithms.h"

namespace PopulationTools{

using GenotypeLikelihoods::TBaseData;
using genometools::Genotype;
using genometools::Base;
using coretools::Probability;

//----------------------------------------------------
//TGenocombinationToBaseMap
//----------------------------------------------------
TGenocombinationToBaseMap::TGenocombinationToBaseMap(){
	using coretools::index;
	using genometools::first;
	using genometools::second;
	for(Genotype g1 = Genotype::min; g1 < Genotype::max; ++g1){
		for(Genotype g2 = Genotype::min; g2 < Genotype::max; ++g2){
			for(Base b = Base::min; b < Base::max; ++b){
				genotypeCombinationHasBase[index(g1)][index(g2)][index(b)] = false;
			}
			genotypeCombinationHasBase[index(g1)][index(g2)][index(first(g1))]  = true;
			genotypeCombinationHasBase[index(g1)][index(g2)][index(second(g1))] = true;
			genotypeCombinationHasBase[index(g1)][index(g2)][index(first(g2))]  = true;
			genotypeCombinationHasBase[index(g1)][index(g2)][index(second(g2))] = true;
		}
	}
};

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
};

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
};

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
};

double TDistanceEuclidian::calculateDistance(const TDistanceData & phi){
	return sqrt(TDistance::calculateDistance(phi));
};

TDistanceUser::TDistanceUser(std::vector<double> vec){
	using coretools::index;
	for(DistancePhi d = DistancePhi::min; d < DistancePhi::max; ++d)
		_distanceWeight[d] = vec[index(d)];
};

//----------------------------------------------------
//TDistanceEstimate
//----------------------------------------------------
TEMforDistanceEstimation::TEMforDistanceEstimation(){
	//prepare storage
	phi.fill(0.0);
	LL = 0.0;
	old_LL = 0.0;
	distance = -1.0;

	//read EM parameters
	using namespace coretools::instances;
	logfile().startIndent("Parameters of EM algorithm:");
	maxNumEMIterations = parameters().get<int>("iterations", 100);
	logfile().list("Will run up to ", maxNumEMIterations, " iterations.");
	epsilonForEM = parameters().get("maxEps", 0.000001);
	logfile().list("Will run EM until deltaLL < ", epsilonForEM, ".");
	logfile().endIndent();

	//set how to calculate distances
//	distanceWeight = new double[9];
	if(parameters().exists("distWeights")){
		logfile().list("Using user-provided distance weights.");
		const auto vec = parameters().get<std::vector<double>>("distWeights");
		if(vec.size() != 9)
			UERROR("Wrong number of distance weights! Required are nine values for aa/aa, aa/ab, ab/aa, aa/bb, ab/ab, ab/ac, aa/bc, ab/cc, ab/cd");

		distanceObject = std::make_unique<TDistanceUser>(vec);

	} else {
		std::string distType = parameters().get<std::string>("distType", "squaredDiff");
		logfile().list("Using distance type '" + distType + "'.");
		if(distType == "probMismatch"){
			distanceObject = std::make_unique<TDistanceProbMismatch>();
		} else if(distType == "squaredDiff"){
			distanceObject = std::make_unique<TDistance>();
		} else if(distType == "euclidian"){
			distanceObject = std::make_unique<TDistanceEuclidian>();
		} else
			UERROR("Unknown distance type '", distType, "'! Use either squaredDiff, euclidian, or probMismatch.");
	}
	logfile().conclude("Using distance weights " + coretools::str::concatenateString(distanceObject->weights(), ", ") + ".");

};

void TEMforDistanceEstimation::guessPi(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2){
	//check sizes are equal
	using coretools::index;
	if(genoQual1.size() != genoQual2.size())
		UERROR("Provided genotype quality vectors are of different size in TEMforDistanceEstimation::guessPi!");

	//just estimate pi as average posterior probability
	pi.fill(0.0);

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
			pi[first(g)] += tmp;
			pi[second(g)] += tmp;
			tmp = (Probability) (*it2)[g] / sum2;
			pi[first(g)] += tmp;
			pi[second(g)] += tmp;
		}
	}

	//normalize
	normalize(pi);
}

void TEMforDistanceEstimation::guessPhi(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2){
	using coretools::index;
	//check sizes are equal
	if(genoQual1.size() != genoQual2.size())
		UERROR("Provided genotype quality vectors are of different size in TEMforDistanceEstimation::guessPhi!");

	//set to zero
	phi.fill(0.0);

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
				phi[distancePhi(g1, g2)] += tmp * ((Probability) (*it2)[g2] / sum2);
			}
		}
	}

	//normalize
	normalize(phi);
}

void TEMforDistanceEstimation::fill_K(TBaseData & thesePi){
	using genometools::Base;
	using DP = DistancePhi;
	//normalizing constant for each phi class
	//case of one base
	K[DP::aa_aa] = 1.0;

	//cases of two bases
	K[DP::ab_ab] = thesePi[Base::A] * thesePi[Base::C]
         + thesePi[Base::A] * thesePi[Base::G]
         + thesePi[Base::A] * thesePi[Base::T]
		 + thesePi[Base::C] * thesePi[Base::G]
		 + thesePi[Base::C] * thesePi[Base::T]
		 + thesePi[Base::G] * thesePi[Base::T];
	K[DP::aa_ab] = 2.0 * K[DP::ab_ab]; //account for AC vs CA
	K[DP::ab_aa] = K[DP::aa_ab];
	K[DP::aa_bb] = K[DP::aa_ab];

	//cases of three bases
	K[DP::aa_bc] = thesePi[Base::C] * thesePi[Base::G] * thesePi[Base::T]
		 + thesePi[Base::A] * thesePi[Base::G] * thesePi[Base::T]
		 + thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::T]
		 + thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::G];
	K[DP::aa_bc] = 3.0 * K[DP::aa_bc]; //account for ways to distribute

	K[DP::ab_cc] = K[DP::aa_bc];
	K[DP::ab_ac] = 2.0 * K[DP::aa_bc]; //twice as many cases than other cases with three bases!

	//case of four bases: each of the 6 cases is equally likely
	//Note: product of pi's cancels out when calculating P_g_given_phi_pi
	K[DP::ab_cd] = 6.0;
};

void TEMforDistanceEstimation::fill_P_g_given_phi_pi(const TDistanceData & thesePhi, TBaseData & thesePi){
	using genometools::Base;
	using coretools::index;
	using DP = DistancePhi;
	using GT = genometools::Genotype;
	//0: case aa/aa (K[0]=1)
	probGeno[index(GT::AA)][index(GT::AA)] = thesePhi[DP::aa_aa] * thesePi[Base::A];
	probGeno[index(GT::CC)][index(GT::CC)] = thesePhi[DP::aa_aa] * thesePi[Base::C];
	probGeno[index(GT::GG)][index(GT::GG)] = thesePhi[DP::aa_aa] * thesePi[Base::G];
	probGeno[index(GT::TT)][index(GT::TT)] = thesePhi[DP::aa_aa] * thesePi[Base::T];

	//1: cases aa/ab
	double tmp = thesePhi[DP::aa_ab] / K[DP::aa_aa];
	double tmp2 = tmp * thesePi[Base::A];
	probGeno[index(GT::AA)][index(GT::AC)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::AA)][index(GT::AG)] = tmp2 * thesePi[Base::G];
	probGeno[index(GT::AA)][index(GT::AT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::C];
	probGeno[index(GT::CC)][index(GT::AC)] = tmp2 * thesePi[Base::A];
	probGeno[index(GT::CC)][index(GT::CG)] = tmp2 * thesePi[Base::G];
	probGeno[index(GT::CC)][index(GT::CT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::G];
	probGeno[index(GT::GG)][index(GT::AG)] = tmp2 * thesePi[Base::A];
	probGeno[index(GT::GG)][index(GT::CG)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::GG)][index(GT::GT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::T];
	probGeno[index(GT::TT)][index(GT::AT)] = tmp2 * thesePi[Base::A];
	probGeno[index(GT::TT)][index(GT::CT)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::TT)][index(GT::GT)] = tmp2 * thesePi[Base::G];

	//2: case ab/aa
	tmp = thesePhi[DP::ab_aa] / K[DP::ab_aa];
	tmp2 = tmp * thesePi[Base::A];
	probGeno[index(GT::AC)][index(GT::AA)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::AG)][index(GT::AA)] = tmp2 * thesePi[Base::G];
	probGeno[index(GT::AT)][index(GT::AA)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::C];
	probGeno[index(GT::AC)][index(GT::CC)] = tmp2 * thesePi[Base::A];
	probGeno[index(GT::CG)][index(GT::CC)] = tmp2 * thesePi[Base::G];
	probGeno[index(GT::CT)][index(GT::CC)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::G];
	probGeno[index(GT::AG)][index(GT::GG)] = tmp2 * thesePi[Base::A];
	probGeno[index(GT::CG)][index(GT::GG)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::GT)][index(GT::GG)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::T];
	probGeno[index(GT::AT)][index(GT::TT)] = tmp2 * thesePi[Base::A];
	probGeno[index(GT::CT)][index(GT::TT)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::GT)][index(GT::TT)] = tmp2 * thesePi[Base::G];

	//3: case aa/bb
	tmp = thesePhi[DP::aa_bb] / K[DP::aa_bb];
	tmp2 = tmp * thesePi[Base::A];
	probGeno[index(GT::AA)][index(GT::CC)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::AA)][index(GT::GG)] = tmp2 * thesePi[Base::G];
	probGeno[index(GT::AA)][index(GT::TT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::C];
	probGeno[index(GT::CC)][index(GT::AA)] = tmp2 * thesePi[Base::A];
	probGeno[index(GT::CC)][index(GT::GG)] = tmp2 * thesePi[Base::G];
	probGeno[index(GT::CC)][index(GT::TT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::G];
	probGeno[index(GT::GG)][index(GT::AA)] = tmp2 * thesePi[Base::A];
	probGeno[index(GT::GG)][index(GT::CC)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::GG)][index(GT::TT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::T];
	probGeno[index(GT::TT)][index(GT::AA)] = tmp2 * thesePi[Base::A];
	probGeno[index(GT::TT)][index(GT::CC)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::TT)][index(GT::GG)] = tmp2 * thesePi[Base::G];

	//4: case ab/ab
	tmp = thesePhi[DP::ab_ab] / K[DP::ab_ab];
	tmp2 = tmp * thesePi[Base::A];
	probGeno[index(GT::AC)][index(GT::AC)] = tmp2 * thesePi[Base::C];
	probGeno[index(GT::AG)][index(GT::AG)] = tmp2 * thesePi[Base::G];
	probGeno[index(GT::AT)][index(GT::AT)] = tmp2 * thesePi[Base::T];
	tmp2 = tmp * thesePi[Base::C];
	probGeno[index(GT::CG)][index(GT::CG)] = tmp2 * thesePi[Base::G];
	probGeno[index(GT::CT)][index(GT::CT)] = tmp2 * thesePi[Base::T];
	probGeno[index(GT::GT)][index(GT::GT)] = tmp * thesePi[Base::G] * thesePi[Base::T];

	//5: case ab/ac
	tmp = thesePhi[DP::ab_ac] / K[DP::ab_ac];
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::G];
	probGeno[index(GT::AC)][index(GT::AG)] = tmp2;
	probGeno[index(GT::AC)][index(GT::CG)] = tmp2;
	probGeno[index(GT::AG)][index(GT::AC)] = tmp2;
	probGeno[index(GT::AG)][index(GT::CG)] = tmp2;
	probGeno[index(GT::CG)][index(GT::AC)] = tmp2;
	probGeno[index(GT::CG)][index(GT::AG)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::T];
	probGeno[index(GT::AC)][index(GT::AT)] = tmp2;
	probGeno[index(GT::AC)][index(GT::CT)] = tmp2;
	probGeno[index(GT::AT)][index(GT::AC)] = tmp2;
	probGeno[index(GT::AT)][index(GT::CT)] = tmp2;
	probGeno[index(GT::CT)][index(GT::AC)] = tmp2;
	probGeno[index(GT::CT)][index(GT::AT)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::G] * thesePi[Base::T];
	probGeno[index(GT::AG)][index(GT::GT)] = tmp2;
	probGeno[index(GT::AG)][index(GT::AT)] = tmp2;
	probGeno[index(GT::AT)][index(GT::AG)] = tmp2;
	probGeno[index(GT::AT)][index(GT::GT)] = tmp2;
	probGeno[index(GT::GT)][index(GT::AG)] = tmp2;
	probGeno[index(GT::GT)][index(GT::AT)] = tmp2;
	tmp2 = tmp * thesePi[Base::C] * thesePi[Base::G] * thesePi[Base::T];
	probGeno[index(GT::CG)][index(GT::CT)] = tmp2;
	probGeno[index(GT::CG)][index(GT::GT)] = tmp2;
	probGeno[index(GT::CT)][index(GT::CG)] = tmp2;
	probGeno[index(GT::CT)][index(GT::GT)] = tmp2;
	probGeno[index(GT::GT)][index(GT::CG)] = tmp2;
	probGeno[index(GT::GT)][index(GT::CT)] = tmp2;

	//6: case aa/bc
	tmp = thesePhi[DP::aa_bc] / K[DP::aa_bc];
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::G];
	probGeno[index(GT::AA)][index(GT::CG)] = tmp2;
	probGeno[index(GT::CC)][index(GT::AG)] = tmp2;
	probGeno[index(GT::GG)][index(GT::AC)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::T];
	probGeno[index(GT::AA)][index(GT::CT)] = tmp2;
	probGeno[index(GT::CC)][index(GT::AT)] = tmp2;
	probGeno[index(GT::TT)][index(GT::AC)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::G] * thesePi[Base::T];
	probGeno[index(GT::AA)][index(GT::GT)] = tmp2;
	probGeno[index(GT::GG)][index(GT::AT)] = tmp2;
	probGeno[index(GT::TT)][index(GT::AG)] = tmp2;
	tmp2 = tmp * thesePi[Base::C] * thesePi[Base::G] * thesePi[Base::T];
	probGeno[index(GT::CC)][index(GT::GT)] = tmp2;
	probGeno[index(GT::GG)][index(GT::CT)] = tmp2;
	probGeno[index(GT::TT)][index(GT::CG)] = tmp2;

	//7: case ab/cc
	tmp = thesePhi[DP::ab_cc] / K[DP::ab_cc];
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::G];
	probGeno[index(GT::AC)][index(GT::GG)] = tmp2;
	probGeno[index(GT::AG)][index(GT::CC)] = tmp2;
	probGeno[index(GT::CG)][index(GT::AA)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::C] * thesePi[Base::T];
	probGeno[index(GT::AC)][index(GT::TT)] = tmp2;
	probGeno[index(GT::AT)][index(GT::CC)] = tmp2;
	probGeno[index(GT::CT)][index(GT::AA)] = tmp2;
	tmp2 = tmp * thesePi[Base::A] * thesePi[Base::G] * thesePi[Base::T];
	probGeno[index(GT::AG)][index(GT::TT)] = tmp2;
	probGeno[index(GT::AT)][index(GT::GG)] = tmp2;
	probGeno[index(GT::GT)][index(GT::AA)] = tmp2;
	tmp2 = tmp * thesePi[Base::C] * thesePi[Base::G] * thesePi[Base::T];
	probGeno[index(GT::CG)][index(GT::TT)] = tmp2;
	probGeno[index(GT::CT)][index(GT::GG)] = tmp2;
	probGeno[index(GT::GT)][index(GT::CC)] = tmp2;

	//8: case ab/cd
	tmp = thesePhi[DP::ab_cd] / K[DP::ab_cd];
	probGeno[index(GT::AC)][index(GT::GT)] = tmp;
	probGeno[index(GT::AG)][index(GT::CT)] = tmp;
	probGeno[index(GT::AT)][index(GT::CG)] = tmp;
	probGeno[index(GT::CG)][index(GT::AT)] = tmp;
	probGeno[index(GT::CT)][index(GT::AG)] = tmp;
	probGeno[index(GT::GT)][index(GT::AC)] = tmp;
};

bool TEMforDistanceEstimation::estimatePhiWithEM(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2){
	//prepare estimates
	using coretools::index;
	using genometools::Genotype;
	using namespace coretools::instances;
	logfile().listFlush("Estimating initial base frequencies pi ...");
	guessPi(genoQual1, genoQual2);
	logfile().done();
	logfile().conclude("Initial pi are A=", pi[Base::A], ", C=", pi[Base::C], ", G=", pi[Base::G], " and T=", pi[Base::T], ".");
	logfile().listFlush("Estimating initial genotype classes phi ...");
	guessPhi(genoQual1, genoQual2);
	logfile().done();
	logfile().conclude("Initial phi are " + coretools::str::concatenateString(phi, ", ") + ".");

	//variables
	double old_LL, LL = 0.0;
	double LL_diff;

	//now run EM
	logfile().startIndent("Estimating phi using an EM algorithm:");
	for(int iter=0; iter<maxNumEMIterations; ++iter){
		logfile().listFlush("Running EM iteration ", iter+1, " ...");
		//save old LL
		old_LL = LL;
		LL = 0.0;

		//calculate P(g|phi, pi)
		fill_K(pi);
		fill_P_g_given_phi_pi(phi, pi);

		//set P_G to zero
		for(int g1 = 0; g1<10; ++g1){
			for(int g2 = 0; g2<10; ++g2){
				P_G[g1][g2] = 0.0;
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
					P_G_one_site[ig1][ig2] = (coretools::Probability) (*it1)[g1] * (coretools::Probability) (*it2)[g2] * probGeno[ig1][ig2];

					//std::cout << g1 << "/" << g2 << ": " <<  phredToLik[genoQual1[s][g1]] << " * " << phredToLik[genoQual2[s][g2]] << " * " << probGeno[g1][g2] << " = " << P_G_one_site[g1][g2] << std::endl;

					sum += P_G_one_site[ig1][ig2];
				}
			}

			//now add to P_G
			for(int g1 = 0; g1<10; ++g1){
				for(int g2 = 0; g2<10; ++g2){
					P_G[g1][g2] += P_G_one_site[g1][g2] / sum;
				}
			}

			//add to LL
			LL += log(sum);
			//std::cout << s << "\t" << sum << " -> " << LL << "\n";
		}

		//update phi
		phi.fill(0.0);

		double sum = 0.0;
		for(Genotype g1 = Genotype::min; g1 < Genotype::max; ++g1){
			for(Genotype g2 = Genotype::min; g2 < Genotype::max; ++g2){
				phi[distancePhi(g1,g2)] += P_G[index(g1)][index(g2)];
				sum += P_G[index(g1)][index(g2)];
			}
		}
		normalize(phi ,sum);

		//update pi
		pi.fill(0.0);
		for(Genotype g1 = Genotype::min; g1 < Genotype::max; ++g1){
			for(Genotype g2 = Genotype::min; g2 < Genotype::max; ++g2){
				for(Base b = Base::min; b < Base::max; ++b){
					if(genoToBaseMap(g1, g2, b)){
						pi[b] += P_G[index(g1)][index(g2)];
					}
				}
			}
		}
		normalize(pi);

		//check if EM converged
		logfile().done();
		//UERROR("done!";
		if(iter > 0 ){
			LL_diff = LL - old_LL;
			logfile().conclude("LL = ", LL, " (deltaLL = ", LL_diff, ").");
			if(LL_diff < epsilonForEM){
				logfile().conclude("EM converged, delatLL = ", LL_diff, " < ", epsilonForEM);
				distance = distanceObject->calculateDistance(phi);
				logfile().conclude("Resulting distance is ", distance);
				logfile().endIndent();
				return true;
			}
		} else
			logfile().conclude("LL = ", LL, ".");
	}
	logfile().warning("EM reached maximum number of iterations (", maxNumEMIterations, ") without converging!");
	distance = distanceObject->calculateDistance(phi);
	logfile().conclude("Resulting distance is ", distance);
	logfile().endIndent();
	return false;
};

//----------------------------------------------------
//TDistanceEstimator
//----------------------------------------------------
TDistanceEstimator::TDistanceEstimator(){
	maxNumEMIterations = 0;
	epsilonForEM = 0.0;
	numGLFs = 0;

	//outputname
	outputName = coretools::instances::parameters().get<std::string>("out", "ATLAS");
	coretools::instances::logfile().list("Writing output files with prefix '" + outputName + "'. (parameter 'out')");
}

void TDistanceEstimator::openGLF(){
	using namespace coretools::instances;
	parameters().fill("glf", GLFNames);
	numGLFs = GLFNames.size();
	if(numGLFs < 2)
		UERROR("At least two GLF files have to be provided to estimate distances!");

	//open files
	logfile().startIndent("Opening GLF files:");
	glfs.reserve(numGLFs);
	for (const auto& name: GLFNames) {
		logfile().listFlush("Opening GLF '", name, "' ...");
		glfs.emplace_back(name);
		logfile().done();
	}
	logfile().endIndent();
}

//------------------------------------------------------------------
void TDistanceEstimator::run(){
	//open all GLF files specified
	openGLF();

	//open EM object
	TEMforDistanceEstimation EM_object;

	//in windows or whole genome?
	long windowLen = coretools::instances::parameters().get("window", -1L);
	if(windowLen < 0)
		estimateDistanceGenomeWide(EM_object);
	else
		estimateDistanceInWindows(EM_object, windowLen);
}

//--------------------------------------------
// Estimation Genome Wide
//--------------------------------------------
void TDistanceEstimator::estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object){
	using namespace coretools::instances;
	logfile().list("Will estimate genetic distances genome wide.");

	//open output file
	std::string filename = outputName + "_distanceEstimates.txt.gz";
	gz::ogzstream out(filename.c_str());
	if(!out)
		UERROR("Failed to open output file '", filename, "'!");

	//write header to output file
	out << "individual1\tindividual2\tnumSitesWithData\tfreqA\tfreqC\tfreqG\tfreqT\tfreqaa/aa\tfreqaa/ab\tfreqab/aa\tfreqaa/bb\tfreqab/ab\tfreqab/ac\tfreqaa/bc\tfreqab/cc\tfreqab/cd\tgeneticDist\n";

	//prepare storage for distance matrix
	std::vector<double> distMatrix;
	distMatrix.resize(numGLFs*numGLFs, 0.);

	//loop over all pairs
	for(int g1=0; g1<(numGLFs-1); ++g1){
		for(int g2 = g1+1; g2 < numGLFs; ++g2){
			logfile().startIndent("Estimating distance between individuals ", g1+1, " (" + GLFNames[g1], ") and ", g2+1, " (", GLFNames[g2], "):");

			//write names to file
			out << GLFNames[g1] << "\t" << GLFNames[g2];

			//run estimation
			estimateDistanceGenomeWide(EM_object, glfs[g1], glfs[g2], out);

			//write to matrix
			distMatrix[g1*numGLFs + g2] = EM_object.distance;
			distMatrix[g1*numGLFs + g2] = EM_object.distance;
			logfile().endIndent();
		}
	}

	out.close();

	//open matrix file
	filename = outputName + "_distanceMatrix.txt";
	std::ofstream distMatrixFile(filename.c_str());
	if(!distMatrixFile)
		UERROR("Failed to open output file '", filename, "'!");

	//write header to matrix file
	distMatrixFile << "/";
	for(int g=0; g<numGLFs; ++g)
		distMatrixFile << "\t" << GLFNames[g];
	distMatrixFile << "\n";

	//write rows
	for(int g1 = 0; g1 < numGLFs; ++g1){
		distMatrixFile << GLFNames[g1];
		for(int g2 = 0; g2 < numGLFs; ++g2)
			distMatrixFile << "\t" << distMatrix[g1*numGLFs + g2];
		distMatrixFile << "\n";
	}

	//close file
	distMatrixFile.close();
};

bool TDistanceEstimator::moveToNextCommonChr(GLF::TGlfReader & g1, GLF::TGlfReader & g2){
	while(g1.refId() != g2.refId() || g1.chrIsHaploid() || g2.chrIsHaploid()){
		//advance the one laging behind
		if(g1.refId() < g2.refId()){
			if(!g1.jumpToNextChr()) return false;
		} else {
			if(!g2.jumpToNextChr()) return false;
		}
	}

	//check names
	if(g1.chr() != g2.chr()){
		UERROR("Chromosome names differ in files ", g1.name(), "' and '", g2.name(), "': '", g1.chr(), "' != '", g2.chr() + "'!");
	}
	if(g1.chrLength() != g2.chrLength()){
		UERROR("Chromosome lengths differ in files ", g1.name(), "' and '", g2.name(), "': '", g1.chrLength(), "' != '", g2.chrLength(), "'!");
	}

	return true;
};

bool TDistanceEstimator::advance(GLF::TGlfReader & g1, GLF::TGlfReader & g2){
	//advance
	if(g2.position() == g1.position()){
		//advance both
		if(!g1.readNext()) return false;
		if(!g2.readNext()) return false;
	} else if(g2.position() < g1.position()){
		//advance g2
		if(!g2.readNext()) return false;
	} else {
		//advance g1
		if(!g1.readNext()) return false;
	}

	//make sure we are on same chromosome
	return(moveToNextCommonChr(g1, g2));
};

void TDistanceEstimator::readCommonSites(GenotypeQualityVector & genoQual1, GenotypeQualityVector & genoQual2, GLF::TGlfReader & g1, GLF::TGlfReader & g2){
	//parse GLFs. Only keep sites where both individuals have data!

	//rewind GLFs
	g1.rewind();
	g2.rewind();

	//if not both are good at least one file reach end. So we are done!
	while(advance(g1, g2)){
		if(g2.position() == g1.position()){
			//add data
			genoQual1.push_back(g1.genotypeLikelihoodsGLF());
			genoQual2.push_back(g2.genotypeLikelihoodsGLF());
		}
	}
};

void TDistanceEstimator::estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object, GLF::TGlfReader & g1, GLF::TGlfReader & g2, gz::ogzstream & out){
	//initialize storage for two windows
	using namespace coretools::instances;
	logfile().listFlush("Reading common sites ...");
	GenotypeQualityVector genoQual1, genoQual2;
	readCommonSites(genoQual1, genoQual2, g1, g2);
	logfile().done();
	logfile().conclude("Read data for ", genoQual1.size(), " sites.");

	//now estimate
	if(genoQual1.size() > 0){
		logfile().startIndent("Estimating genetic distance:");
		EM_object.estimatePhiWithEM(genoQual1, genoQual2);
		writeDistanceEstimates(out, genoQual1.size(), EM_object);
		logfile().endIndent();
	} else {
		logfile().conclude("Not enough data to estimate distance.");
		writeDistanceEstimatesNoData(out);
	}

	//clean up memory
	logfile().listFlush("Cleaning up memory ...");
	logfile().done();
};

//--------------------------------------------
// Estimation in windows
//--------------------------------------------
void TDistanceEstimator::estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, uint32_t windowLen){
	using namespace coretools::instances;
	logfile().list("Will estimate genetic distance in windows of length ", windowLen, ".");
	if(windowLen < 100)
		UERROR("Window size must be at least 100bp!");

	//loop over all pairs
	for(int g1=0; g1<(numGLFs-1); ++g1){
		for(int g2 = g1+1; g2 < numGLFs; ++g2){
			logfile().startIndent("Estimating distance between individuals ", g1+1, " (", GLFNames[g1], ") and ", g2+1, " (", GLFNames[g2], "):");

			//output file
			std::string filename = outputName + "_" + GLFNames[g1] + "_" + GLFNames[g2] + "_distanceEstimates.txt.gz";
			logfile().list("Will write estimates to file '" + filename + "'.");

			//rewind GLFs
			glfs[g1].rewind();
			glfs[g2].rewind();

			//now run estimation
			estimateDistanceInWindows(EM_object, filename, glfs[g1], glfs[g2], windowLen);

			logfile().endIndent();
		}
	}
};

void TDistanceEstimator::estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, std::string filename, GLF::TGlfReader & g1, GLF::TGlfReader & g2, uint32_t windowLen){
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
	gz::ogzstream out(filename.c_str());
	if(!out) UERROR("Failed to open file '", filename, "' for writing!");

	//write header to output file
	out << "chr\twindowStart\twindowEnd\tnumSitesWithData\tfreqA\tfreqC\tfreqG\tfreqT\tfreq00_00\tfreq00_01\tfreq01_00\tfreq00_11\tfreq01_01\tfreq01_02\tfreq00_12\tfreq01_22\tfreq01_23\tgeneticDist\n";

	//prepare variables
	uint32_t curRefId;
	std::string curChr;
	long curChrLen;
	long windowStart;
	long windowEnd;

	int numSitesWithData = 100;

	//parse GLFs in windows
	using namespace coretools::instances;
	logfile().startIndent("Will estimate distance in windows of size ", windowLen, ":");
	while(!g1.eof() && !g2.eof()){
		//move to new chromosome
		curRefId = g1.refId();
		curChr = g1.chr();
		curChrLen = g1.chrLength();
		windowStart = 0;
		windowEnd = windowLen;

		logfile().startNumbering("Chromosome " + curChr + ":");

		//parse all windows of chromosome
		while(windowStart < curChrLen && !g1.eof() && !g2.eof()){
			logfile().number("Window [", windowStart, ", ", windowEnd, ")");
			logfile().addIndent();

			//read data
			isGood1 = g1.readNextWindow(genoQual1, curRefId, windowStart, windowEnd);
			if(isGood1 || g1.eof()){
				isGood2 = g2.readNextWindow(genoQual2, curRefId, windowStart, windowEnd);
				if(isGood2 || g2.eof()){
					//estimate distance
					EM_object.estimatePhiWithEM(genoQual1, genoQual2);

					//write to file
					writeDistanceEstimates(out, curChr, windowStart, windowEnd, numSitesWithData, EM_object);


				} else writeDistanceEstimatesNoData(out, curChr, windowStart, windowEnd);
			} else writeDistanceEstimatesNoData(out, curChr, windowStart, windowEnd);

/*
			//print data
			for(int i=0; i<windowLen; ++i){
				std::cout << "Sample 1: " << genoQual1[i][0];
				for(int g=1; g<10; ++g){
					std::cout << "," << genoQual1[i][g];
				}
				std::cout << "\tSample 2: " << genoQual2[i][0];
				for(int g=1; g<10; ++g){
					std::cout << "," << genoQual2[i][g];
				}
				std::cout << std::endl;
			}
*/

			//move window
			windowStart = windowEnd;
			windowEnd = windowStart + windowLen;
			logfile().endIndent();
		}
		logfile().endNumbering();
	}

	logfile().endIndent();
};

//--------------------------------------------
// Writing estimates
//--------------------------------------------
void TDistanceEstimator::writeDistanceEstimates(gz::ogzstream & out, std::string & chr, uint32_t windowStart, uint32_t windowEnd, uint32_t numsitesWithData, TEMforDistanceEstimation & EM_object){
	out << chr << "\t" << windowStart + 1 << "\t" << windowEnd; //internal position is zero-based
	writeDistanceEstimates(out, numsitesWithData, EM_object);
};

void TDistanceEstimator::writeDistanceEstimates(gz::ogzstream & out, int numsitesWithData, TEMforDistanceEstimation & EM_object){
	using coretools::index;
	out << "\t" << numsitesWithData;
	//write pi
	for(Base b = Base::min; b < Base::max; ++b){
		out << "\t" << EM_object.pi[b];
	}
	//write phi
	for(DistancePhi p = DistancePhi::min; p < DistancePhi::max; ++p){
		out << "\t" << EM_object.phi[p];
	}
	//write distance
	out << "\t" << EM_object.distance;
	out << "\n";
};

void TDistanceEstimator::writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, uint32_t windowStart, uint32_t windowEnd){
	out << chr << "\t" << windowStart + 1 << "\t" << windowEnd << "\t"; //internal position is zero-based
	writeDistanceEstimatesNoData(out);
};

void TDistanceEstimator::writeDistanceEstimatesNoData(gz::ogzstream & out){
	out << "\t0";
	for(int i=0; i<14; ++i)
		out << "\t-";
	out << "\n";
};

}; //end namespace
