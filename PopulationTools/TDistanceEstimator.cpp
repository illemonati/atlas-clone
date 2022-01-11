/*
 * TDistanceCalculator.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#include "TDistanceEstimator.h"

namespace PopulationTools{

//------------------------------------------------
// DistancePhi
//------------------------------------------------
std::string DistancePhi::_toString[nn_nn + 1] = {"aa_aa", "aa_ab", "ab_aa", "aa_bb", "ab_ab", "ab_ac", "aa_bc", "ab_cc", "ab_cd", "nn_nn"};

std::ostream& operator<<(std::ostream& os, const DistancePhi & Phi){
	os << (std::string) Phi;
	return os;
};

//----------------------------------------------------
//TGenocombinationToBaseMap
//----------------------------------------------------
TGenocombinationToBaseMap::TGenocombinationToBaseMap(){
	for(Genotype g1 = Genotype::min(); g1 < Genotype::max(); ++g1){
		for(Genotype g2 = Genotype::min(); g2 < Genotype::max(); ++g2){
			for(Base b = Base::min(); b < Base::max(); ++b){
				genotypeCombinationHasBase[g1.get()][g2.get()][b.get()] = false;
			}

			genotypeCombinationHasBase[g1.get()][g2.get()][g1.firstAllele().get() ] = true;
			genotypeCombinationHasBase[g1.get()][g2.get()][g1.secondAllele().get()] = true;
			genotypeCombinationHasBase[g1.get()][g2.get()][g2.firstAllele().get() ] = true;
			genotypeCombinationHasBase[g1.get()][g2.get()][g2.secondAllele().get()] = true;
		}
	}
};

//----------------------------------------------------
//TDistanceClass
//----------------------------------------------------
TDistance::TDistance(){
	//squared difference between genotypes
	_distanceWeight[aa_aa] = 0.0;
	_distanceWeight[ab_aa] = 1.0;
	_distanceWeight[aa_ab] = 1.0;
	_distanceWeight[aa_bb] = 4.0;
	_distanceWeight[ab_ab] = 0.0;
	_distanceWeight[ab_ac] = 1.0;
	_distanceWeight[aa_bc] = 4.0;
	_distanceWeight[ab_cc] = 4.0;
	_distanceWeight[ab_cd] = 4.0;
};

double TDistance::calculateDistance(const TDistanceData & phi){
	return    phi[aa_aa] * _distanceWeight[aa_aa]
		    + phi[ab_aa] * _distanceWeight[ab_aa]
			+ phi[aa_ab] * _distanceWeight[aa_ab]
			+ phi[aa_bb] * _distanceWeight[aa_bb]
			+ phi[ab_ab] * _distanceWeight[ab_ab]
			+ phi[ab_ac] * _distanceWeight[ab_ac]
			+ phi[aa_bc] * _distanceWeight[aa_bc]
			+ phi[ab_cc] * _distanceWeight[ab_cc]
			+ phi[ab_cd] * _distanceWeight[ab_cd];
};

TDistanceProbMismatch::TDistanceProbMismatch(){
	//probability that a random allele from each individual is different
	_distanceWeight[aa_aa] = 0.0;
	_distanceWeight[ab_aa] = 0.5;
	_distanceWeight[aa_ab] = 0.5;
	_distanceWeight[aa_bb] = 1.0;
	_distanceWeight[ab_ab] = 0.5;
	_distanceWeight[ab_ac] = 0.75;
	_distanceWeight[aa_bc] = 1.0;
	_distanceWeight[ab_cc] = 1.0;
	_distanceWeight[ab_cd] = 1.0;
};

double TDistanceEuclidian::calculateDistance(const TDistanceData & phi){
	return sqrt(TDistance::calculateDistance(phi));
};

TDistanceUser::TDistanceUser(std::vector<double> vec){
	for(DistancePhi d = DistancePhi::min(); d < DistancePhi::max(); ++d)
		_distanceWeight[d] = vec[d.get()];
};

//----------------------------------------------------
//TDistanceEstimate
//----------------------------------------------------
TEMforDistanceEstimation::TEMforDistanceEstimation(coretools::TLog* Logfile, coretools::TParameters & params){
	logfile = Logfile;

	//prepare storage
	phi.set(0.0);
	LL = 0.0;
	old_LL = 0.0;

	probGeno = new double*[10];
	P_G = new double*[10];
	P_G_one_site = new double*[10];
	for(int i=0; i<10; ++i){
		probGeno[i] = new double[10];
		P_G[i] = new double[10];
		P_G_one_site[i] = new double[10];
	}

	//other variables
	distance = -1.0;

	//read EM parameters
	logfile->startIndent("Parameters of EM algorithm:");
	maxNumEMIterations = params.getParameterWithDefault<int>("iterations", 100);
	logfile->list("Will run up to ", maxNumEMIterations, " iterations.");
	epsilonForEM = params.getParameterWithDefault("maxEps", 0.000001);
	logfile->list("Will run EM until deltaLL < ", epsilonForEM, ".");
	logfile->endIndent();

	//set how to calculate distances
//	distanceWeight = new double[9];
	if(params.parameterExists("distWeights")){
		logfile->list("Using user-provided distance weights.");
		std::vector<double> vec;
		std::vector<std::string> tmp;
		params.fillParameterIntoContainer("distWeights", tmp, ',');
		coretools::str::repeatIndexes(tmp, vec);
		if(vec.size() != 9)
			throw "Wrong number of distance weights! Required are nine values for 00/00, 00/01, 01/00, 00/11, 01/01, 01/02, 00/12, 01/22, 01/23";

		distanceObject = new TDistanceUser(vec);

	} else {
		std::string distType = params.getParameterWithDefault<std::string>("distType", "squaredDiff");
		logfile->list("Using distance type '" + distType + "'.");
		if(distType == "probMismatch"){
			distanceObject = new TDistanceProbMismatch();
		} else if(distType == "squaredDiff"){
			distanceObject = new TDistance();
		} else if(distType == "euclidian"){
			distanceObject = new TDistanceEuclidian();
		} else
			throw "Unknown distance type '" + distType + "'! Use probMismatch.";
	}
	logfile->conclude("Using distance weights " + coretools::str::concatenateString(distanceObject->weights(), ", ") + ".");
};

void TEMforDistanceEstimation::guessPi(std::vector<HighPrecisionPhredIntProbability*> & genoQual1, std::vector<HighPrecisionPhredIntProbability*> & genoQual2){
	//check sizes are equal
	if(genoQual1.size() != genoQual2.size())
		throw "Provided genotype quality vectors are of different size in TEMforDistanceEstimation::guessPi!";

	//just estimate pi as average posterior probability
	pi.set(0.0);

	//now loop over sites
	auto it1 = genoQual1.begin();
	auto it2 = genoQual2.begin();
	for(; it1 != genoQual1.end(); ++it1, ++it2){
		double sum1 = 0.0;
		double sum2 = 0.0;
		for(Genotype g = Genotype::min(); g < Genotype::max(); ++g){
			sum1 += (Probability) (*it1)[g.get()];
			sum2 += (Probability) (*it2)[g.get()];
		}
		for(Genotype g = Genotype::min(); g < Genotype::max(); ++g){
			double tmp = (Probability) (*it1)[g.get()] / sum1;
			pi[g.firstAllele()] += tmp;
			pi[g.secondAllele()] += tmp;
			tmp = (Probability) (*it2)[g.get()] / sum2;
			pi[g.firstAllele()] += tmp;
			pi[g.secondAllele()] += tmp;
		}
	}

	//normalize
	pi.normalize();
}

void TEMforDistanceEstimation::guessPhi(std::vector<HighPrecisionPhredIntProbability*> & genoQual1, std::vector<HighPrecisionPhredIntProbability*> & genoQual2){
	//check sizes are equal
	if(genoQual1.size() != genoQual2.size())
		throw "Provided genotype quality vectors are of different size in TEMforDistanceEstimation::guessPhi!";

	//set to zero
	phi.set(0.0);

	//now loop over sites and add posterior probs
	auto it1 = genoQual1.begin();
	auto it2 = genoQual2.begin();
	for(; it1 != genoQual1.end(); ++it1, ++it2){
		double sum1 = 0.0; double sum2 = 0.0;
		for(Genotype g = Genotype::min(); g < Genotype::max(); ++g){
			sum1 += (Probability) (*it1)[g.get()];
			sum2 += (Probability) (*it2)[g.get()];
		}
		for(Genotype g1 = Genotype::min(); g1 < Genotype::max(); ++g1){
			double tmp = ((Probability) (*it1)[g1.get()] / sum1);
			for(Genotype g2 = Genotype::min(); g2 < Genotype::max(); ++g2){
				phi(g1, g2) += tmp * ((Probability) (*it2)[g2.get()] / sum2);
			}
		}
	}

	//normalize
	phi.normalize();
}

void TEMforDistanceEstimation::fill_K(TBaseData & thesePi){
	//normalizing constant for each phi class
	//case of one base
	K[aa_aa] = 1.0;

	//cases of two bases
	K[ab_ab] = thesePi[genometools::A] * thesePi[genometools::C]
         + thesePi[genometools::A] * thesePi[genometools::G]
         + thesePi[genometools::A] * thesePi[genometools::T]
		 + thesePi[genometools::C] * thesePi[genometools::G]
		 + thesePi[genometools::C] * thesePi[genometools::T]
		 + thesePi[genometools::G] * thesePi[genometools::T];
	K[aa_ab] = 2.0 * K[ab_ab]; //account for AC vs CA
	K[ab_aa] = K[aa_ab];
	K[aa_bb] = K[aa_ab];

	//cases of three bases
	K[aa_bc] = thesePi[genometools::C] * thesePi[genometools::G] * thesePi[genometools::T]
		 + thesePi[genometools::A] * thesePi[genometools::G] * thesePi[genometools::T]
		 + thesePi[genometools::A] * thesePi[genometools::C] * thesePi[genometools::T]
		 + thesePi[genometools::A] * thesePi[genometools::C] * thesePi[genometools::G];
	K[aa_bc] = 3.0 * K[aa_bc]; //account for ways to distribute

	K[ab_cc] = K[aa_bc];
	K[ab_ac] = 2.0 * K[aa_bc]; //twice as many cases than other cases with three bases!

	//case of four bases: each of the 6 cases is equally likely
	//Note: product of pi's cancels out when calculating P_g_given_phi_pi
	K[ab_cd] = 6.0;
};

void TEMforDistanceEstimation::fill_P_g_given_phi_pi(const TDistanceData & thesePhi, TBaseData & thesePi){
	//0: case aa/aa (K[0]=1)
	probGeno[genometools::AA][genometools::AA] = thesePhi[aa_aa] * thesePi[genometools::A];
	probGeno[genometools::CC][genometools::CC] = thesePhi[aa_aa] * thesePi[genometools::C];
	probGeno[genometools::GG][genometools::GG] = thesePhi[aa_aa] * thesePi[genometools::G];
	probGeno[genometools::TT][genometools::TT] = thesePhi[aa_aa] * thesePi[genometools::T];

	//1: cases aa/ab
	double tmp = thesePhi[aa_ab] / K[aa_aa];
	double tmp2 = tmp * thesePi[genometools::A];
	probGeno[genometools::AA][genometools::AC] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::AA][genometools::AG] = tmp2 * thesePi[genometools::G];
	probGeno[genometools::AA][genometools::AT] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::C];
	probGeno[genometools::CC][genometools::AC] = tmp2 * thesePi[genometools::A];
	probGeno[genometools::CC][genometools::CG] = tmp2 * thesePi[genometools::G];
	probGeno[genometools::CC][genometools::CT] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::G];
	probGeno[genometools::GG][genometools::AG] = tmp2 * thesePi[genometools::A];
	probGeno[genometools::GG][genometools::CG] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::GG][genometools::GT] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::T];
	probGeno[genometools::TT][genometools::AT] = tmp2 * thesePi[genometools::A];
	probGeno[genometools::TT][genometools::CT] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::TT][genometools::GT] = tmp2 * thesePi[genometools::G];

	//2: case ab/aa
	tmp = thesePhi[ab_aa] / K[ab_aa];
	tmp2 = tmp * thesePi[genometools::A];
	probGeno[genometools::AC][genometools::AA] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::AG][genometools::AA] = tmp2 * thesePi[genometools::G];
	probGeno[genometools::AT][genometools::AA] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::C];
	probGeno[genometools::AC][genometools::CC] = tmp2 * thesePi[genometools::A];
	probGeno[genometools::CG][genometools::CC] = tmp2 * thesePi[genometools::G];
	probGeno[genometools::CT][genometools::CC] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::G];
	probGeno[genometools::AG][genometools::GG] = tmp2 * thesePi[genometools::A];
	probGeno[genometools::CG][genometools::GG] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::GT][genometools::GG] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::T];
	probGeno[genometools::AT][genometools::TT] = tmp2 * thesePi[genometools::A];
	probGeno[genometools::CT][genometools::TT] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::GT][genometools::TT] = tmp2 * thesePi[genometools::G];

	//3: case aa/bb
	tmp = thesePhi[aa_bb] / K[aa_bb];
	tmp2 = tmp * thesePi[genometools::A];
	probGeno[genometools::AA][genometools::CC] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::AA][genometools::GG] = tmp2 * thesePi[genometools::G];
	probGeno[genometools::AA][genometools::TT] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::C];
	probGeno[genometools::CC][genometools::AA] = tmp2 * thesePi[genometools::A];
	probGeno[genometools::CC][genometools::GG] = tmp2 * thesePi[genometools::G];
	probGeno[genometools::CC][genometools::TT] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::G];
	probGeno[genometools::GG][genometools::AA] = tmp2 * thesePi[genometools::A];
	probGeno[genometools::GG][genometools::CC] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::GG][genometools::TT] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::T];
	probGeno[genometools::TT][genometools::AA] = tmp2 * thesePi[genometools::A];
	probGeno[genometools::TT][genometools::CC] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::TT][genometools::GG] = tmp2 * thesePi[genometools::G];

	//4: case ab/ab
	tmp = thesePhi[ab_ab] / K[ab_ab];
	tmp2 = tmp * thesePi[genometools::A];
	probGeno[genometools::AC][genometools::AC] = tmp2 * thesePi[genometools::C];
	probGeno[genometools::AG][genometools::AG] = tmp2 * thesePi[genometools::G];
	probGeno[genometools::AT][genometools::AT] = tmp2 * thesePi[genometools::T];
	tmp2 = tmp * thesePi[genometools::C];
	probGeno[genometools::CG][genometools::CG] = tmp2 * thesePi[genometools::G];
	probGeno[genometools::CT][genometools::CT] = tmp2 * thesePi[genometools::T];
	probGeno[genometools::GT][genometools::GT] = tmp * thesePi[genometools::G] * thesePi[genometools::T];

	//5: case ab/ac
	tmp = thesePhi[ab_ac] / K[ab_ac];
	tmp2 = tmp * thesePi[genometools::A] * thesePi[genometools::C] * thesePi[genometools::G];
	probGeno[genometools::AC][genometools::AG] = tmp2;
	probGeno[genometools::AC][genometools::CG] = tmp2;
	probGeno[genometools::AG][genometools::AC] = tmp2;
	probGeno[genometools::AG][genometools::CG] = tmp2;
	probGeno[genometools::CG][genometools::AC] = tmp2;
	probGeno[genometools::CG][genometools::AG] = tmp2;
	tmp2 = tmp * thesePi[genometools::A] * thesePi[genometools::C] * thesePi[genometools::T];
	probGeno[genometools::AC][genometools::AT] = tmp2;
	probGeno[genometools::AC][genometools::CT] = tmp2;
	probGeno[genometools::AT][genometools::AC] = tmp2;
	probGeno[genometools::AT][genometools::CT] = tmp2;
	probGeno[genometools::CT][genometools::AC] = tmp2;
	probGeno[genometools::CT][genometools::AT] = tmp2;
	tmp2 = tmp * thesePi[genometools::A] * thesePi[genometools::G] * thesePi[genometools::T];
	probGeno[genometools::AG][genometools::GT] = tmp2;
	probGeno[genometools::AG][genometools::AT] = tmp2;
	probGeno[genometools::AT][genometools::AG] = tmp2;
	probGeno[genometools::AT][genometools::GT] = tmp2;
	probGeno[genometools::GT][genometools::AG] = tmp2;
	probGeno[genometools::GT][genometools::AT] = tmp2;
	tmp2 = tmp * thesePi[genometools::C] * thesePi[genometools::G] * thesePi[genometools::T];
	probGeno[genometools::CG][genometools::CT] = tmp2;
	probGeno[genometools::CG][genometools::GT] = tmp2;
	probGeno[genometools::CT][genometools::CG] = tmp2;
	probGeno[genometools::CT][genometools::GT] = tmp2;
	probGeno[genometools::GT][genometools::CG] = tmp2;
	probGeno[genometools::GT][genometools::CT] = tmp2;

	//6: case aa/bc
	tmp = thesePhi[aa_bc] / K[aa_bc];
	tmp2 = tmp * thesePi[genometools::A] * thesePi[genometools::C] * thesePi[genometools::G];
	probGeno[genometools::AA][genometools::CG] = tmp2;
	probGeno[genometools::CC][genometools::AG] = tmp2;
	probGeno[genometools::GG][genometools::AC] = tmp2;
	tmp2 = tmp * thesePi[genometools::A] * thesePi[genometools::C] * thesePi[genometools::T];
	probGeno[genometools::AA][genometools::CT] = tmp2;
	probGeno[genometools::CC][genometools::AT] = tmp2;
	probGeno[genometools::TT][genometools::AC] = tmp2;
	tmp2 = tmp * thesePi[genometools::A] * thesePi[genometools::G] * thesePi[genometools::T];
	probGeno[genometools::AA][genometools::GT] = tmp2;
	probGeno[genometools::GG][genometools::AT] = tmp2;
	probGeno[genometools::TT][genometools::AG] = tmp2;
	tmp2 = tmp * thesePi[genometools::C] * thesePi[genometools::G] * thesePi[genometools::T];
	probGeno[genometools::CC][genometools::GT] = tmp2;
	probGeno[genometools::GG][genometools::CT] = tmp2;
	probGeno[genometools::TT][genometools::CG] = tmp2;

	//7: case ab/cc
	tmp = thesePhi[ab_cc] / K[ab_cc];
	tmp2 = tmp * thesePi[genometools::A] * thesePi[genometools::C] * thesePi[genometools::G];
	probGeno[genometools::AC][genometools::GG] = tmp2;
	probGeno[genometools::AG][genometools::CC] = tmp2;
	probGeno[genometools::CG][genometools::AA] = tmp2;
	tmp2 = tmp * thesePi[genometools::A] * thesePi[genometools::C] * thesePi[genometools::T];
	probGeno[genometools::AC][genometools::TT] = tmp2;
	probGeno[genometools::AT][genometools::CC] = tmp2;
	probGeno[genometools::CT][genometools::AA] = tmp2;
	tmp2 = tmp * thesePi[genometools::A] * thesePi[genometools::G] * thesePi[genometools::T];
	probGeno[genometools::AG][genometools::TT] = tmp2;
	probGeno[genometools::AT][genometools::GG] = tmp2;
	probGeno[genometools::GT][genometools::AA] = tmp2;
	tmp2 = tmp * thesePi[genometools::C] * thesePi[genometools::G] * thesePi[genometools::T];
	probGeno[genometools::CG][genometools::TT] = tmp2;
	probGeno[genometools::CT][genometools::GG] = tmp2;
	probGeno[genometools::GT][genometools::CC] = tmp2;

	//8: case ab/cd
	tmp = thesePhi[ab_cd] / K[ab_cd];
	probGeno[genometools::AC][genometools::GT] = tmp;
	probGeno[genometools::AG][genometools::CT] = tmp;
	probGeno[genometools::AT][genometools::CG] = tmp;
	probGeno[genometools::CG][genometools::AT] = tmp;
	probGeno[genometools::CT][genometools::AG] = tmp;
	probGeno[genometools::GT][genometools::AC] = tmp;
};

bool TEMforDistanceEstimation::estimatePhiWithEM(std::vector<HighPrecisionPhredIntProbability*> & genoQual1, std::vector<HighPrecisionPhredIntProbability*> & genoQual2){
	//prepare estimates
	logfile->listFlush("Estimating initial base frequencies pi ...");
	guessPi(genoQual1, genoQual2);
	logfile->done();
	logfile->conclude("Initial pi are A=", pi[genometools::A], ", C=", pi[genometools::C], ", G=", pi[genometools::G], " and T=", pi[genometools::T], ".");
	logfile->listFlush("Estimating initial genotype classes phi ...");
	guessPhi(genoQual1, genoQual2);
	logfile->done();
	logfile->conclude("Initial phi are " + coretools::str::concatenateString(phi, ", ") + ".");

	//variables
	double old_LL, LL = 0.0;
	double LL_diff;

	//now run EM
	logfile->startIndent("Estimating phi using an EM algorithm:");
	for(int iter=0; iter<maxNumEMIterations; ++iter){
		logfile->listFlush("Running EM iteration ", iter+1, " ...");
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
		std::vector<HighPrecisionPhredIntProbability*>::iterator it1 = genoQual1.begin();
		std::vector<HighPrecisionPhredIntProbability*>::iterator it2 = genoQual2.begin();
		for(; it1 != genoQual1.end(); ++it1, ++it2){
			//calculate P_G per site
			double sum = 0.0;
			for(int g1 = 0; g1<10; ++g1){
				for(int g2 = 0; g2<10; ++g2){
					P_G_one_site[g1][g2] = (coretools::Probability) (*it1)[g1] * (coretools::Probability) (*it2)[g2] * probGeno[g1][g2];

					//std::cout << g1 << "/" << g2 << ": " <<  phredToLik[genoQual1[s][g1]] << " * " << phredToLik[genoQual2[s][g2]] << " * " << probGeno[g1][g2] << " = " << P_G_one_site[g1][g2] << std::endl;

					sum += P_G_one_site[g1][g2];
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
		phi.set(0.0);

		double sum = 0.0;
		for(Genotype g1 = Genotype::min(); g1 < Genotype::max(); ++g1){
			for(Genotype g2 = Genotype::min(); g2 < Genotype::max(); ++g2){
				phi(g1,g2) += P_G[g1.get()][g2.get()];
				sum += P_G[g1.get()][g2.get()];
			}
		}
		phi.normalize(sum);

		//update pi
		pi.set(0.0);
		for(Genotype g1 = Genotype::min(); g1 < Genotype::max(); ++g1){
			for(Genotype g2 = Genotype::min(); g2 < Genotype::max(); ++g2){
				for(Base b = Base::min(); b < Base::max(); ++b){
					if(genoToBaseMap(g1, g2, b)){
						pi[b] += P_G[g1.get()][g2.get()];
					}
				}
			}
		}
		pi.normalize();

		//check if EM converged
		logfile->done();
		//throw "done!";
		if(iter > 0 ){
			LL_diff = LL - old_LL;
			logfile->conclude("LL = ", LL, " (deltaLL = ", LL_diff, ").");
			if(LL_diff < epsilonForEM){
				logfile->conclude("EM converged, delatLL = ", LL_diff, " < ", epsilonForEM);
				distance = distanceObject->calculateDistance(phi);
				logfile->conclude("Resulting distance is ", distance);
				logfile->endIndent();
				return true;
			}
		} else
			logfile->conclude("LL = ", LL, ".");
	}
	logfile->warning("EM reached maximum number of iterations (", maxNumEMIterations, ") without converging!");
	distance = distanceObject->calculateDistance(phi);
	logfile->conclude("Resulting distance is ", distance);
	logfile->endIndent();
	return false;
};

//----------------------------------------------------
//TDistanceEstimator
//----------------------------------------------------
TDistanceEstimator::TDistanceEstimator(coretools::TLog* Logfile, coretools::TParameters & params){
	logfile = Logfile;
	maxNumEMIterations = 0;
	epsilonForEM = 0.0;
	numGLFs = 0;
	glfs = NULL;
	readersOpened = false;

	//outputname
	outputName = params.getParameterWithDefault<std::string>("out", "ATLAS");
	logfile->list("Writing output files with prefix '" + outputName + "'. (parameter 'out')");
}

void TDistanceEstimator::printGLF(coretools::TParameters & params){
	//test first to parse GLF files
	std::string glf = params.getParameter<std::string>("glf");
	GLF::TGlfReader reader(glf);

	//print file
	reader.printToEnd();
}

void TDistanceEstimator::openGLF(coretools::TParameters & params){
	params.fillParameterIntoContainer("glf", GLFNames, ',');
	numGLFs = GLFNames.size();
	if(numGLFs < 2)
		throw "At least two GLF files have to be provided to estimate distances!";

	//open files
	glfs = new GLF::TGlfReader[numGLFs];
	readersOpened = true;
	logfile->startIndent("Opening GLF files:");
	int g = 0;
	for(std::vector<std::string>::iterator it=GLFNames.begin(); it != GLFNames.end(); ++it, ++g){
		logfile->listFlush("Opening GLF '" + *it + "' ...");
		glfs[g].open(*it);
		logfile->done();
	}
	logfile->endIndent();
}

void TDistanceEstimator::closeGLF(){
	if(readersOpened){
		//close all glf handlers
		for(int g=0; g<numGLFs; ++g)
			glfs[g].close();

		delete[] glfs;
		GLFNames.clear();
		numGLFs = 0;
		readersOpened = false;
	}
}

//------------------------------------------------------------------
void TDistanceEstimator::estimateDistances(coretools::TParameters & params){
	//open all GLF files specified
	openGLF(params);

	//open EM object
	TEMforDistanceEstimation EM_object(logfile, params);

	//in windows or whole genome?
	long windowLen = params.getParameterWithDefault("window", -1L);
	if(windowLen < 0)
		estimateDistanceGenomeWide(EM_object);
	else
		estimateDistanceInWindows(EM_object, windowLen);

	//close files
	closeGLF();
}

//--------------------------------------------
// Estimation Genome Wide
//--------------------------------------------
void TDistanceEstimator::estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object){
	logfile->list("Will estimate genetic distances genome wide.");

	//open output file
	std::string filename = outputName + "_distanceEstimates.txt.gz";
	gz::ogzstream out(filename.c_str());
	if(!out)
		throw "Failed to open output file '" + filename + "'!";

	//write header to output file
	out << "individual1\tindividual2\tnumSitesWithData\tfreqA\tfreqC\tfreqG\tfreqT\tfreq00_00\tfreq00_01\tfreq01_00\tfreq00_11\tfreq01_01\tfreq01_02\tfreq00_12\tfreq01_22\tfreq01_23\tgeneticDist\n";

	//prepare storage for distance matrix
	double** distMatrix = new double*[numGLFs];
	for(int g=0; g<numGLFs; ++g){
		distMatrix[g] = new double[numGLFs];
		distMatrix[g][g] = 0.0;
	}

	//loop over all pairs
	for(int g1=0; g1<(numGLFs-1); ++g1){
		for(int g2 = g1+1; g2 < numGLFs; ++g2){
			logfile->startIndent("Estimating distance between individuals ", g1+1, " (" + GLFNames[g1], ") and ", g2+1, " (", GLFNames[g2], "):");

			//write names to file
			out << GLFNames[g1] << "\t" << GLFNames[g2];

			//run estimation
			estimateDistanceGenomeWide(EM_object, glfs[g1], glfs[g2], out);

			//write to matrix
			distMatrix[g1][g2] = EM_object.distance;
			distMatrix[g2][g1] = EM_object.distance;
			logfile->endIndent();
		}
	}

	out.close();

	//open matrix file
	filename = outputName + "_distanceMatrix.txt";
	std::ofstream distMatrixFile(filename.c_str());
	if(!distMatrixFile)
		throw "Failed to open output file '" + filename + "'!";

	//write header to matrix file
	distMatrixFile << "/";
	for(int g=0; g<numGLFs; ++g)
		distMatrixFile << "\t" << GLFNames[g];
	distMatrixFile << "\n";

	//write rows
	for(int g1 = 0; g1 < numGLFs; ++g1){
		distMatrixFile << GLFNames[g1];
		for(int g2 = 0; g2 < numGLFs; ++g2)
			distMatrixFile << "\t" << distMatrix[g1][g2];
		distMatrixFile << "\n";
	}

	//close file
	distMatrixFile.close();
};

bool TDistanceEstimator::moveToNextCommonChr(GLF::TGlfReader & g1, GLF::TGlfReader & g2){
	while(g1.refId() != g2.refId() || g2.chrIsHaploid()){
		//advance the one laging behind
		if(g1.refId() < g2.refId()){
			if(!g1.jumpToNextChr()) return false;
		} else {
			if(!g2.jumpToNextChr()) return false;
		}
	}

	//check names
	if(g1.chr() != g2.chr()){
		throw "Chromosome names differ in files " + g1.name() + "' and '" + g2.name() + "': '" + g1.chr() + "' != '" + g2.chr() + "'!";
	}
	if(g1.chrLength() != g2.chrLength()){
		throw coretools::str::toString("Chromosome lengths differ in files " + g1.name() + "' and '" + g2.name() + "': '",  g1.chrLength(), "' != '", g2.chrLength(), "'!");
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

void TDistanceEstimator::readCommonSites(std::vector<HighPrecisionPhredIntProbability*> & genoQual1, std::vector<HighPrecisionPhredIntProbability*> & genoQual2, GLF::TGlfReader & g1, GLF::TGlfReader & g2){
	//parse GLFs. Only keep sites where both individuals have data!

	//rewind GLFs
	g1.rewind();
	g2.rewind();

	//if not both are good at least one file reach end. So we are done!
	while(advance(g1, g2)){
		if(g2.position() == g1.position()){
			//add data
			genoQual1.push_back(new HighPrecisionPhredIntProbability[10]);
			genoQual2.push_back(new HighPrecisionPhredIntProbability[10]);
			g1.fillGenotypeLikelihoodsGLF(*genoQual1.rbegin());
			g2.fillGenotypeLikelihoodsGLF(*genoQual2.rbegin());
		}
	}
};

void TDistanceEstimator::estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object, GLF::TGlfReader & g1, GLF::TGlfReader & g2, gz::ogzstream & out){
	//initialize storage for two windows
	logfile->listFlush("Reading common sites ...");
	std::vector<HighPrecisionPhredIntProbability*> genoQual1, genoQual2;
	readCommonSites(genoQual1, genoQual2, g1, g2);
	logfile->done();
	logfile->conclude("Read data for ", genoQual1.size(), " sites.");

	//now estimate
	if(genoQual1.size() > 0){
		logfile->startIndent("Estimating genetic distance:");
		EM_object.estimatePhiWithEM(genoQual1, genoQual2);
		writeDistanceEstimates(out, genoQual1.size(), EM_object);
		logfile->endIndent();
	} else {
		logfile->conclude("Not enough data to estimate distance.");
		writeDistanceEstimatesNoData(out);
	}

	//clean up memory
	logfile->listFlush("Cleaning up memory ...");
	for(auto& it : genoQual1){
		delete[] it;
	}
	for(auto& it : genoQual2){
		delete[] it;
	}
	logfile->done();
};

//--------------------------------------------
// Estimation in windows
//--------------------------------------------
void TDistanceEstimator::estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, uint32_t windowLen){
	logfile->list("Will estimate genetic distance in windows of length ", windowLen, ".");
	if(windowLen < 100)
		throw "Window size must be at least 100bp!";

	//loop over all pairs
	for(int g1=0; g1<(numGLFs-1); ++g1){
		for(int g2 = g1+1; g2 < numGLFs; ++g2){
			logfile->startIndent("Estimating distance between individuals ", g1+1, " (", GLFNames[g1], ") and ", g2+1, " (", GLFNames[g2], "):");

			//output file
			std::string filename = outputName + "_" + GLFNames[g1] + "_" + GLFNames[g2] + "_distanceEstimates.txt.gz";
			logfile->list("Will write estimates to file '" + filename + "'.");

			//rewind GLFs
			glfs[g1].rewind();
			glfs[g2].rewind();

			//now run estimation
			estimateDistanceInWindows(EM_object, filename, glfs[g1], glfs[g2], windowLen);

			logfile->endIndent();
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
	std::vector<HighPrecisionPhredIntProbability*> genoQual1, genoQual2;
	genoQual1.resize(windowLen);
	genoQual2.resize(windowLen);
	for(size_t i=0; i<windowLen; ++i){
		genoQual1[i] = new HighPrecisionPhredIntProbability[10];
		genoQual2[i] = new HighPrecisionPhredIntProbability[10];
	}

	//open output file
	gz::ogzstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

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
	logfile->startIndent("Will estimate distance in windows of size ", windowLen, ":");
	while(!g1.eof() && !g2.eof()){
		//move to new chromosome
		curRefId = g1.refId();
		curChr = g1.chr();
		curChrLen = g1.chrLength();
		windowStart = 0;
		windowEnd = windowLen;

		logfile->startNumbering("Chromosome " + curChr + ":");

		//parse all windows of chromosome
		while(windowStart < curChrLen && !g1.eof() && !g2.eof()){
			logfile->number("Window [", windowStart, ", ", windowEnd, ")");
			logfile->addIndent();

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
			logfile->endIndent();
		}
		logfile->endNumbering();
	}

	//clean up memory
	for(size_t i=0; i<windowLen; ++i){
		delete[] genoQual1[i];
		delete[] genoQual2[i];
	}

	logfile->endIndent();
};

//--------------------------------------------
// Writing estimates
//--------------------------------------------
void TDistanceEstimator::writeDistanceEstimates(gz::ogzstream & out, std::string & chr, uint32_t windowStart, uint32_t windowEnd, uint32_t numsitesWithData, TEMforDistanceEstimation & EM_object){
	out << chr << "\t" << windowStart + 1 << "\t" << windowEnd; //internal position is zero-based
	writeDistanceEstimates(out, numsitesWithData, EM_object);
};

void TDistanceEstimator::writeDistanceEstimates(gz::ogzstream & out, int numsitesWithData, TEMforDistanceEstimation & EM_object){
	out << "\t" << numsitesWithData;
	//write pi
	for(Base b = Base::min(); b < Base::max(); ++b){
		out << "\t" << EM_object.pi[b];
	}
	//write phi
	for(DistancePhi p = DistancePhi::min(); p < DistancePhi::max(); ++p){
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
