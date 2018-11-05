/*
 * TDistanceCalculator.cpp
 *
 *  Created on: Jul 25, 2017
 *      Author: phaentu
 */

#include "TDistanceEstimator.h"

//-------------------------------------------------
//TPhiToGenoMap
//-------------------------------------------------
TGenoToPhiMap::TGenoToPhiMap(){
	//--------------------------
	//mapping genotypes to phi
	//--------------------------
	genoToPhiMap = new int*[10];
	for(int i=0; i<10; ++i)
		genoToPhiMap[i] = new int[10];

	//case aa/aa
	genoToPhiMap[AA][AA] = 0;
	genoToPhiMap[CC][CC] = 0;
	genoToPhiMap[GG][GG] = 0;
	genoToPhiMap[TT][TT] = 0;

	Genotype G1, G2;
	int b1, b2;
	for(b1=0; b1<4; ++b1){
		G1 = genoMap.getGenotype(b1, b1);
		genoToPhiMap[G1][G1] = 0;
	}

	//case aa/ab and ab/aa AND aa/bb
	for(b1=0; b1<4; ++b1){
		for(b2=0; b2<4; ++b2){
			if(b2 != b1){
				//aa/ab and ab/aa
				G1 = genoMap.getGenotype(b1, b1);
				G2 = genoMap.getGenotype(b1, b2);
				genoToPhiMap[G1][G2] = 1;
				genoToPhiMap[G2][G1] = 2;

				//aa/bb
				G1 = genoMap.getGenotype(b1, b1);
				G2 = genoMap.getGenotype(b2, b2);
				genoToPhiMap[G1][G2] = 3;
			}
		}
	}

	//case ab/ab
	for(b1=0; b1<3; ++b1){
		for(b2=b1+1; b2<4; ++b2){
			G1 = genoMap.getGenotype(b1, b2);
			genoToPhiMap[G1][G1] = 4;
		}
	}

	//case ab/ac
	int b3;
	for(b1=0; b1<4; ++b1){
		for(b2=0; b2<4; ++b2){
			if(b2 != b1){
				for(b3=0; b3<4; ++b3){
					if(b3 != b1 && b3 != b2){
						G1 = genoMap.getGenotype(b1, b2);
						G2 = genoMap.getGenotype(b1, b3);
						genoToPhiMap[G1][G2] = 5;
					}
				}
			}
		}
	}

	//case aa/bc and ab/cc
	for(b1=0; b1<4; ++b1){
		for(b2=0; b2<3; ++b2){
			for(b3=b2+1; b3<4; ++b3){
				if(b2 != b1 && b3 != b1){
					G1 = genoMap.getGenotype(b1, b1);
					G2 = genoMap.getGenotype(b2, b3);
					genoToPhiMap[G1][G2] = 6;
					genoToPhiMap[G2][G1] = 7;
				}
			}
		}
	}

	//case ab/cd
	genoToPhiMap[AC][GT] = 8;
	genoToPhiMap[AG][CT] = 8;
	genoToPhiMap[AT][CG] = 8;
	genoToPhiMap[CG][AT] = 8;
	genoToPhiMap[CT][AG] = 8;
	genoToPhiMap[GT][AC] = 8;


	//test (seemed correct on 24/8/2017)
	/*
	for(int g1 = 0; g1<10; ++g1){
		for(int g2 = 0; g2<10; ++g2){
			std::cout << genoMap.getGenotypeString(g1) << "/" << genoMap.getGenotypeString(g2) << " -> " << genoToPhiMap[g1][g2] << std::endl;
		}
	}
	*/

};

//----------------------------------------------------
//TGenocombinationToBaseMap
//----------------------------------------------------
TGenocombinationToBaseMap::TGenocombinationToBaseMap(){
	genotypeCombinationHasBase = new bool**[10];
	int g1, g2, b;
	for(g1 = 0; g1<10; ++g1){
		genotypeCombinationHasBase[g1] = new bool*[10];
		for(g2 = 0; g2<10; ++g2){
			genotypeCombinationHasBase[g1][g2] = new bool[4];
			for(b=0; b<4; ++b)
				genotypeCombinationHasBase[g1][g2][b] = false;

			genotypeCombinationHasBase[g1][g2][genoMap.genotypeToBase[g1][0]] = true;
			genotypeCombinationHasBase[g1][g2][genoMap.genotypeToBase[g1][1]] = true;
			genotypeCombinationHasBase[g1][g2][genoMap.genotypeToBase[g2][0]] = true;
			genotypeCombinationHasBase[g1][g2][genoMap.genotypeToBase[g2][1]] = true;
		}
	}
};
//----------------------------------------------------
//TDistanceClass
//----------------------------------------------------
TDistance::TDistance(){
	distanceWeight = new double[9];

	//squared difference between genotypes
	distanceWeight[0] = 0.0; //case aa/aa
	distanceWeight[1] = 1.0; //case ab/aa
	distanceWeight[2] = 1.0; //case aa/ab
	distanceWeight[3] = 4.0; //case aa/bb
	distanceWeight[4] = 0.0; //case ab/ab
	distanceWeight[5] = 1.0; //case ab/ac
	distanceWeight[6] = 4.0; //case aa/bc
	distanceWeight[7] = 4.0; //case ab/cc
	distanceWeight[8] = 4.0; //case ab/cd
}

double TDistance::calculateDistance(double* phi){
	double distance = 0.0;
	for(int i=0; i<9; ++i)
		distance += phi[i] * distanceWeight[i];
	return distance;
}

TDistanceProbMismatch::TDistanceProbMismatch(){
	//probability that a random allele from each individual is different
	distanceWeight[0] = 0.0; //case aa/aa
	distanceWeight[1] = 0.5; //case ab/aa
	distanceWeight[2] = 0.5; //case aa/ab
	distanceWeight[3] = 1.0; //case aa/bb
	distanceWeight[4] = 0.5; //case ab/ab
	distanceWeight[5] = 0.75; //case ab/ac
	distanceWeight[6] = 1.0; //case aa/bc
	distanceWeight[7] = 1.0; //case ab/cc
	distanceWeight[8] = 1.0; //case ab/cd
}

double TDistanceEuclidian::calculateDistance(double* phi){
	return sqrt(TDistance::calculateDistance(phi));
}

TDistanceUser::TDistanceUser(std::vector<double> vec){
	for(int i=0; i<9; ++i)
		distanceWeight[i] = vec[i];
}

//----------------------------------------------------
//TDistanceEstimate
//----------------------------------------------------
TEMforDistanceEstimation::TEMforDistanceEstimation(TLog* Logfile, TParameters & params){
	logfile = Logfile;

	//prepare storage
	phi = new double[9];
	for(int i=0; i<9; ++i){
		phi[i] = 0.0;
	}
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
	K = new double[9];

	//other variables
	distance = -1.0;

	//read EM parameters
	logfile->startIndent("Parameters of EM algorithm:");
	maxNumEMIterations = params.getParameterIntWithDefault("iterations", 100);
	logfile->list("Will run up to " + toString(maxNumEMIterations) + " iterations.");
	epsilonForEM = params.getParameterDoubleWithDefault("maxEps", 0.000001);
	logfile->list("Will run EM until deltaLL < " + toString(epsilonForEM) + ".");
	logfile->endIndent();

	//set how to calculate distances
//	distanceWeight = new double[9];
	if(params.parameterExists("distWeights")){
		logfile->list("Using user-provided distance weights.");
		std::vector<double> vec;
		std::vector<std::string> tmp;
		params.fillParameterIntoVector("distWeights", tmp, ',');
		repeatIndexes(tmp, vec);
		if(vec.size() != 9)
			throw "Wrong number of distance weights! Required are nine values for 00/00, 00/01, 01/00, 00/11, 01/01, 01/02, 00/12, 01/22, 01/23";
		distanceObject = new TDistanceUser(vec);
	} else {
		std::string distType = params.getParameterStringWithDefault("distType", "squaredDiff");
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
	logfile->conclude("Using distance weights " + concatenateString(distanceObject->distanceWeight, 9, ", ") + ".");
};

void TEMforDistanceEstimation::guessPi(std::vector<uint8_t*> & genoQual1, std::vector<uint8_t*> & genoQual2){
	//check sizes are equal
	if(genoQual1.size() != genoQual2.size())
		throw "Provided genotype quality vectors are of different size in TEMforDistanceEstimation::guessPi!";

	//just estimate pi as average posterior probability
	pi.clear();
	double sum1, sum2;

	//now loop over sites
	std::vector<uint8_t*>::iterator it1 = genoQual1.begin();
	std::vector<uint8_t*>::iterator it2 = genoQual2.begin();
	for(; it1 != genoQual1.end(); ++it1, ++it2){
		sum1 = 0.0; sum2 = 0.0;
		for(int i=0; i<10; ++i){
			sum1 += phredToLik[(*it1)[i]];
			sum2 += phredToLik[(*it2)[i]];
		}
		for(int i=0; i<10; ++i){
			double tmp = phredToLik[(*it1)[i]] / sum1;
			pi.addNoRef(genoMap.genotypeToBase[i][0], tmp);
			pi.addNoRef(genoMap.genotypeToBase[i][1], tmp);
			tmp = phredToLik[(*it2)[i]] / sum2;
			pi.addNoRef(genoMap.genotypeToBase[i][0], tmp);
			pi.addNoRef(genoMap.genotypeToBase[i][1], tmp);
		}
	}

	//normalize
	pi.normalize();
}

void TEMforDistanceEstimation::guessPhi(std::vector<uint8_t*> & genoQual1, std::vector<uint8_t*> & genoQual2){
	//check sizes are equal
	if(genoQual1.size() != genoQual2.size())
		throw "Provided genotype quality vectors are of different size in TEMforDistanceEstimation::guessPhi!";

	//set to zero
	for(int i=0; i<9; ++i)
		phi[i] = 0.0;

	//now loop over sites and add posterior probs
	double sum1, sum2;
	std::vector<uint8_t*>::iterator it1 = genoQual1.begin();
	std::vector<uint8_t*>::iterator it2 = genoQual2.begin();
	for(; it1 != genoQual1.end(); ++it1, ++it2){
		sum1 = 0.0; sum2 = 0.0;
		for(int i=0; i<10; ++i){
			sum1 += phredToLik[(*it1)[i]];
			sum2 += phredToLik[(*it2)[i]];
		}
		for(int g1 = 0; g1<10; ++g1){
			double tmp = (phredToLik[(*it1)[g1]] / sum1);
			for(int g2 = 0; g2<10; ++g2){
				phi[genoToPhiMap(g1,g2)] += tmp * (phredToLik[(*it2)[g2]] / sum2);
			}
		}
	}

	//normalize
	sum1 = 0.0;
	for(int i=0; i<9; ++i)
		sum1 += phi[i];
	for(int i=0; i<9; ++i)
		phi[i] /= sum1;
}

void TEMforDistanceEstimation::fill_K(TBaseFrequencies  & thesePi){
	//normalizing constant for each phi class
	//case of one base
	K[0] = 1.0;

	//cases of two bases
	K[4] = 0.0;
	for(int i=0; i<3; ++i){
		for(int j=i; j<4; ++j)
			K[4] += thesePi[i] * thesePi[j];
	}
	K[1] = 2.0 * K[4]; //account for AC vs CA
	K[2] = K[1];
	K[3] = K[1];

	//cases of three bases
	K[6] = thesePi[1] * thesePi[2] * thesePi[3]
		 + thesePi[0] * thesePi[2] * thesePi[3]
		 + thesePi[0] * thesePi[1] * thesePi[3]
		 + thesePi[0] * thesePi[1] * thesePi[2];
	K[6] = 3.0 * K[6]; //account for ways to distribute

	K[7] = K[6];
	K[5] = 2.0 * K[6]; //twice as many cases than other cases with three bases!

	//case of four bases: each of the 6 cases is equally likely
	K[8] = 6.0;
}

void TEMforDistanceEstimation::fill_P_g_given_phi_pi(double* thesePhi, TBaseFrequencies & thesePi){
	//case aa/aa
	probGeno[AA][AA] = thesePhi[0] * thesePi[A];
	probGeno[CC][CC] = thesePhi[0] * thesePi[C];
	probGeno[GG][GG] = thesePhi[0] * thesePi[G];
	probGeno[TT][TT] = thesePhi[0] * thesePi[T];

	//cases aa/ab
	double tmp = thesePhi[1] / K[1];
	double tmp2 = tmp * thesePi[A];
	probGeno[AA][AC] = tmp2 * thesePi[C];
	probGeno[AA][AG] = tmp2 * thesePi[G];
	probGeno[AA][AT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[C];
	probGeno[CC][AC] = tmp2 * thesePi[A];
	probGeno[CC][CG] = tmp2 * thesePi[G];
	probGeno[CC][CT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[G];
	probGeno[GG][AG] = tmp2 * thesePi[A];
	probGeno[GG][CG] = tmp2 * thesePi[C];
	probGeno[GG][GT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[T];
	probGeno[TT][AT] = tmp2 * thesePi[A];
	probGeno[TT][CT] = tmp2 * thesePi[C];
	probGeno[TT][GT] = tmp2 * thesePi[G];

	//case ab/aa
	tmp = thesePhi[2] / K[2];
	tmp2 = tmp * thesePi[A];
	probGeno[AC][AA] = tmp2 * thesePi[C];
	probGeno[AG][AA] = tmp2 * thesePi[G];
	probGeno[AT][AA] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[C];
	probGeno[AC][CC] = tmp2 * thesePi[A];
	probGeno[CG][CC] = tmp2 * thesePi[G];
	probGeno[CT][CC] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[G];
	probGeno[AG][GG] = tmp2 * thesePi[A];
	probGeno[CG][GG] = tmp2 * thesePi[C];
	probGeno[GT][GG] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[T];
	probGeno[AT][TT] = tmp2 * thesePi[A];
	probGeno[CT][TT] = tmp2 * thesePi[C];
	probGeno[GT][TT] = tmp2 * thesePi[G];

	//case aa/bb
	tmp = thesePhi[3] / K[3];
	tmp2 = tmp * thesePi[A];
	probGeno[AA][CC] = tmp2 * thesePi[C];
	probGeno[AA][GG] = tmp2 * thesePi[G];
	probGeno[AA][TT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[C];
	probGeno[CC][AA] = tmp2 * thesePi[A];
	probGeno[CC][GG] = tmp2 * thesePi[G];
	probGeno[CC][TT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[G];
	probGeno[GG][AA] = tmp2 * thesePi[A];
	probGeno[GG][CC] = tmp2 * thesePi[C];
	probGeno[GG][TT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[T];
	probGeno[TT][AA] = tmp2 * thesePi[A];
	probGeno[TT][CC] = tmp2 * thesePi[C];
	probGeno[TT][GG] = tmp2 * thesePi[G];

	//case ab/ab
	tmp = thesePhi[4] / K[4];
	tmp2 = tmp * thesePi[A];
	probGeno[AC][AC] = tmp2 * thesePi[C];
	probGeno[AG][AG] = tmp2 * thesePi[G];
	probGeno[AT][AT] = tmp2 * thesePi[T];
	tmp2 = tmp * thesePi[C];
	probGeno[CG][CG] = tmp2 * thesePi[G];
	probGeno[CT][CT] = tmp2 * thesePi[T];
	probGeno[GT][GT] = tmp * thesePi[G] * thesePi[T];

	//case ab/ac
	tmp = thesePhi[5] / K[5];
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[G];
	probGeno[AC][AG] = tmp2;
	probGeno[AC][CG] = tmp2;
	probGeno[AG][AC] = tmp2;
	probGeno[AG][CG] = tmp2;
	probGeno[CG][AC] = tmp2;
	probGeno[CG][AG] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[T];
	probGeno[AC][AT] = tmp2;
	probGeno[AC][CT] = tmp2;
	probGeno[AT][AC] = tmp2;
	probGeno[AT][CT] = tmp2;
	probGeno[CT][AC] = tmp2;
	probGeno[CT][AT] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[G] * thesePi[T];
	probGeno[AG][GT] = tmp2;
	probGeno[AG][AT] = tmp2;
	probGeno[AT][AG] = tmp2;
	probGeno[AT][GT] = tmp2;
	probGeno[GT][AG] = tmp2;
	probGeno[GT][AT] = tmp2;
	tmp2 = tmp * thesePi[C] * thesePi[G] * thesePi[T];
	probGeno[CG][CT] = tmp2;
	probGeno[CG][GT] = tmp2;
	probGeno[CT][CG] = tmp2;
	probGeno[CT][GT] = tmp2;
	probGeno[GT][CG] = tmp2;
	probGeno[GT][CT] = tmp2;

	//case aa/bc
	tmp = thesePhi[6] / K[6];
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[G];
	probGeno[AA][CG] = tmp2;
	probGeno[CC][AG] = tmp2;
	probGeno[GG][AC] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[T];
	probGeno[AA][CT] = tmp2;
	probGeno[CC][AT] = tmp2;
	probGeno[TT][AC] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[G] * thesePi[T];
	probGeno[AA][GT] = tmp2;
	probGeno[GG][AT] = tmp2;
	probGeno[TT][AG] = tmp2;
	tmp2 = tmp * thesePi[C] * thesePi[G] * thesePi[T];
	probGeno[CC][GT] = tmp2;
	probGeno[GG][CT] = tmp2;
	probGeno[TT][CG] = tmp2;

	//case ab/cc
	tmp = thesePhi[7] / K[7];
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[G];
	probGeno[AC][GG] = tmp2;
	probGeno[AG][CC] = tmp2;
	probGeno[CG][AA] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[C] * thesePi[T];
	probGeno[AC][TT] = tmp2;
	probGeno[AT][CC] = tmp2;
	probGeno[CT][AA] = tmp2;
	tmp2 = tmp * thesePi[A] * thesePi[G] * thesePi[T];
	probGeno[AG][TT] = tmp2;
	probGeno[AT][GG] = tmp2;
	probGeno[GT][AA] = tmp2;
	tmp2 = tmp * thesePi[C] * thesePi[G] * thesePi[T];
	probGeno[CG][TT] = tmp2;
	probGeno[CT][GG] = tmp2;
	probGeno[GT][CC] = tmp2;

	//case ab/cd
	tmp = thesePhi[8] / K[8];
	probGeno[AC][GT] = tmp;
	probGeno[AG][CT] = tmp;
	probGeno[AT][CG] = tmp;
	probGeno[CG][AT] = tmp;
	probGeno[CT][AG] = tmp;
	probGeno[GT][AC] = tmp;
}

bool TEMforDistanceEstimation::estimatePhiWithEM(std::vector<uint8_t*> & genoQual1, std::vector<uint8_t*> & genoQual2){
	//prepare estimates
	logfile->listFlush("Estimating initial base frequencies pi ...");
	guessPi(genoQual1, genoQual2);
	logfile->done();
	logfile->conclude("Initial pi are A=" + toString(pi.freq[0]) + ", C=" + toString(pi.freq[1]) + ", G=" + toString(pi.freq[2]) + " and T=" + toString(pi.freq[3]) + ".");
	logfile->listFlush("Estimating initial genotype classes phi ...");
	guessPhi(genoQual1, genoQual2);
	logfile->done();
	logfile->conclude("Initial phi are " + concatenateString(phi, 9, ", ") + ".");

	//variables
	double old_LL, LL = 0.0;
	double LL_diff;

	//now run EM
	logfile->startIndent("Estimating phi using an EM algorithm:");
	for(int iter=0; iter<maxNumEMIterations; ++iter){
		logfile->listFlush("Running EM iteration " + toString(iter+1) + " ...");
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
		std::vector<uint8_t*>::iterator it1 = genoQual1.begin();
		std::vector<uint8_t*>::iterator it2 = genoQual2.begin();
		for(; it1 != genoQual1.end(); ++it1, ++it2){
			//calculate P_G per site
			double sum = 0.0;
			for(int g1 = 0; g1<10; ++g1){
				for(int g2 = 0; g2<10; ++g2){
					P_G_one_site[g1][g2] = phredToLik[(*it1)[g1]] * phredToLik[(*it2)[g2]] * probGeno[g1][g2];

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
		for(int i=0; i<9; ++i)
			phi[i] = 0.0;

		double sum = 0.0;
		for(int g1 = 0; g1<10; ++g1){
			for(int g2 = 0; g2<10; ++g2){
				phi[genoToPhiMap(g1,g2)] += P_G[g1][g2];
				sum += P_G[g1][g2];
			}
		}

		for(int i=0; i<9; ++i)
			phi[i] /= sum;

		//update pi
		pi.clear();
		for(int g1 = 0; g1<10; ++g1){
			for(int g2 = 0; g2<10; ++g2){
				for(int b=0; b<4; ++b)
					if(genoToBaseMap(g1, g2, b)) pi.add((Base) b, P_G[g1][g2]);
			}
		}
		pi.normalize();

		//check if EM converged
		logfile->done();
		//throw "done!";
		if(iter > 0 ){
			LL_diff = LL - old_LL;
			logfile->conclude("LL = " + toString(LL) + " (deltaLL = " + toString(LL_diff) + ").");
			if(LL_diff < epsilonForEM){
				logfile->conclude("EM converged, delatLL = " + toString(LL_diff) + " < " + toString(epsilonForEM));
				distance = distanceObject->calculateDistance(phi);
				logfile->conclude("Resulting distance is " + toString(distance));
				logfile->endIndent();
				return true;
			}
		} else
			logfile->conclude("LL = " + toString(LL) + ".");
	}
	logfile->warning("EM reached maximum number of iterations (" + toString(maxNumEMIterations) + ") without converging!");
	distance = distanceObject->calculateDistance(phi);
	logfile->conclude("Resulting distance is " + toString(distance));
	logfile->endIndent();
	return false;
};

//----------------------------------------------------
//TDistanceEstimator
//----------------------------------------------------
TDistanceEstimator::TDistanceEstimator(TLog* Logfile, TParameters & params){
	logfile = Logfile;
	maxNumEMIterations = 0;
	epsilonForEM = 0.0;
	numGLFs = 0;
	glfs = NULL;
	readersOpened = false;

	//outputname
	outputName = params.getParameterStringWithDefault("out", "ATLAS");
	logfile->list("Writing output files with prefix '" + outputName + "'.");
}

void TDistanceEstimator::openGLF(TParameters & params){
	params.fillParameterIntoVector("glf", GLFNames, ',');
	numGLFs = GLFNames.size();
	if(numGLFs < 2)
		throw "At least two GLF files have to be provided to estimate distances!";

	//open files
	glfs = new TGlfReader[numGLFs];
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
void TDistanceEstimator::estimateDistances(TParameters & params){
	//open all GLF files specified
	openGLF(params);

	//open EM object
	TEMforDistanceEstimation EM_object(logfile, params);

	//in windows or whole genome?
	long windowLen = params.getParameterLongWithDefault("window", -1);
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
			logfile->startIndent("Estimating distance between individuals " + toString(g1+1) + " (" + GLFNames[g1] + ") and " + toString(g2+1) + " (" + GLFNames[g2] + "):");

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
}

bool TDistanceEstimator::moveToNextCommonChr(TGlfReader & g1, TGlfReader & g2){
	std::string chr1 = g1.chr();
	eraseAllOccurences(chr1,"chr");
	std::string chr2 = g2.chr();
	eraseAllOccurences(chr2,"chr");

	while(chr1 != chr2){
		//advance the one lagging behind
		if(g1.chrNumber() < g2.chrNumber()){
			if(!g1.jumpToNextChr()) return false;
		} else if(g1.chrNumber() > g2.chrNumber()){
			if(!g2.jumpToNextChr()) return false;
		} else
			throw "Different chromosomes in files " + g1.name() + "' and '" + g2.name() + "'!";

		chr1 = g1.chr();
		eraseAllOccurences(chr1,"chr");
		chr2 = g2.chr();
		eraseAllOccurences(chr2,"chr");
	}

	return true;
}


bool TDistanceEstimator::advance(TGlfReader & g1, TGlfReader & g2){
	//advance
	if(g2.position == g1.position){
		//advance both
		if(!g1.readNext()) return false;
		if(!g2.readNext()) return false;
	} else if(g2.position < g1.position){
		//advance g2
		if(!g2.readNext()) return false;
	} else {
		//advance g1
		if(!g1.readNext()) return false;
	}

	//make sure we are on same chromosome
	return(moveToNextCommonChr(g1, g2));
};

void TDistanceEstimator::readCommonSites(std::vector<uint8_t*> & genoQual1, std::vector<uint8_t*> & genoQual2, TGlfReader & g1, TGlfReader & g2){
	//parse GLFs. Only keep sites where both individuals have data!

	//rewind GLFs
	g1.rewind();
	g2.rewind();

	//if not both are good at least one file reach end. So we are done!
	while(advance(g1, g2)){

		//std::cout << "@\t" << g1.chr() << ":" << g1.position << "\t" << g2.chr() << ":" << g2.position << std::endl;

		if(g2.position == g1.position){
			//add data
			genoQual1.push_back(new uint8_t[10]);
			genoQual2.push_back(new uint8_t[10]);
			g1.fillGenotypeQualities(*genoQual1.rbegin());
			g2.fillGenotypeQualities(*genoQual2.rbegin());
		}
	}
};

void TDistanceEstimator::estimateDistanceGenomeWide(TEMforDistanceEstimation & EM_object, TGlfReader & g1, TGlfReader & g2, gz::ogzstream & out){
	//initialize storage for two windows
	logfile->listFlush("Reading common sites ...");
	std::vector<uint8_t*> genoQual1, genoQual2;
	readCommonSites(genoQual1, genoQual2, g1, g2);
	logfile->done();
	logfile->conclude("Read data for " + toString(genoQual1.size()) + " sites.");

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
	for(std::vector<uint8_t*>::iterator it=genoQual1.begin(); it!=genoQual1.end(); ++it)
		delete[] *it;
	for(std::vector<uint8_t*>::iterator it=genoQual2.begin(); it!=genoQual2.end(); ++it)
		delete[] *it;
	genoQual1.clear();
	genoQual2.clear();
	logfile->done();
};


//--------------------------------------------
// Estimation in windows
//--------------------------------------------
void TDistanceEstimator::estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, long windowLen){
	logfile->list("Will estimate genetic distance in windows of length " + toString(windowLen) + ".");
	if(windowLen < 100)
		throw "Window size must be at least 100bp!";

	//loop over all pairs
	for(int g1=0; g1<(numGLFs-1); ++g1){
		for(int g2 = g1+1; g2 < numGLFs; ++g2){
			logfile->startIndent("Estimating distance between individuals " + toString(g1+1) + " (" + GLFNames[g1] + ") and " + toString(g2+1) + " (" + GLFNames[g2] + "):");

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
}

void TDistanceEstimator::estimateDistanceInWindows(TEMforDistanceEstimation & EM_object, std::string filename, TGlfReader & g1, TGlfReader & g2, long windowLen){
	//initialize variables
	bool keepReading = true;
	bool isGood1 = true;
	bool isGood2 = true;

	//initialize storage for two windows
	//TODO: share across pairs? Do all pairs at once?
	std::vector<uint8_t*> genoQual1, genoQual2;
	genoQual1.resize(windowLen);
	genoQual2.resize(windowLen);
	for(int i=0; i<windowLen; ++i){
		genoQual1[i] = new uint8_t[10];
		genoQual2[i] = new uint8_t[10];
	}

	//open output file
	gz::ogzstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";

	//write header to output file
	out << "chr\twindowStart\twindowEnd\tnumSitesWithData\tfreqA\tfreqC\tfreqG\tfreqT\tfreq00_00\tfreq00_01\tfreq01_00\tfreq00_11\tfreq01_01\tfreq01_02\tfreq00_12\tfreq01_22\tfreq01_23\tgeneticDist\n";

	//prepare variables
	std::string curChr;
	long curChrLen;
	long windowStart;
	long windowEnd;

	int numSitesWithData = 100;

	//parse GLFs in windows
	logfile->startIndent("Will estimate distance in windows of size " + toString(windowLen) + ":");
	while(keepReading){
		//move to new chromosome
		curChr = g1.chr();
		curChrLen = g1.chrLength();
		windowStart = 1;
		windowEnd = windowLen+1;

		logfile->startNumbering("Chromosome " + curChr + ":");

		//parse all windows of chromosome
		while(windowStart < curChrLen){
			logfile->number("Window [" + toString(windowStart) + ", " + toString(windowEnd) + ")");
			logfile->addIndent();
			//read data
			isGood1 = g1.readNextWindow(genoQual1, curChr, windowStart, windowEnd);
			if(isGood1 || g1.eof()){
				isGood2 = g2.readNextWindow(genoQual2, curChr, windowStart, windowEnd);
				if(isGood2){
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

		if(g1.eof() && g2.eof())
			keepReading = false;
	}

	//clean up memory
	for(int i=0; i<windowLen; ++i){
		delete[] genoQual1[i];
		delete[] genoQual2[i];
	}

	logfile->endIndent();

}


//--------------------------------------------
// Writing estimates
//--------------------------------------------
void TDistanceEstimator::writeDistanceEstimates(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd, int numsitesWithData, TEMforDistanceEstimation & EM_object){
	out << chr << "\t" << windowStart << "\t" << windowEnd;
	writeDistanceEstimates(out, numsitesWithData, EM_object);
}

void TDistanceEstimator::writeDistanceEstimates(gz::ogzstream & out, int numsitesWithData, TEMforDistanceEstimation & EM_object){
	out << "\t" << numsitesWithData;
	//write pi
	for(int b=0; b<4; ++b)
		out << "\t" << EM_object.pi[b];
	//write phi
	for(int b=0; b<9; ++b)
		out << "\t" << EM_object.phi[b];
	//write distance
	out << "\t" << EM_object.distance;
	out << "\n";
}

void TDistanceEstimator::writeDistanceEstimatesNoData(gz::ogzstream & out, std::string & chr, long & windowStart, long & windowEnd){
	out << chr << "\t" << windowStart << "\t" << windowEnd << "\t";
	writeDistanceEstimatesNoData(out);
}

void TDistanceEstimator::writeDistanceEstimatesNoData(gz::ogzstream & out){
	out << "\t0";
	for(int i=0; i<14; ++i)
		out << "\t-";
	out << "\n";
}
















