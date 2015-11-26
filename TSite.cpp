/*
 * TBase.cpp
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#include "TSite.h"

//-------------------------------------------------------
//TSite
//-------------------------------------------------------
void TSite::clear(){
	if(hasData){
		for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it)
			delete *it;
		bases.clear();
		hasData = false;
	}
};

void TSiteDiploid::add(char & base, char & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup){
	if(base == 'A') bases.push_back(new TBaseDiploidA(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == 'C') bases.push_back(new TBaseDiploidC(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == 'G') bases.push_back(new TBaseDiploidG(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else bases.push_back(new TBaseDiploidT(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	hasData = true;
};
void TSiteHaploid::add(char & base, char & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup){
	if(base == 'A') bases.push_back(new TBaseHaploidA(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == 'C') bases.push_back(new TBaseHaploidC(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == 'G') bases.push_back(new TBaseHaploidG(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else bases.push_back(new TBaseHaploidT(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	hasData = true;
};

void TSite::addToBaseFrequencies(TBaseFrequencies & frequencies){
	double weight = 1.0 / bases.size();
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		(*it)->addToBaseFrequencies(frequencies, weight);
	}
}

void TSite::calcEmissionProbabilities(){
	//assumes that emission probabilities were calculated for TBase!

	//do in log if coverage is high
	if(bases.size() < 50){
		for(int i=0; i<numGenotypes; ++i){
			emissionProbabilities[i] = 1.0;
		}
		for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
			for(int i=0; i<numGenotypes; ++i){
				emissionProbabilities[i] *= (*it)->getEmissionProbability(i);
			}
		}
	} else {
		for(int i=0; i<numGenotypes; ++i){
			emissionProbabilities[i] = 0.0;
		}
		for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
			for(int i=0; i<numGenotypes; ++i){
				emissionProbabilities[i] += log((*it)->getEmissionProbability(i));
			}
		}
		//now standardize before delog
		double max = emissionProbabilities[0];
		for(int i=1; i<numGenotypes; ++i){
			if(emissionProbabilities[i] > max) max = emissionProbabilities[i];
		}
		for(int i=0; i<numGenotypes; ++i){
			emissionProbabilities[i] = exp(emissionProbabilities[i] - max);
		}
	}
}

std::string TSite::getBases(){
	if(bases.size()==0) return "-";
	std::string b = "";
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		b += (*it)->getBase();
	}
	return b;
}

std::string TSite::getEmissionProbs(){
	std::string b = toString(emissionProbabilities[0]);
	for(int i=1; i<numGenotypes; ++i){
		b += "\t" + toString(emissionProbabilities[i]);
	}
	return b;
}

void TSite::calculateP_g(double* genotypeProbabilities){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<numGenotypes; ++i){
		P_g[i] =  emissionProbabilities[i] * genotypeProbabilities[i];
		sum += P_g[i];
	}
	for(int i=0; i<10; ++i){
		P_g[i] /= sum;
	}
}

double TSite::calculateWeightedSumOfEmissionProbs(double* weights){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<numGenotypes; ++i){
		sum += emissionProbabilities[i] * weights[i];
	}
	return sum;
}

double TSite::calculateLogLikelihood(double* genotypeProbabilities){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<numGenotypes; ++i){
		sum +=  emissionProbabilities[i] * genotypeProbabilities[i];
	}
	return log(sum);
}

//-----------------------------------------------------------------------
//MLE Callers
//-----------------------------------------------------------------------
void TSite::calculateNormalizedGenotypeLikelihoods(TRandomGenerator & randomGenerator, double & quality, double & maxGenotypeProb, int & MLGenotype){
	//calculate phred-scaled likelihoods and find max
	maxGenotypeProb = 100000.0;
	quality = 100000.0;
	std::vector<int> MLEs;
	for(int i=0; i<numGenotypes; ++i){
		emissionProbabilities[i] = makePhredByRef(emissionProbabilities[i]);
		if(emissionProbabilities[i] < maxGenotypeProb){
			MLGenotype = i;
			quality = maxGenotypeProb;
			maxGenotypeProb = emissionProbabilities[i];
			MLEs.clear();
			MLEs.push_back(i);
		} else if(emissionProbabilities[i] == maxGenotypeProb){
			MLEs.push_back(i);
			quality = emissionProbabilities[i];
		} else if(emissionProbabilities[i] < quality){
			quality = emissionProbabilities[i];
		}
	}

	//select best allele at random if there are multiple options
	MLGenotype = MLEs[randomGenerator.pickOne(MLEs.size())];
}

void TSite::callMLEGenotype(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool printRef){
	if(hasData){
		//print reference allele
		if(printRef) out << "\t" << referenceBase;

		//print coverage (and read bases)
		out << "\t" << bases.size();
		//out << "\t" << getBases(); //printing data for debugging

		//calc normalized likelihoods
		double quality, maxGenotypeProb;
		int MLGenotype;
		calculateNormalizedGenotypeLikelihoods(randomGenerator, quality, maxGenotypeProb, MLGenotype);

		//now print normalized (max = 0)
		for(int i=0; i<numGenotypes; ++i){
			out << "\t" << round(emissionProbabilities[i] - maxGenotypeProb);
		}

		//add MLE genotype and quality = second smallest phred-scaled likelihood (like GATK)
		out << "\t" << genoMap.getGenotypeString(MLGenotype);
		out << "\t" << round(quality - maxGenotypeProb);
	} else {
		out << "\t0";
		for(int i=0; i<numGenotypes; ++i) out << "\t-";
		out << "\t-\t0";
	}
}

void TSite::callMLEGenotypeVCF(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool printRef){
	if(hasData){
		//print reference allele
		out << "\t.\t" << referenceBase;
		//out << "\t(" << getBases() << ")"; //printing data for debugging

		//calc normalized likelihoods
		double quality, maxGenotypeProb;
		int MLGenotype;
		calculateNormalizedGenotypeLikelihoods(randomGenerator, quality, maxGenotypeProb, MLGenotype);

		//find alternative alleles
		std::string genoVCF;
		std::string PL;
		if(referenceBase != 'N') PL = toString(round(emissionProbabilities[genoMap.getGenotype(referenceBase, referenceBase)] - maxGenotypeProb)); //for PL field in VCF
		std::string geno = genoMap.getGenotypeString(MLGenotype);

		if(geno[0] != referenceBase){
			if(geno[1] != referenceBase){
				if(referenceBase != 'N'){
					PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(referenceBase, geno[0])] - maxGenotypeProb));
					PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(geno[0], geno[0])] - maxGenotypeProb));
				}
				if(geno[0] == geno[1]){
					out << "\t" << geno[0];
					genoVCF = "1/1";
				} else {
					out << "\t" << geno;
					genoVCF = "1/2";
					if(referenceBase != 'N'){
						PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(referenceBase, geno[1])] - maxGenotypeProb));
						PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(geno[0], geno[1])] - maxGenotypeProb));
						PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(geno[1], geno[1])] - maxGenotypeProb));
					}
				}
			} else {
				out << "\t" << geno[0];
				genoVCF = "0/1";
				if(referenceBase != 'N'){
					PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(referenceBase, geno[0])] - maxGenotypeProb));
					PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(geno[0], geno[0])] - maxGenotypeProb));
				}
			}
		} else if(geno[1] != referenceBase){
			out << "\t" << geno[1];
			genoVCF = "0/1";
			if(referenceBase != 'N'){
				PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(referenceBase, geno[1])] - maxGenotypeProb));
				PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(geno[1], geno[1])] - maxGenotypeProb));
			}
		} else {
			//both are ref -> let's find the second most likely genotype
			double maxPostProb = -1.0;
			std::vector<int> secondMostLikely;
			for(int i=0; i<numGenotypes; ++i){
				if(i != MLGenotype){
					if(emissionProbabilities[i] > maxPostProb){
						maxPostProb = emissionProbabilities[i];
						secondMostLikely.clear();
						secondMostLikely.push_back(i);
					} else if(emissionProbabilities[i] == maxPostProb){
						secondMostLikely.push_back(i);
					}
				}
			}
			//select best allele at random if there are multiple options
			std::string genoSecond = genoMap.getGenotypeString(secondMostLikely[randomGenerator.pickOne(secondMostLikely.size())]);

			//now use this one to decide on alternative allele
			if(genoSecond[0] != referenceBase){
				out << "\t" << genoSecond[0];
				if(referenceBase != 'N'){
					PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(referenceBase, genoSecond[0])] - maxGenotypeProb));
					PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(genoSecond[0], genoSecond[0])] - maxGenotypeProb));
				}
			} else {
				out << "\t" << genoSecond[1];
				if(referenceBase != 'N'){
					PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(referenceBase, genoSecond[1])] - maxGenotypeProb));
					PL +=  "," + toString(round(emissionProbabilities[genoMap.getGenotype(genoSecond[1], genoSecond[1])] - maxGenotypeProb));
				}
			}
			genoVCF = "0/0";
		}

		//print quality
		out << "\t" << round(quality);

		//print (no) filter
		out << "\t.";

		//print info fields: coverage and all normalized likelihoods
		out << "\tDP=" << bases.size();
		out << ";GG=" << round(emissionProbabilities[0] - maxGenotypeProb);
		for(int i=1; i<numGenotypes; ++i){
			out << "," << round(emissionProbabilities[i] - maxGenotypeProb);
		}

		//print format and genotype field
		if(referenceBase != 'N') out << "\tGT:PL\t" << genoVCF << ":" << PL;
		else out << "\tGT\t" << genoVCF;
	} else {
		out << "\t0";
		for(int i=0; i<numGenotypes; ++i) out << "\t-";
		out << "\t-\t0";
	}
}

void TSiteDiploid::calculatePhredScaledGenotypeLikelihoodsKnownAlleles(TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* phredEmissionProbs, double & quality, double & maxGenotypeProb, int & MLGenotype){
	//which genotypes?
	int genotypes[3];
	genotypes[0] = genoMap.getGenotype(referenceBase, referenceBase);
	genotypes[1] = genoMap.getGenotype(referenceBase, alt);
	genotypes[2] = genoMap.getGenotype(alt, alt);

	//calculate phred-scaled likelihoods and find max
	maxGenotypeProb = 100000.0;
	quality = 100000.0;
	std::vector<int> MLEs;

	for(int j=0; j<3; ++j){
		phredEmissionProbs[j] = makePhredByRef(emissionProbabilities[genotypes[j]]);
		if(phredEmissionProbs[j] < maxGenotypeProb){
			MLGenotype = j;
			quality = maxGenotypeProb;
			maxGenotypeProb = phredEmissionProbs[j];
			MLEs.clear();
			MLEs.push_back(j);
		} else if(phredEmissionProbs[j] == maxGenotypeProb){
			MLEs.push_back(j);
			quality = phredEmissionProbs[j];
		} else if(phredEmissionProbs[j] < quality){
			quality = phredEmissionProbs[j];
		}
	}

	//select best allele at random if there are multiple options
	MLGenotype = MLEs[randomGenerator.pickOne(MLEs.size())];
}

void TSiteDiploid::callMLEGenotypeKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){
	if(hasData){
		//print reference allele
		out << "\t" << referenceBase << "\t" << alt;

		//print coverage (and read bases)
		out << "\t" << bases.size();
		//out << "\t" << getBases(); //printing data for debugging

		//calc normalized likelihoods
		double quality, maxGenotypeProb;
		int MLGenotype;
		double phredEmissionProbs[3];
		calculatePhredScaledGenotypeLikelihoodsKnownAlleles(genoMap, alt, randomGenerator, phredEmissionProbs, quality, maxGenotypeProb, MLGenotype);

		//now print normalized (max = 0)
		for(int i=0; i<3; ++i){
			out << "\t" << round(phredEmissionProbs[i] - maxGenotypeProb);
		}

		//add MLE genotype and quality = second smallest phred-scaled likelihood (like GATK)
		out << "\t" << genoMap.getGenotypeString(MLGenotype);
		out << "\t" << round(quality - maxGenotypeProb);
	} else {
		out << "\t" << referenceBase << "\t" << alt;
		for(int i=0; i<3; ++i) out << "\t-";
		out << "\t-\t0";
	}
}


void TSiteDiploid::callMLEGenotypeVCFKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){
	if(hasData){
		//print reference allele
		out << "\t.\t" << referenceBase << "\t" << alt;

		//calc normalized likelihoods
		double quality, maxGenotypeProb;
		int MLGenotype;
		double phredEmissionProb[3];
		calculatePhredScaledGenotypeLikelihoodsKnownAlleles(genoMap, alt, randomGenerator, phredEmissionProb, quality, maxGenotypeProb, MLGenotype);

		//print quality
		out << "\t" << round(quality);

		//print (no) filter
		out << "\t.";

		//print info fields: coverage and all normalized likelihoods
		out << "\tDP=" << bases.size();

		//print format and genotype field
		out << "\tGT:PL\t";
		if(MLGenotype == 0) out << "0/0";
		else if(MLGenotype == 1) out << "0/1";
		else out << "1/1";
		out << ":" << round(phredEmissionProb[0] - maxGenotypeProb) << "," << round(phredEmissionProb[1] - maxGenotypeProb) << "," << round(phredEmissionProb[2] - maxGenotypeProb);
	} else {
		out << "\t0";
		for(int i=0; i<numGenotypes; ++i) out << "\t-";
		out << "\t-\t0";
	}
}

//-----------------------------------------------------------------------
//Bayesian Callers
//-----------------------------------------------------------------------
void TSite::calculateGenotypePosteriorProbabilities(double* pGenotype, TRandomGenerator & randomGenerator, double* postProb, int & MAP){
	double tot = 0.0;

	for(int i=0; i<numGenotypes; ++i){
		postProb[i] = emissionProbabilities[i] * pGenotype[i];
		tot += postProb[i];
	}

	double maxPostProb = -1.0;
	std::vector<int> MAPs;
	for(int i=0; i<numGenotypes; ++i){
		postProb[i] /= tot;
		if(postProb[i] > maxPostProb){
			maxPostProb = postProb[i];
			MAPs.clear();
			MAPs.push_back(i);
		} else if(postProb[i] == maxPostProb){
			MAPs.push_back(i);
		}
	}

	//select best allele at random if there are multiple options
	MAP = MAPs[randomGenerator.pickOne(MAPs.size())];
}

void TSite::callBayesianGenotype(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool printRef){
	if(hasData){
		//print reference allele
		if(printRef) out << "\t" << referenceBase;

		//print coverage (and read bases)
		out << "\t" << bases.size();
		//out << "\t" << getBases(); //printing data for debugging

		//calculate posterior probability for each genotype
		double postProb[numGenotypes];
		int MAPGenotype;
		calculateGenotypePosteriorProbabilities(pGenotype, randomGenerator, postProb, MAPGenotype);

		//print out phred-scaled posteriors
		for(int i=0; i<numGenotypes; ++i){
			out << "\t" << round(makePhredByRef(postProb[i]));
			//out << "\t" << postProb[i];
		}

		//add MAP genotype and quality
		out << "\t" << genoMap.getGenotypeString(MAPGenotype);
		out << "\t" << round(makePhred(1.0 - postProb[MAPGenotype]));
	} else {
		out << "\t0";
		for(int i=0; i<numGenotypes; ++i) out << "\t-";
		out << "\t-\t0";
	}
}

void TSite::callBayesianGenotypeVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out){
	if(hasData){
		//print reference allele
		out << "\t.\t" << referenceBase;
		//out << "\t(" << getBases() << ")"; //printing data for debugging

		//calculate posterior probability for each genotype
		double postProb[numGenotypes];
		int MAPGenotype;
		calculateGenotypePosteriorProbabilities(pGenotype, randomGenerator, postProb, MAPGenotype);

		//find alternative allele
		//if MAP genotype contains non ref allele, these are the alternatives
		std::string genoVCF;
		std::string GP;
		if(referenceBase != 'N') GP = toString(round(postProb[genoMap.getGenotype(referenceBase, referenceBase)])); //for GP field in VCF
		std::string geno = genoMap.getGenotypeString(MAPGenotype);
		if(geno[0] != referenceBase){
			if(geno[1] != referenceBase){
				if(referenceBase != 'N'){
					GP +=  "," + toString(round(postProb[genoMap.getGenotype(referenceBase, geno[0])]));
					GP +=  "," + toString(round(postProb[genoMap.getGenotype(geno[0], geno[0])]));
				}
				if(geno[0] == geno[1]){
					out << "\t" << geno[0];
					genoVCF = "1/1";
				} else {
					out << "\t" << geno;
					genoVCF = "1/2";
					if(referenceBase != 'N'){
						GP +=  "," + toString(round(postProb[genoMap.getGenotype(referenceBase, geno[1])]));
						GP +=  "," + toString(round(postProb[genoMap.getGenotype(geno[0], geno[1])]));
						GP +=  "," + toString(round(postProb[genoMap.getGenotype(geno[1], geno[1])]));
					}
				}
			} else {
				out << "\t" << geno[0];
				genoVCF = "0/1";
				if(referenceBase != 'N'){
					GP +=  "," + toString(round(postProb[genoMap.getGenotype(referenceBase, geno[0])]));
					GP +=  "," + toString(round(postProb[genoMap.getGenotype(geno[0], geno[0])]));
				}
			}
		} else if(geno[1] != referenceBase){
			out << "\t" << geno[1];
			genoVCF = "0/1";
			if(referenceBase != 'N'){
				GP +=  "," + toString(round(postProb[genoMap.getGenotype(referenceBase, geno[1])]));
				GP +=  "," + toString(round(postProb[genoMap.getGenotype(geno[1], geno[1])]));
			}
		} else {
			//both are ref -> let's find the second most likely genotype
			double maxPostProb = -1.0;
			std::vector<int> secondMostLikely;
			for(int i=0; i<numGenotypes; ++i){
				if(i != MAPGenotype){
					if(postProb[i] > maxPostProb){
						maxPostProb = postProb[i];
						secondMostLikely.clear();
						secondMostLikely.push_back(i);
					} else if(postProb[i] == maxPostProb){
						secondMostLikely.push_back(i);
					}
				}
			}
			//select best allele at random if there are multiple options
			std::string genoSecond = genoMap.getGenotypeString(secondMostLikely[randomGenerator.pickOne(secondMostLikely.size())]);

			//now use this one to decide on alternativ allele
			if(genoSecond[0] != referenceBase){
				out << "\t" << genoSecond[0];
				if(referenceBase != 'N'){
					GP +=  "," + toString(round(postProb[genoMap.getGenotype(referenceBase, geno[0])]));
					GP +=  "," + toString(round(postProb[genoMap.getGenotype(geno[0], geno[0])]));
				}
			} else {
				out << "\t" << genoSecond[1];
				if(referenceBase != 'N'){
					GP +=  "," + toString(round(postProb[genoMap.getGenotype(referenceBase, geno[1])]));
					GP +=  "," + toString(round(postProb[genoMap.getGenotype(geno[1], geno[1])]));
				}
			}
			genoVCF = "0/0";
		}

		//print quality
		out << "\t" << round(makePhred(1.0 - postProb[MAPGenotype]));

		//print (no) filter
		out << "\t.";

		//print info fields: coverage and all posterior probabilities
		out << "\tDP=" << bases.size() << ";PP=" << round(makePhred(postProb[0]));
		for(int i=1; i<numGenotypes; ++i){
			out << "," << round(makePhred(postProb[i]));
		}

		//print format and genotype field
		if(referenceBase != 'N') out << "\tGT:GP\t" << genoVCF << ":" << GP;
		else out << "\tGT\t" << genoVCF;
	} else {
		out << "\t.\t" << referenceBase << "\t.\t.\t.\tDP=0\tGT\t.";
	}
}

void TSiteDiploid::calculateGenotypePosteriorProbabilitiesKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* postProb, int & MAP){
	//which genotypes?
	int genotypes[3];
	genotypes[0] = genoMap.getGenotype(referenceBase, referenceBase);
	genotypes[1] = genoMap.getGenotype(referenceBase, alt);
	genotypes[2] = genoMap.getGenotype(alt, alt);

	double tot = 0.0;
	for(int i=0; i<3; ++i){
		postProb[i] = emissionProbabilities[genotypes[i]] * pGenotype[genotypes[i]];
		tot += postProb[i];
	}

	double maxPostProb = -1.0;
	std::vector<int> MAPs;
	for(int i=0; i<3; ++i){
		postProb[i] /= tot;
		if(postProb[i] > maxPostProb){
			maxPostProb = postProb[i];
			MAPs.clear();
			MAPs.push_back(i);
		} else if(postProb[i] == maxPostProb){
			MAPs.push_back(i);
		}
	}

	//select best allele at random if there are multiple options
	MAP = MAPs[randomGenerator.pickOne(MAPs.size())];
}


void TSiteDiploid::callBayesianGenotypeKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){
	if(hasData){
		//print reference allele
		out << "\t" << referenceBase << "\t" << alt;

		//print coverage (and read bases)
		out << "\t" << bases.size();
		//out << "\t" << getBases(); //printing data for debugging

		//calculate posterior probability for each genotype
		double postProb[3];
		int MAPGenotype;
		calculateGenotypePosteriorProbabilitiesKnownAlleles(pGenotype, genoMap, alt, randomGenerator, postProb, MAPGenotype);

		//print out phred-scaled posteriors
		for(int i=0; i<3; ++i){
			out << "\t" << round(makePhredByRef(postProb[i]));
		}

		//add MAP genotype and quality
		if(MAPGenotype == 0) out << "\t" << referenceBase << referenceBase;
		else if(MAPGenotype == 1) out << "\t" << referenceBase << alt;
		else out << "\t" << alt << alt;
		out << "\t" << round(makePhred(1.0 - postProb[MAPGenotype]));
	} else {
		out << "\t0";
		for(int i=0; i<3; ++i) out << "\t-";
		out << "\t-\t0";
	}
}

void TSiteDiploid::callBayesianGenotypeVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){
	//just consider known alleles: ref and alt
	if(hasData){
		//print reference allele
		out << "\t.\t" << referenceBase << "\t" << alt;
		//out << "\t(" << getBases() << ")"; //printing data for debugging

		//calculate posterior probability for the genotypes RR, AR and AA (R = ref, A = alt)
		double postProb[3];
		int MAPGenotype;
		calculateGenotypePosteriorProbabilitiesKnownAlleles(pGenotype, genoMap, alt, randomGenerator, postProb, MAPGenotype);

		//set genotype and GP string
		std::string genoVCF;
		if(MAPGenotype == 0) genoVCF = "0/0";
		else if(MAPGenotype == 1) genoVCF = "0/1";
		else genoVCF = "1/1";
		std::string GP =  toString(round(postProb[0])) + "," + toString(round(postProb[1])) + "," + toString(round(postProb[2]));

		//print quality
		out << "\t" << round(makePhred(1.0 - postProb[MAPGenotype]));

		//print (no) filter
		out << "\t.";

		//print info fields: coverage and all posterior probabilities
		out << "\tDP=" << bases.size();

		//print format and genotype field
		if(referenceBase != 'N') out << "\tGT:GP\t" << genoVCF << ":" << GP;
		else out << "\tGT\t" << genoVCF;
	} else {
		out << "\t.\t" << referenceBase << "\t" << alt << "\t.\t.\tDP=0\tGT\t.";
	}
}

//-----------------------------------------------------------------------
//Allele Presence Callers
//-----------------------------------------------------------------------
void TSiteDiploid::calculatePosteriorOnAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, double* postProbAllele, int & MAP){
	//calculate posterior probability for each genotype
	double postProb[numGenotypes];
	double tot = 0.0;

	for(int i=0; i<numGenotypes; ++i){
		postProb[i] = emissionProbabilities[i] * pGenotype[i];
		tot += postProb[i];
	}

	//make sums for different bases
	Genotype g;
	double maxProb = -1.0;
	std::vector<int> MAPs;
	for(int i=0; i<4; ++i){
		postProbAllele[i] = 0.0;
		for(int j=0; j<4; ++j){
			g = genoMap.getGenotype(i, j);
			postProbAllele[i] += postProb[g];
		}
		postProbAllele[i] /= tot;
		if(postProbAllele[i] > maxProb){
			maxProb = postProbAllele[i];
			MAPs.clear();
			MAPs.push_back(i);
		} else if(postProbAllele[i] == maxProb){
			MAPs.push_back(i);
		}
	}

	//select best allele at random if there are multiple options
	MAP = MAPs[randomGenerator.pickOne(MAPs.size())];
}

void TSiteDiploid::callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool printRef){
	if(hasData){
		//print ref base, coverage (and read bases)
		if(printRef) out << "\t" << referenceBase;
		out << "\t" << bases.size();
		//out << "\t" << getBases(); //printing data for debugging

		//calculate posterior probability for each genotype
		double postProbAllele[4];
		int MAPAllele;
		calculatePosteriorOnAllelePresence(pGenotype, genoMap, randomGenerator, postProbAllele, MAPAllele);

		//now print
		for(int i=0; i<4; ++i){
			out << "\t" << round(makePhredByRef(postProbAllele[i]));
		}

		//add chosen allele and quality = 1 - posterior probability
		out << "\t" << genoMap.getBaseAsChar(MAPAllele);
		out << "\t" << round(makePhred(1.0 - postProbAllele[MAPAllele]));
		//out << "\t" << quality << " -> " << maxProb;
	} else {
		out << "\t0\t-\t-\t-\t-\t-\t0";
	}
}

void TSiteDiploid::callAllelePresenceVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out){
	if(hasData){
		//print reference allele
		out << "\t.\t" << referenceBase;
		//out << "\t(" << getBases() << ")"; //printing data for debugging

		//calculate posterior probability for each genotype
		double postProbAllele[4];
		int MAPAllele;
		calculatePosteriorOnAllelePresence(pGenotype, genoMap, randomGenerator, postProbAllele, MAPAllele);

		//print alternative allele
		std::string genoVCF;
		char base = genoMap.getBaseAsChar(MAPAllele);

		if(base == referenceBase){
			//find second most likely base
			std::vector<int> secondBase;
			double maxProb = -1.0;
			for(int i=0; i<4; ++i){
				if(i != MAPAllele){
					if(postProbAllele[i] > maxProb){
						maxProb = postProbAllele[i];
						secondBase.clear();
						secondBase.push_back(i);
					} else if(postProbAllele[i] == maxProb){
						secondBase.push_back(i);
					}
				}
			}
			//select alternative allele at random if there are multiple options
			out << "\t" << genoMap.getBaseAsChar(secondBase[randomGenerator.pickOne(secondBase.size())]);
			genoVCF = "0|0";
		} else {
			out << "\t" << base;
			genoVCF = "1|1";
		}

		//print quality
		out << "\t" << round(makePhred(1.0 - postProbAllele[MAPAllele]));

		//print (no) filter
		out << "\t.";

		//print info fields: coverage and all posterior probabilities
		out << "\tDP=" << bases.size() << ";PP=" << round(makePhred(postProbAllele[0]));
		for(int i=1; i<4; ++i){
			out << "," << round(makePhred(postProbAllele[i]));
		}

		//print format and genotype field
		out << "\tGT\t" << genoVCF;

	} else {
		out << "\t.\t" << referenceBase << "\t.\t.\t.\tDP=0\tGT\t.";
	}
}

void TSiteDiploid::calculatePosteriorOnAllelePresenceKnownAlleles(double* pGenotype, char & alt, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, double* postProbAllele, int & MAP){
	//calculate posterior probability for each genotype
	double postProb[3];
	double tot = 0.0;

	//which genotypes?
	int genotypes[3];
	genotypes[0] = genoMap.getGenotype(referenceBase, referenceBase);
	genotypes[1] = genoMap.getGenotype(referenceBase, alt);
	genotypes[2] = genoMap.getGenotype(alt, alt);

	for(int i=0; i<3; ++i){
		postProb[i] = emissionProbabilities[genotypes[i]] * pGenotype[genotypes[i]];
		tot += postProb[i];
	}

	//standardize
	for(int i=0; i<3; ++i){
		postProb[i] /= tot;
	}

	//make sums for different bases
	postProbAllele[0] = postProb[0] + postProb[1]; // ref/ref and ref/alt
	postProbAllele[1] = postProb[1] + postProb[2]; // ref/alt and alt/alt

	if(postProbAllele[0] > postProbAllele[1]) MAP = 0;
	else if(postProbAllele[0] < postProbAllele[1]) MAP = 1;
	else MAP = randomGenerator.pickOne(2);
}

void TSiteDiploid::callAllelePresenceKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){
	if(hasData){
		//print ref base, alt base, coverage (and read bases)
		out << "\t" << referenceBase << "\t" << alt;
		out << "\t" << bases.size();
		//out << "\t" << getBases(); //printing data for debugging

		//calculate posterior probability for each possible
		double postProbAllele[2];
		int MAPAllele;
		calculatePosteriorOnAllelePresenceKnownAlleles(pGenotype, alt, genoMap, randomGenerator, postProbAllele, MAPAllele);

		//now print
		out << "\t" << round(makePhredByRef(postProbAllele[0])) << "\t" << round(makePhredByRef(postProbAllele[1])); //ref and then alt

		//add chosen allele and quality = 1 - posterior probability
		if(MAPAllele == 0) out << "\t" << referenceBase;
		else out << "\t" << alt;
		out << "\t" << round(makePhred(1.0 - postProbAllele[MAPAllele]));
	} else {
		out << "\t0\t-\t-\t-\t-\t-\t0";
	}
}

void TSiteDiploid::callAllelePresenceVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt){
	if(hasData){
		//print reference and alternative allele
		out << "\t.\t" << referenceBase << "\t" << alt;
		//out << "\t(" << getBases() << ")"; //printing data for debugging

		//calculate posterior probability for each genotype
		double postProbAllele[2];
		int MAPAllele;
		calculatePosteriorOnAllelePresenceKnownAlleles(pGenotype, alt, genoMap, randomGenerator, postProbAllele, MAPAllele);

		//print quality
		out << "\t" << round(makePhred(1.0 - postProbAllele[MAPAllele]));

		//print (no) filter
		out << "\t.";

		//print info fields: coverage and all posterior probabilities
		out << "\tDP=" << bases.size() << ";PP=" << round(makePhred(postProbAllele[0])) << "," << round(makePhred(postProbAllele[1]));

		//print chosen genotype
		std::string genoVCF;
		if(MAPAllele == 0) out << "\tGT\t0|0";
		else out << "\tGT\t1|1";
	} else {
		out << "\t.\t" << referenceBase << "\t" << alt << "\t.\t.\tDP=0\tGT\t.";
	}
}
