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

void TSiteDiploid::add(char & base, char & quality, int PosInRead, int pos5, int pos3, BaseContext & Context, int & ReadGroup){
	if(base == 'A') bases.push_back(new TBaseDiploidA(quality, PosInRead, pos5, pos3, Context, ReadGroup));
	else if(base == 'C') bases.push_back(new TBaseDiploidC(quality, PosInRead, pos5, pos3, Context, ReadGroup));
	else if(base == 'G') bases.push_back(new TBaseDiploidG(quality, PosInRead, pos5, pos3, Context, ReadGroup));
	else bases.push_back(new TBaseDiploidT(quality, PosInRead, pos5, pos3, Context, ReadGroup));
	hasData = true;
};
void TSiteHaploid::add(char & base, char & quality, int PosInRead, int pos5, int pos3, BaseContext & Context, int & ReadGroup){
	if(base == 'A') bases.push_back(new TBaseHaploidA(quality, PosInRead, pos5, pos3, Context, ReadGroup));
	else if(base == 'C') bases.push_back(new TBaseHaploidC(quality, PosInRead, pos5, pos3, Context, ReadGroup));
	else if(base == 'G') bases.push_back(new TBaseHaploidG(quality, PosInRead, pos5, pos3, Context, ReadGroup));
	else bases.push_back(new TBaseHaploidT(quality, PosInRead, pos5, pos3, Context, ReadGroup));
	hasData = true;
};

void TSite::addToBaseFrequencies(TBaseFrequencies & frequencies){
	double weight = 1.0 / bases.size();
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		(*it)->addToBaseFrequencies(frequencies, weight);
	}
}

void TSite::calcEmissionProbabilities(){
	for(int i=0; i<numGenotypes; ++i){
		emissionProbabilities[i] = 1.0;
	}
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		for(int i=0; i<numGenotypes; ++i){
			emissionProbabilities[i] *= (*it)->getEmissionProbability(i);
		}
	}
}

void TSite::callMLEGenotype(TGenotypeMap & genoMap, gz::ogzstream & out){
	if(hasData){
		//print coverage (and read bases)
		out << "\t" << bases.size();
		//out << "\t" << getBases(); //printing data for debugging

		//calculate phred-scaled likelihoods and find max
		//calcEmissionProbabilities(pmdObject);
		double maxGenotypeProb = 100000.0;
		int MLGenotype = 0;
		double quality = 100000.0;
		for(int i=0; i<numGenotypes; ++i){
			emissionProbabilities[i] = -10.0 * log10(emissionProbabilities[i]);
			if(emissionProbabilities[i] < maxGenotypeProb){
				MLGenotype = i;
				quality = maxGenotypeProb;
				maxGenotypeProb = emissionProbabilities[i];
			} else if(emissionProbabilities[i] < quality){
				quality = emissionProbabilities[i];
			}
		}

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

void TSiteDiploid::callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, gz::ogzstream & out){
	if(hasData){
		//print coverage (and read bases)
		out << "\t" << bases.size();
		out << "\t" << getBases(); //printing data for debugging

		//calculate posterior probability for each genotype
		double postProb[numGenotypes];
		double tot = 0.0;
		for(int i=0; i<numGenotypes; ++i){
			postProb[i] += emissionProbabilities[i] * pGenotype[i];
			tot += postProb[i];
		}
		for(int i=0; i<numGenotypes; ++i){
			postProb[i] /= tot;
		}

		//make sums for different bases
		double postProbAllele[4];
		Genotype g;
		double maxProb = 0.0;
		int maxAllele = 0;
		for(int i=0; i<4; ++i){
			postProb[i] = 0.0;
			for(int j=0; j<4; ++j){
				g = genoMap.getGenotype(i, j);
				postProbAllele[i] += postProb[g];
			}
			if(postProb[i] > maxProb){
				maxAllele = i;
				maxProb = postProb[i];
			}
			//phred scale
			postProb[i] = -10.0 * log10(postProb[i]);
		}

		//quality = phred(1 - prob)
		double quality = round(-10.0 * log10(1.0 - maxProb));

		//now print
		for(int i=0; i<4; ++i){
			out << "\t" << round(postProb[i]);
		}

		//add MLE genotype and quality = second smallest phred-scaled likelihood (like GATK)
		out << "\t" << genoMap.getBaseAsChar(maxAllele);
		//out << "\t" << round(-10.0 * log10(quality));
		out << "\t" << quality << " -> " << maxProb;
	} else {
		out << "\t0\t-\t-\t-\t-\t-\t0";
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

