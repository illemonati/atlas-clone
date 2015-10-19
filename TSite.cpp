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
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it)
		delete *it;
	bases.clear();
	hasData = false;
};

void TSiteDiploid::add(char & base, char & quality, int PosInRead, int pos5, int pos3, BaseContext & Context, int & ReadGroup){
	if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){
		if((int) quality > 33){
			if(base == 'A') bases.push_back(new TBaseDiploidA(quality, PosInRead, pos5, pos3, Context, ReadGroup));
			else if(base == 'C') bases.push_back(new TBaseDiploidC(quality, PosInRead, pos5, pos3, Context, ReadGroup));
			else if(base == 'G') bases.push_back(new TBaseDiploidG(quality, PosInRead, pos5, pos3, Context, ReadGroup));
			else bases.push_back(new TBaseDiploidT(quality, PosInRead, pos5, pos3, Context, ReadGroup));
		}
		hasData = true;
	}
};
void TSiteHaploid::add(char & base, char & quality, int PosInRead, int pos5, int pos3, BaseContext & Context, int & ReadGroup){
	if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){
		if((int) quality > 33){
			if(base == 'A') bases.push_back(new TBaseHaploidA(quality, PosInRead, pos5, pos3, Context, ReadGroup));
			else if(base == 'C') bases.push_back(new TBaseHaploidC(quality, PosInRead, pos5, pos3, Context, ReadGroup));
			else if(base == 'G') bases.push_back(new TBaseHaploidG(quality, PosInRead, pos5, pos3, Context, ReadGroup));
			else bases.push_back(new TBaseHaploidT(quality, PosInRead, pos5, pos3, Context, ReadGroup));
		}
		hasData = true;
	}
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
		int MLGenotype;
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

void TSiteDiploid::callMLEAllelePresence(TGenotypeMap & genoMap, gz::ogzstream & out){

	//TODO: not correct this way. Need Weighted sum, but where to get the weights from?

	if(hasData){
		//print coverage (and read bases)
		out << "\t" << bases.size();
		out << "\t" << getBases(); //printing data for debugging

		//calculate likelihoods for allele presence
		double likelihoods[4];
		likelihoods[0] = emissionProbabilities[0] + emissionProbabilities[1] + emissionProbabilities[2] + emissionProbabilities[3]; //A
		likelihoods[1] = emissionProbabilities[1] + emissionProbabilities[4] + emissionProbabilities[5] + emissionProbabilities[6]; //C
		likelihoods[2] = emissionProbabilities[2] + emissionProbabilities[5] + emissionProbabilities[7] + emissionProbabilities[8]; //G
		likelihoods[3] = emissionProbabilities[3] + emissionProbabilities[6] + emissionProbabilities[8] + emissionProbabilities[9]; //T

		//calculate phred-scaled likelihoods and find max
		//calcEmissionProbabilities(pmdObject);
		double maxGenotypeProb = 0.0;
		int MLAllele;
		for(int i=0; i<4; ++i){
			if(likelihoods[i] > maxGenotypeProb){
				MLAllele = i;
				maxGenotypeProb = likelihoods[i];
			}
		}

		//quality = 1 - prob
		double quality = 1.0 - likelihoods[MLAllele];

		//phred scale and normalize
		for(int i=0; i<4; ++i){
			likelihoods[i] = -10.0 * log10(likelihoods[i]);
		}

		//now print normalized (max = 0)
		for(int i=0; i<4; ++i){
			out << "\t" << round(likelihoods[i] - likelihoods[MLAllele]);
		}

		//add MLE genotype and quality = second smallest phred-scaled likelihood (like GATK)
		out << "\t" << genoMap.getGenotypeString(MLAllele);
		//out << "\t" << round(-10.0 * log10(quality));
		out << "\t" << quality;
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

