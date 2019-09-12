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
		bases.clear();
		hasData = false;
		referenceBase = 'N';
	}
};

void TSite::stealFromOther(TSite* other){
	//this function extracts all data from the other object and sets it to empty
	hasData = other->hasData;
	if(hasData){
		//copy data
		referenceBase = other->referenceBase;
		for(int i=0; i<numGenotypes; ++i){
			emissionProbabilities[i] = other->emissionProbabilities[i];
		}
		//copy pointers to bases, BUT NOT BASES
		for(std::vector<TBase*>::iterator it = other->bases.begin(); it!=other->bases.end(); ++it){
			bases.push_back(*it);
		}
		//remove pointers from other site
		other->bases.clear();
		other->hasData = false;
	}
}
void TSite::add(TBase* base){
	bases.push_back(base);
	hasData = true;
};


/*void TSite::add(Base & base, int & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup){
	if(base == A) bases.push_back(new TBase(base, quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == C) bases.push_back(new TBase(base, quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == G) bases.push_back(new TBase(base, quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else bases.push_back(new TBase(base, quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	hasData = true;
};*/

void TSite::addToBaseFrequencies(TBaseFrequencies & frequencies){
	if(hasData){
		static double weight = 1.0 / bases.size();
		for(baseIterator = bases.begin(); baseIterator!=bases.end(); ++baseIterator){
			(*baseIterator)->addToBaseFrequencies(frequencies, weight);
		}
	}
};

void TSite::calcEmissionProbabilities(double* vec){
	//do in log if coverage is high
	if(bases.size() < 50){
		//initialize
		for(int i=0; i<numGenotypes; ++i)
			vec[i] = 1.0;
		//multiply over emission probabilities of bases
		for(baseIterator = bases.begin(); baseIterator!=bases.end(); ++baseIterator)
			(*baseIterator)->addToEmissionProb(vec);

	} else {
		//initialize
		for(int i=0; i<numGenotypes; ++i)
			vec[i] = 0.0;
		//sum over log(emission probability) of bases
		for(baseIterator = bases.begin(); baseIterator!=bases.end(); ++baseIterator)
			(*baseIterator)->addToEmissionProbLog(vec);
		//now standardize before delog
		double max = vec[0];
		for(int i=1; i<numGenotypes; ++i)
			if(vec[i] > max) max = vec[i];
		for(int i=0; i<numGenotypes; ++i)
			vec[i] = exp(vec[i] - max);
	}
};

void TSite::calcEmissionProbabilities(){
	calcEmissionProbabilities(emissionProbabilities);
};

std::string TSite::getBases(){
	if(!hasData) return "-";
	std::string b = "";
	for(baseIterator = bases.begin(); baseIterator!=bases.end(); ++baseIterator)
		b += getBaseAsChar((*baseIterator)->getBaseAsEnum());
	return b;
};

int TSite::depth(){
	if(!hasData) return 0;
	return bases.size();
};

int TSite::refDepth(){
	if(!hasData) return 0;
	if(referenceBase == 'N') return 0;
	int counter = 0;
	for(unsigned int i=0; i<bases.size(); ++i){
		if(getBaseAsChar(bases[i]->getBaseAsEnum()) == referenceBase) ++counter;
	}
	return counter;
};

std::string TSite::getEmissionProbs(){
	std::string b;
	if(!hasData){
		b = "1";
		for(int i=1; i<numGenotypes; ++i){
			b += "\t1";
		}
	} else {
		b = toString(emissionProbabilities[0]);
		for(int i=1; i<numGenotypes; ++i){
			b += "\t" + toString(emissionProbabilities[i]);
		}
	}
	return b;
};

void TSite::calculateP_g(double* genotypeProbabilities, double* P_g){
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


void TSite::countAlleles(int* alleleCounts){
	alleleCounts[0] = 0;
	alleleCounts[1] = 0;
	alleleCounts[2] = 0;
	alleleCounts[3] = 0;

	for(TBase* it : bases)
		++alleleCounts[it->getBaseAsEnum()];
};

void TSite::countAllelesForImbalance(long**** siteImbalance){
	//calculate and return imbalance
	int b[4] = {0};
	for(TBase* it : bases){
		++b[it->getBaseAsEnum()];
	}
	++siteImbalance[b[0]][b[1]][b[2]][b[3]];
};

void TSite::countMates(int* mateCounts){
	mateCounts[0] = 0;
	mateCounts[1] = 0;

	for(TBase* it : bases)
		++mateCounts[it->isSecondMate];
};

void TSite::countFwdRev(int* frCounts){
	frCounts[0] = 0;
	frCounts[1] = 0;

	for(TBase* it : bases)
		++frCounts[it->isReverseStrand];
};

void TSite::contextStats(int** contextCounts, TQualityMap & qualMap){
	for(TBase* it : bases){
		int q = qualMap.errorToPhredInt(it->errorRate);
		int c = it->context;
		++contextCounts[q][c];
	}
};

void TSite::printPileup(gz::ogzstream & out){
	out << "\t" << referenceBase;
	out << "\t" << depth() << "\t" << refDepth();
	out << "\t" << getBases() << "\t" << getEmissionProbs();
}

void TSite::printPileupToScreen(){
	std::cout << "\t" << referenceBase;
	std::cout << "\t" << depth() << "\t" << refDepth();
	std::cout << "\t" << getBases() << "\t" << getEmissionProbs();
}


//-----------------------------------------------------------------------
// Genotype likelihoods
//-----------------------------------------------------------------------
void TSite::calculateDiploidPhredScaledGenotypeLikelihoods(double* phredGLs){
	if(hasData){
		int tmp;
		//calculate phred-scaled likelihoods and find max
		for(int i=0; i<numGenotypes; ++i){
			phredGLs[i] = makePhredByRef(emissionProbabilities[i]);
		}
	} else {
		for(int i=0; i<numGenotypes; ++i)
			phredGLs[i] = 0.0;
	}
};

void TSite::calculateNormalizedGenotypeLikelihoodsAndQuality(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled,  double & quality, double & maxGenotypeProb, int & MLGenotype){
	//calculate phred-scaled likelihoods and find max
	maxGenotypeProb = 100000.0;
	quality = 100000.0;
	std::vector<int> MLEs;
	std::vector<int>::iterator it;
	for(int i=0; i<numGenotypes; ++i){
		emissionProbabilitiesPhredScaled[i] = makePhredByRef(emissionProbabilities[i]);
		if(emissionProbabilitiesPhredScaled[i] < maxGenotypeProb){
			MLGenotype = i;
			quality = maxGenotypeProb;
			maxGenotypeProb = emissionProbabilitiesPhredScaled[i];
			MLEs.clear();
			MLEs.push_back(i);
		} else if(emissionProbabilitiesPhredScaled[i] == maxGenotypeProb){
			MLEs.push_back(i);
			quality = emissionProbabilitiesPhredScaled[i];
		} else if(emissionProbabilitiesPhredScaled[i] < quality){
			quality = emissionProbabilitiesPhredScaled[i];
		}
	}
	//select best allele at random if there are multiple options
	MLGenotype = MLEs[randomGenerator.pickOne(MLEs.size())];
	quality = quality - maxGenotypeProb;
};

double TSite::calculatePHomozygous(double* pGenotype){
	//calculate posterior probability for each genotype
	double postProb[numGenotypes];
	double tot = 0.0;

	for(int i=0; i<numGenotypes; ++i){
		postProb[i] = emissionProbabilities[i] * pGenotype[i];
		tot += postProb[i];
	}

	//make sum for all homozygous genotypes
	return (postProb[AA] + postProb[CC] + postProb[GG] + postProb[TT]) / tot;
};


























