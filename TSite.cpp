/*
 * TBase.cpp
 *
 *  Created on: May 9, 2015
 *      Author: wegmannd
 */

#include "TSite.h"


//-------------------------------------------------------
//TBase
//-------------------------------------------------------
void TBaseA::fillEmissionProbabilities(TPMD* pmdCT, TPMD* pmdGA){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;
	double pmd = pmdGA->getProb(pos3);

	emissionProbabilities.set(AA, oneMinusError);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - pmd) * errorOneThird + pmd * oneMinusError);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, 0.5 - errorOneThird);
	emissionProbabilities.set(AG, ((1.0 + pmd) * oneMinusError + (1.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, (pmd * oneMinusError + (2.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(CG));
};

void TBaseC::fillEmissionProbabilities(TPMD* pmdCT, TPMD* pmdGA){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;
	double pmd = pmdCT->getProb(pos5);

	emissionProbabilities.set(AA , errorOneThird);
	emissionProbabilities.set(CC , (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(GG , errorOneThird);
	emissionProbabilities.set(TT , errorOneThird);
	emissionProbabilities.set(AC , ((1.0 - pmd) * oneMinusError + (1.0 + pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG , errorOneThird);
	emissionProbabilities.set(AT , errorOneThird);
	emissionProbabilities.set(CG , emissionProbabilities.get(AC));
	emissionProbabilities.set(CT , emissionProbabilities.get(AC));
	emissionProbabilities.set(GT , errorOneThird);
};

void TBaseG::fillEmissionProbabilities(TPMD* pmdCT, TPMD* pmdGA){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;
	double pmd = pmdGA->getProb(pos3);

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, errorOneThird);
	emissionProbabilities.set(GG, (1.0 - pmd) * oneMinusError + pmd * errorOneThird);
	emissionProbabilities.set(TT, errorOneThird);
	emissionProbabilities.set(AC, errorOneThird);
	emissionProbabilities.set(AG, ((1.0 - pmd) * oneMinusError + (1.0 + pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AT, errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AG));
	emissionProbabilities.set(CT, errorOneThird);
	emissionProbabilities.set(GT, emissionProbabilities.get(AG));

};

void TBaseT::fillEmissionProbabilities(TPMD* pmdCT, TPMD* pmdGA){
	//pre-calculate all emission probabilities given the error rate and the distance from either end of the read
	double errorOneThird = errorRate / 3.0;
	double oneMinusError = 1.0 - errorRate;
	double pmd = pmdCT->getProb(pos5);

	emissionProbabilities.set(AA, errorOneThird);
	emissionProbabilities.set(CC, (1.0 - pmd) * errorOneThird + pmd * oneMinusError);
	emissionProbabilities.set(GG, errorOneThird);
	emissionProbabilities.set(TT, oneMinusError);
	emissionProbabilities.set(AC, (pmd * oneMinusError + (2.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(AG, errorOneThird);
	emissionProbabilities.set(AT, 0.5 - errorOneThird);
	emissionProbabilities.set(CG, emissionProbabilities.get(AC));
	emissionProbabilities.set(CT, ((1.0 + pmd) * oneMinusError + (1.0 - pmd) * errorOneThird) / 2.0);
	emissionProbabilities.set(GT, 0.5 - errorOneThird);
};

//-------------------------------------------------------
//TPMD
//-------------------------------------------------------
void TPMD::initializeFunction(std::string & pmdString){
	//parse string to get model.  options are
	// none
	// Skoglund[lambda,c]
	// Veeramah[a,b,c]
	std::string example = "Use either Skoglund[p,c] or Veeramah[a,b,c]";

	//check if it is none
	if(pmdString == "none"){
		myFunction = new TPMDFunction();
	} else {
		std::string::size_type pos = pmdString.find_first_of('[');
		if(pos == std::string::npos) throw "Can not initialize post mortem damage function: wrong format!\n" + example;
		std::string name = pmdString.substr(0,pos);

		//prepare first value
		std::string tmp = pmdString.substr(pos+1, pmdString.length() - pos - 1);
		pos = tmp.find_first_of(',');
		if(pos == std::string::npos) throw "Can not initialize post mortem damage function: wrong format!\n" + example;
		double first = atof(tmp.substr(0, pos).c_str());

		//switch between functions
		if(name == "Skoglund"){
			//get lambda and
			if(first < 0.0) throw "Can not initialize Skoglund function with lambda < 0!";
			if(first > 1.0) throw "Can not initialize Skoglund function with lambda > 1!";
			double c = atof(tmp.substr(pos+1).c_str());
			if(c < 0.0) throw "Can not initialize Skoglund function with c < 0!";
			myFunction = new TPMDSkoglund(first, c);
		} else if(name == "Veeramah"){
			//get a, b and c
			if(first < 0.0) throw "Can not initialize Veeramah function with a < 0!";

			//get b
			tmp = tmp.substr(pos+1);
			pos = tmp.find_first_of(',');
			if(pos == std::string::npos) throw "Can not initialize post mortem damage function: wrong format!\n" + example;
			double b = atof(tmp.substr(0, pos).c_str());
			if(b < 0.0) throw "Can not initialize Veeramah function with b < 0!";

			//get c
			double c = atof(tmp.substr(pos+1).c_str());
			if(c < 0.0) throw "Can not initialize Veeramah function with c < 0!";

			//test if it can be > 1
			if(first + c > 1) throw "Can not initialize Veeramah function with a + c > 1!";

			//initialze
			myFunction = new TPMDVeeramah(first, b, c);
		} else throw "Can not initialize post mortem damage function: wrong name!\n" + example;
	}
	functionInitialized = true;
}

//-------------------------------------------------------
//TSite
//-------------------------------------------------------
void TSite::add(char & base, char & quality, int pos5, int pos3, TPMD* pmdCT, TPMD* pmdGA){
	if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){
		double error = pow(10.0, - ((double)quality - 33.0)/10.0);
		if(error < 1.0){
			if(base == 'A') bases.push_back(new TBaseA(error, pos5, pos3, pmdCT, pmdGA));
			else if(base == 'C') bases.push_back(new TBaseC(error, pos5, pos3, pmdCT, pmdGA));
			else if(base == 'G') bases.push_back(new TBaseG(error, pos5, pos3, pmdCT, pmdGA));
			else bases.push_back(new TBaseT(error, pos5, pos3, pmdCT, pmdGA));
		}
		hasData = true;
	}
};

std::string TSite::getBases(){
	if(bases.size()==0) return "-";
	std::string b = "";
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		b += (*it)->getBase();
	}
	return b;
}

void TSite::addToBaseFrequencies(TBaseFrequencies & frequencies){
	double weight = 1.0 / bases.size();
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		(*it)->addToBaseFrequencies(frequencies, weight);
	}
}

void TSite::calcEmissionProbabilities(){
	//calculate normalized genotype probabilities according to Bayes rule
	for(int i=0; i<10; ++i){
		emissionProbabilities[i] = 1.0;
		for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
			emissionProbabilities[i] *= (*it)->getEmissionProbability(i);
		}
	}
}

std::string TSite::getEmissionProbs(){
	std::string b = toString(emissionProbabilities[0]);
	for(int i=1; i<10; ++i){
		b += "\t" + toString(emissionProbabilities[i]);
	}
	return b;
}

void TSite::calculateP_g(double* genotypeProbabilities){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<10; ++i){
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
	for(int i=0; i<10; ++i){
		sum += emissionProbabilities[i] * weights[i];
	}
	return sum;
}

double TSite::calculateLogLikelihood(double* genotypeProbabilities){
	//calculate normalized genotype probabilities according to Bayes rule
	double sum = 0.0;
	for(int i=0; i<10; ++i){
		sum +=  emissionProbabilities[i] * genotypeProbabilities[i];
	}
	return log(sum);
}

