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

void TSiteDiploid::add(Base & base, int & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup){
	if(base == A) bases.push_back(new TBaseDiploidA(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == C) bases.push_back(new TBaseDiploidC(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == G) bases.push_back(new TBaseDiploidG(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else bases.push_back(new TBaseDiploidT(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	hasData = true;
};
void TSiteHaploid::add(Base & base, int & quality, int PosInRead, int PosInReadRev, double thisPMD_CT, double thisPMD_GA, BaseContext & Context, int & ReadGroup){
	if(base == A) bases.push_back(new TBaseHaploidA(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == C) bases.push_back(new TBaseHaploidC(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else if(base == G) bases.push_back(new TBaseHaploidG(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	else bases.push_back(new TBaseHaploidT(quality, PosInRead, PosInReadRev, thisPMD_CT, thisPMD_GA, Context, ReadGroup));
	hasData = true;
};

void TSite::addToBaseFrequencies(TBaseFrequencies & frequencies){
	if(hasData){
		static double weight = 1.0 / bases.size();
		for(baseIterator = bases.begin(); baseIterator!=bases.end(); ++baseIterator){
			(*baseIterator)->addToBaseFrequencies(frequencies, weight);
		}
	}
};

void TSite::calcEmissionProbabilities(double* vec){
	//assumes that emission probabilities were calculated for TBase!

	//do in log if coverage is high
	if(bases.size() < 50){
		for(int i=0; i<numGenotypes; ++i){
			vec[i] = 1.0;
		}
		for(baseIterator = bases.begin(); baseIterator!=bases.end(); ++baseIterator){
			for(int i=0; i<numGenotypes; ++i){
				vec[i] *= (*baseIterator)->getEmissionProbability(i);
			}
		}
	} else {
		for(int i=0; i<numGenotypes; ++i){
			vec[i] = 0.0;
		}
		for(baseIterator = bases.begin(); baseIterator!=bases.end(); ++baseIterator){
			for(int i=0; i<numGenotypes; ++i){
				vec[i] += log((*baseIterator)->getEmissionProbability(i));
			}
		}
		//now standardize before delog
		double max = vec[0];
		for(int i=1; i<numGenotypes; ++i){
			if(vec[i] > max) max = vec[i];
		}
		for(int i=0; i<numGenotypes; ++i){
			vec[i] = exp(vec[i] - max);
		}
	}
}

void TSite::calcEmissionProbabilities(){
	calcEmissionProbabilities(emissionProbabilities);
}

std::string TSite::getBases(){
	if(!hasData) return "-";
	std::string b = "";
	for(baseIterator = bases.begin(); baseIterator!=bases.end(); ++baseIterator){
		b += (*baseIterator)->getBase();
	}
	return b;
}

int TSite::depth(){
	if(!hasData) return 0;
	return bases.size();
};

int TSite::refDepth(){
	if(!hasData) return 0;
	if(referenceBase == 'N') return 0;
	int counter = 0;
	for(int i=0; i<bases.size(); ++i){
		if(bases[i]->getBase() == referenceBase) ++counter;
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
}

void TSite::calculateP_g(double* & genotypeProbabilities, double* & P_g){
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


void TSite::countAlleles(long**** siteImbalance){
	//calculate and return imbalance
	int b[4] = {0};
	for(std::vector<TBase*>::iterator it = bases.begin(); it!=bases.end(); ++it){
		++b[(*it)->getBaseAsEnum()];
	}
	++siteImbalance[b[0]][b[1]][b[2]][b[3]];
}

void TSite::printPileup(gz::ogzstream & out){
	out << "\t" << referenceBase;
	out << "\t" << depth() << "\t" << refDepth();
	out << "\t" << getBases() << "\t" << getEmissionProbs();
}


//-----------------------------------------------------------------------
//MLE Callers
//-----------------------------------------------------------------------
void TSite::calculateNormalizedGenotypeLikelihoods(uint8_t* normalizedGL, uint32_t & maxLL){
	if(hasData){
		int tmp;
		//calculate phred-scaled likelihoods and find max
		double maxGenotypeProb = 100000.0;
		double* emissionProbabilitiesPhredScaled = new double[numGenotypes];
		for(int i=0; i<numGenotypes; ++i){
			emissionProbabilitiesPhredScaled[i] = makePhredByRef(emissionProbabilities[i]);
			if(emissionProbabilitiesPhredScaled[i] < maxGenotypeProb)
				maxGenotypeProb = emissionProbabilitiesPhredScaled[i];
		}
		for(int i=0; i<numGenotypes; ++i){
			tmp = round(emissionProbabilitiesPhredScaled[i] - maxGenotypeProb);
			if(tmp > 255) tmp = 255;
			normalizedGL[i] = tmp;
		}
		delete[] emissionProbabilitiesPhredScaled;
		maxLL = round(maxGenotypeProb);
	} else {
		for(int i=0; i<numGenotypes; ++i)
			normalizedGL[i] = 0;
		maxLL = 0;
	}
}

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
}

void TSite::findSecondMostLikelyGenotype(TRandomGenerator & randomGenerator, double* emissionProbabilitiesPhredScaled, TGenotypeMap & genoMap, int MLGenotype, std::string & genoSecond){
	double maxPostProb = 1000.0;
	std::vector<int> secondMostLikely;
	for(int i=0; i<numGenotypes; ++i){
		if(i != MLGenotype){
			if(emissionProbabilitiesPhredScaled[i] < maxPostProb){
				maxPostProb = emissionProbabilitiesPhredScaled[i];
				secondMostLikely.clear();
				secondMostLikely.push_back(i);
			} else if(emissionProbabilitiesPhredScaled[i] == maxPostProb){
				secondMostLikely.push_back(i);
			}
		}
	}
	//select best allele at random if there are multiple options
	genoSecond = genoMap.getGenotypeString(secondMostLikely[randomGenerator.pickOne(secondMostLikely.size())]);
}

void TSite::findSecondMostProbableGenotype(TRandomGenerator & randomGenerator, double* postProb, TGenotypeMap & genoMap, int MAPGenotype, std::string & genoSecond){
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
	genoSecond = genoMap.getGenotypeString(secondMostLikely[randomGenerator.pickOne(secondMostLikely.size())]);
}

int TSite::getAltAlleleBasedOnSecondMostLikelyGenotype(std::string & geno, std::string & genoSecond, TRandomGenerator & randomGenerator){
	int altAllele;
	if(genoSecond[0] != geno[0]){
		if(genoSecond[1] != geno[0]){
			if(randomGenerator.getRand() < 0.5) altAllele = 0;
			else altAllele = 1;
		} else altAllele = 0;
	} else altAllele = 1;
	return altAllele;
}

void TSite::callMLEGenotype(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out){
	out << "\t" << referenceBase;

	if(hasData){
		//print reference allele

		//print coverage (and read bases)
		out << "\t" << bases.size();
		//out << "\t" << getBases(); //printing data for debugging

		//calc normalized likelihoods
		double quality, maxGenotypeProb;
		int MLGenotype;
		double* emissionProbabilitiesPhredScaled = new double[numGenotypes];
		calculateNormalizedGenotypeLikelihoodsAndQuality(randomGenerator, emissionProbabilitiesPhredScaled, quality, maxGenotypeProb, MLGenotype);

		//now print normalized (max = 0)
		for(int i=0; i<numGenotypes; ++i){
			out << "\t" << round(emissionProbabilitiesPhredScaled[i] - maxGenotypeProb);
		}

		//add MLE genotype and quality = second smallest phred-scaled likelihood (like GATK)
		out << "\t" << genoMap.getGenotypeString(MLGenotype);
		out << "\t" << round(quality);
		delete[] emissionProbabilitiesPhredScaled;
	} else {
		out << "\t0";
		for(int i=0; i<numGenotypes; ++i) out << "\t-";
		out << "\t-\t0";
	}

}

void TSite::callMLEGenotypeVCF(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool gVCF, bool noAltIfHomoRef, std::string & basesString){
	//if you have alleles R, A, B, C then the order of the PL is: RR, RA, AA | RB, AB, BB | RC, AC, BC, CC

	if(hasData){
		//print reference allele
		out << "\t.\t" << referenceBase;

		//allelic depth
		int R_AD=0, A_AD=0, B_AD=0, C_AD=0;

		//calc normalized likelihoods
		double quality, maxGenotypeProb;
		int MLGenotype;
		double* emissionProbabilitiesPhredScaled = new double[numGenotypes];
		calculateNormalizedGenotypeLikelihoodsAndQuality(randomGenerator, emissionProbabilitiesPhredScaled, quality, maxGenotypeProb, MLGenotype);

		//find alternative alleles
		std::string genoVCF;
		std::string PL, AD, AI, rest;
		if(referenceBase != 'N') PL = toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, referenceBase)] - maxGenotypeProb)); //for PL field in VCF
		std::string geno = genoMap.getGenotypeString(MLGenotype);

		if(geno[0] != referenceBase){
			if(geno[1] != referenceBase){
				if(referenceBase != 'N'){
					PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, geno[0])] - maxGenotypeProb));
					PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(geno[0], geno[0])] - maxGenotypeProb));
				}

				if(geno[0] == geno[1]){
					out << "\t" << geno[0];
					genoVCF = "1/1";

					//calculate AD for R, A and B
					for(unsigned int i=0; i<basesString.size(); ++i){
						if(basesString[i] == referenceBase) ++R_AD;
						else if(basesString[i] == geno[0]) ++A_AD;
						else rest += basesString[i];
					}
					AD = toString(R_AD) + ',' + toString(A_AD);

					if(gVCF){
						out << ",<NON_REF>";
						std::string genoSecond = "";
						findSecondMostLikelyGenotype(randomGenerator, emissionProbabilitiesPhredScaled, genoMap, MLGenotype, genoSecond);

						//now use second most likely genotype one to decide on alternative allele
						int altAllele = getAltAlleleBasedOnSecondMostLikelyGenotype(geno, genoSecond, randomGenerator);

						//additional PL
						if(referenceBase != 'N'){
							PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, genoSecond[altAllele])] - maxGenotypeProb));
							PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(geno[0], genoSecond[altAllele])] - maxGenotypeProb));
							PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(genoSecond[altAllele], genoSecond[altAllele])] - maxGenotypeProb));
						}

						//additional AD
						if(rest.size() == 0) AD += ",0";
						else {
							for(unsigned int i=0; i<rest.size(); ++i) if (rest[i] == genoSecond[altAllele]) ++B_AD;
							AD += ',' + toString(B_AD);
						}
					}
				}
				else {
					out << "\t" << geno[0] << ',' << geno[1];
					genoVCF = "1/2";

					if(referenceBase != 'N'){
						PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, geno[1])] - maxGenotypeProb));
						PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(geno[0], geno[1])] - maxGenotypeProb));
						PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(geno[1], geno[1])] - maxGenotypeProb));
					}

					//calculate AD for R, A and B
					for(unsigned int i=0; i<basesString.size(); ++i){
						if(basesString[i] == referenceBase) ++R_AD;
						else if(basesString[i] == geno[0]) ++A_AD;
						else if(basesString[i] == geno[1]) ++B_AD;
						else rest += basesString[i];
					}
					AD = toString(R_AD) + ',' + toString(A_AD) + ',' + toString(B_AD);

					if(gVCF){
						out << ",<NON_REF>";
						std::string genoSecond = "";
						findSecondMostLikelyGenotype(randomGenerator, emissionProbabilitiesPhredScaled, genoMap, MLGenotype, genoSecond);

						//now use second most likely genotype one to decide on alternative allele
						int altAllele = getAltAlleleBasedOnSecondMostLikelyGenotype(geno, genoSecond, randomGenerator);

						//additional PL
						if(referenceBase != 'N'){
							PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, genoSecond[altAllele])] - maxGenotypeProb));
							PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(geno[0], genoSecond[altAllele])] - maxGenotypeProb));
							PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(genoSecond[1], genoSecond[altAllele])] - maxGenotypeProb));
							PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(genoSecond[altAllele], genoSecond[altAllele])] - maxGenotypeProb));
						}

						//additional AD
						if(rest.size() == 0) AD += ",0";
						else {
							for(unsigned int i=0; i<rest.size(); ++i)	if (rest[i] == genoSecond[altAllele]) ++C_AD;
							AD += ',' + toString(C_AD);
						}
					}
				}
			} else { //geno[1]=ref
				out << "\t" << geno[0];
				genoVCF = "0/1";
				if(referenceBase != 'N'){
					PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, geno[0])] - maxGenotypeProb));
					PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(geno[0], geno[0])] - maxGenotypeProb));
				}

				//calculate AD for R, A and B
				for(unsigned int i=0; i<basesString.size(); ++i){
					if(basesString[i] == referenceBase) ++R_AD;
					else if(basesString[i] == geno[0]) ++A_AD;
					else rest += basesString[i];
				}
				AD = toString(R_AD) + ',' + toString(A_AD);

				if(gVCF){
					out << ",<NON_REF>";
					std::string genoSecond = "";
					findSecondMostLikelyGenotype(randomGenerator, emissionProbabilitiesPhredScaled, genoMap, MLGenotype, genoSecond);

					//now use second most likely genotype one to decide on alternative allele
					int altAllele = getAltAlleleBasedOnSecondMostLikelyGenotype(geno, genoSecond, randomGenerator);

					//additional PL
					if(referenceBase != 'N'){
						PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, genoSecond[altAllele])] - maxGenotypeProb));
						PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(geno[0], genoSecond[altAllele])] - maxGenotypeProb));
						PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(genoSecond[altAllele], genoSecond[altAllele])] - maxGenotypeProb));
					}

					//additional AD
					if(rest.size() == 0) AD += ",0";
					else{
						for(unsigned int i=0; i<rest.size(); ++i) if (rest[i] == genoSecond[altAllele]) ++B_AD;
						AD += ',' + toString(B_AD);
					}
				}
			}
		} else if(geno[1] != referenceBase){  //geno[0]=ref
			out << "\t" << geno[1];
			genoVCF = "0/1";
			if(referenceBase != 'N'){
				PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, geno[1])] - maxGenotypeProb));
				PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(geno[1], geno[1])] - maxGenotypeProb));
			}

			//calculate AD for R, A and B
			for(unsigned int i=0; i<basesString.size(); ++i){
				if(basesString[i] == referenceBase) ++R_AD;
				else if(basesString[i] == geno[1]) ++A_AD;
				else rest += basesString[i];
			}
			AD = toString(R_AD) + ',' + toString(A_AD);

			if(gVCF){
				out << ",<NON_REF>";
				std::string genoSecond = "";
				findSecondMostLikelyGenotype(randomGenerator, emissionProbabilitiesPhredScaled, genoMap, MLGenotype, genoSecond);

				//now use second most likely genotype one to decide on alternative allele
				int altAllele = getAltAlleleBasedOnSecondMostLikelyGenotype(geno, genoSecond, randomGenerator);

				//additional PL
				if(referenceBase != 'N'){
					PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, genoSecond[altAllele])] - maxGenotypeProb));
					PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(geno[1], genoSecond[altAllele])] - maxGenotypeProb));
					PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(genoSecond[altAllele], genoSecond[altAllele])] - maxGenotypeProb));
				}

				//additional AD
				if(rest.size() == 0) AD += ",0";
				else{
					for(unsigned int i=0; i<rest.size(); ++i)	if (rest[i] == genoSecond[altAllele]) ++B_AD;
					AD += ',' + toString(B_AD);
				}
			}

		} else {
			//both are ref -> let's find the second most likely genotype
			genoVCF = "0/0";
			std::string genoSecond = "";
			findSecondMostLikelyGenotype(randomGenerator, emissionProbabilitiesPhredScaled, genoMap, MLGenotype, genoSecond);

			//now use second most likely genotype one to decide on alternative allele
			int altAllele = getAltAlleleBasedOnSecondMostLikelyGenotype(geno, genoSecond, randomGenerator);

			if(referenceBase != 'N' && !noAltIfHomoRef){
				PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(referenceBase, genoSecond[altAllele])] - maxGenotypeProb));
				PL +=  "," + toString(round(emissionProbabilitiesPhredScaled[genoMap.getGenotype(genoSecond[altAllele], genoSecond[altAllele])] - maxGenotypeProb));
			}

			//calculate AD for R and A
			for(unsigned int i=0; i<basesString.size(); ++i){
				if(basesString[i] == referenceBase) ++R_AD;
				else if(basesString[i] == genoSecond[altAllele]) ++A_AD;
			}
			AD = toString(R_AD);
			if (!noAltIfHomoRef) AD += ',' + toString(A_AD);

			//what to print in alt field
			if(gVCF) out << "\t<NON_REF>";
			else if (noAltIfHomoRef) out << "\t."; //for programs like vcf-tools that are confused by alt alleles when homozygous ref
			else {
				if(genoSecond[0] != referenceBase) out << "\t" << genoSecond[0];
				else out << "\t" << genoSecond[1];
			}
		}

		//prepare to calculate AI for info field
		std::vector<double> vec;
		char del=',';
		fillVectorFromString(AD,vec,del);
		AI = toString(randomGenerator.binomPValue(vec[0], vec[1]));
		//old GATK version
	//	AI = toString(vec[0] / bases.size());
	//	for(unsigned int i=1; i<vec.size(); ++i){
	//		AI += ',' + toString(vec[i] / bases.size());
	//	}
		//print (no) variant quality and (no) filter
		out << "\t.\t.";

		//print info fields: coverage
		out << "\tDP=" << bases.size(); // << ";bases=" << basesString;

		//print format and genotype and all normalized likelihoodsfield
		if(genoVCF.size() == 0){
			std::cout << "ref=" << referenceBase << " "<< geno[0] << geno[1] << std::endl;
		}
		if(referenceBase != 'N'){
			if(!gVCF){
				out << "\tGT:AD:AI:DP:GQ:PL:GG\t" <<  genoVCF << ':' << AD << ":" << AI << ":" <<  bases.size() << ":" << round(quality) << ':' << PL << ':'<< round(emissionProbabilitiesPhredScaled[0] - maxGenotypeProb);
				for(int i=1; i<numGenotypes; ++i){
					out << "," << round(emissionProbabilitiesPhredScaled[i] - maxGenotypeProb);
				}
			} else 	out << "\tGT:AD:AI:DP:GQ:PL\t" <<  genoVCF << ':' << AD << ":" << AI << ":" <<  bases.size() << ":" << round(quality) << ':' << PL;
		}
		else{
			out << "\tGT:DP:GQ\t" << genoVCF << ":" <<  bases.size() << ':' << round(quality);
		}
		delete[] emissionProbabilitiesPhredScaled;
	} else {
		//hasData is false
		//if(gVCF) out << "\t.\t" << referenceBase << "\t.\t.\t.\t.\tGT:DP\t./.:0";
		out << "\t.\t" << referenceBase << "\t.\t.\t.\t.\tGT:DP:GQ\t./.:0:0";
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

void TSiteDiploid::calculateGenotypeLikelihoodsKnownAlleles(TGenotypeMap & genoMap, char & alt, TRandomGenerator & randomGenerator, double* emissionProbs, double & sumEmissionProbs, int & pos){
	//which genotypes?
	int genotypes[3];
	genotypes[0] = genoMap.getGenotype(referenceBase, referenceBase);
	genotypes[1] = genoMap.getGenotype(referenceBase, alt);
	genotypes[2] = genoMap.getGenotype(alt, alt);

	//calculate likelihoods and their sum

	for(int j=0; j<3; ++j){
		emissionProbs[j] = emissionProbabilities[genotypes[j]];
		sumEmissionProbs += emissionProbs[j];
	}
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
		out << "\t" << genoMap.getGenotypeStringKnownAlleles(MLGenotype, referenceBase, alt);
		out << "\t" << round(quality - maxGenotypeProb);
	} else {
		out << "\t" << referenceBase << "\t" << alt << "\t" << 0;
		for(int i=0; i<3; ++i) out << "\t-";
		out << "\t-\t0";
	}
}

void TSiteDiploid::callMLEGenotypeKnownAllelesBeagle(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, std::string & chr, int & pos, long & start, bool & printOnlyGL){
	//print reference allele
	if(hasData){
		if(!printOnlyGL) out << chr << "_" << pos + start + 1 << "\t" << referenceBase << "\t" << alt;
		//calc normalized likelihoods
		double sumEmissionProbs = 0;
		double emissionProbs[3];
		calculateGenotypeLikelihoodsKnownAlleles(genoMap, alt, randomGenerator, emissionProbs, sumEmissionProbs, pos);

		//now print normalized (max = 0)
		out << std::setprecision(6);
		out << emissionProbs[0] / sumEmissionProbs;
		out << "\t" << emissionProbs[1] / sumEmissionProbs;
		out << "\t" << emissionProbs[2] / sumEmissionProbs;
	}
	else if(!printOnlyGL) out << chr << "_" << pos + start + 1 << "\tN\tN\t0.333\t0.333\t0.333";
	else out << "\t0.333\t0.333\t0.333";
	out << "\n";
}


void TSiteDiploid::callMLEGenotypeVCFKnownAlleles(TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, bool noAltIfHomoRef, std::string & basesString){
	if(hasData){
		//calc normalized likelihoods
		double quality, maxGenotypeProb;
		int MLGenotype, R_AD=0, A_AD=0;
		double phredEmissionProb[3];
		calculatePhredScaledGenotypeLikelihoodsKnownAlleles(genoMap, alt, randomGenerator, phredEmissionProb, quality, maxGenotypeProb, MLGenotype);
		//calculate AD
		for(unsigned int i=0; i<basesString.size(); ++i){
			if(basesString[i] == referenceBase) ++R_AD;
			else if(basesString[i] == alt) ++A_AD;
		}

		//print reference and alt allele
		out << "\t.\t" << referenceBase;
		if(noAltIfHomoRef && MLGenotype == 0) out << "\t.";
		else out << "\t" << alt;

		//print (no) variant quality and (no) filter
		out << "\t.\t.";

		//print info fields: coverage and all normalized likelihoods
		out << "\tDP=" << bases.size();

		//print format and genotype field
		out << "\tGT:AD:DP:GQ:PL\t";
		if(MLGenotype == 0) out << "0/0";
		else if(MLGenotype == 1) out << "0/1";
		else out << "1/1";
		out << ":" << toString(R_AD) << "," << toString(A_AD) << ":" <<  bases.size() << ":" << round(quality) << ':' << round(phredEmissionProb[0] - maxGenotypeProb) << "," << round(phredEmissionProb[1] - maxGenotypeProb) << "," << round(phredEmissionProb[2] - maxGenotypeProb);
	} else {
		out << "\t.\t" << referenceBase << "\t.\t.\t.\t.\tGT:DP:GQ\t./.:0:0";
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

void TSite::callBayesianGenotype(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out){
	//print reference allele
	out << "\t" << referenceBase;
	if(hasData){
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

void TSite::callBayesianGenotypeVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool noAltIfHomoRef, bool printPP, bool onlyPhredGP, std::string & basesString){
	if(hasData){
		//variables for allelic depth
		int R_AD=0, A_AD=0;
		std::string AD, AI, rest;

		//print reference allele
		out << "\t.\t" << referenceBase;
		//out << "\t(" << getBases() << ")"; //printing data for debugging

		//calculate posterior probabilities
		double postProb[numGenotypes];
		int MAPGenotype;
		calculateGenotypePosteriorProbabilities(pGenotype, randomGenerator, postProb, MAPGenotype);

//		//calc normalized likelihoods
//		double quality, maxGenotypeProb;
//		int MLGenotype;
//		double* emissionProbabilitiesPhredScaled = new double[numGenotypes];
//		calculateNormalizedGenotypeLikelihoodsAndQuality(randomGenerator, emissionProbabilitiesPhredScaled, quality, maxGenotypeProb, MLGenotype);

		//find alternative allele
		//if MAP genotype contains non ref allele, these are the alternatives
		std::string genoVCF;
		std::string GP;
		if(referenceBase != 'N'){
			if(onlyPhredGP)
				GP = toString(makePhred(postProb[genoMap.getGenotype(referenceBase, referenceBase)])); //for GP field in VCF
			else
				GP = toString(round(makePhred(1 - postProb[genoMap.getGenotype(referenceBase, referenceBase)]))); //for GP field in VCF
		}
		std::string geno = genoMap.getGenotypeString(MAPGenotype);
		if(geno[0] != referenceBase){
			if(geno[1] != referenceBase){
				if(referenceBase != 'N'){
					if(onlyPhredGP){
						GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(referenceBase, geno[0])]));
						GP +=  "," + toString(makePhred(1 - postProb[genoMap.getGenotype(geno[0], geno[0])]));
					} else {
						GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(referenceBase, geno[0])])));
						GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(geno[0], geno[0])])));
					}
				}
				if(geno[0] == geno[1]){
					out << "\t" << geno[0];
					genoVCF = "1/1";
				} else {
					out << "\t" << geno[0] << ',' << geno[1];
					genoVCF = "1/2";
					if(referenceBase != 'N'){
						if(onlyPhredGP){
							GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(referenceBase, geno[1])]));
							GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(geno[0], geno[1])]));
							GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(geno[1], geno[1])]));
						} else {
							GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(referenceBase, geno[1])])));
							GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(geno[0], geno[1])])));
							GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(geno[1], geno[1])])));
						}
					}
				}
				//calculate AD for R and A
				for(unsigned int i=0; i<basesString.size(); ++i){
					if(basesString[i] == referenceBase) ++R_AD;
					else if(basesString[i] == geno[0]) ++A_AD;
					else rest += basesString[i];

				}
				AD = toString(R_AD) + ',' + toString(A_AD);


			} else {
				out << "\t" << geno[0];
				genoVCF = "0/1";
				if(referenceBase != 'N'){
					if(onlyPhredGP){
						GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(referenceBase, geno[0])]));
						GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(geno[0], geno[0])]));
					} else {
						GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(referenceBase, geno[0])])));
						GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(geno[0], geno[0])])));
					}
				}

				//calculate AD for R and A
				for(unsigned int i=0; i<basesString.size(); ++i){
					if(basesString[i] == referenceBase) ++R_AD;
					else if(basesString[i] == geno[0]) ++A_AD;

				}
				AD = toString(R_AD) + ',' + toString(A_AD);
			}
		} else if(geno[1] != referenceBase){
			out << "\t" << geno[1];
			genoVCF = "0/1";
			if(referenceBase != 'N'){
				if(onlyPhredGP){
					GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(referenceBase, geno[1])]));
					GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(geno[1], geno[1])]));
				} else {
					GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(referenceBase, geno[1])])));
					GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(geno[1], geno[1])])));
				}
			}
			//calculate AD for R and A
			for(unsigned int i=0; i<basesString.size(); ++i){
				if(basesString[i] == referenceBase) ++R_AD;
				else if(basesString[i] == geno[1]) ++A_AD;
			}
			AD = toString(R_AD) + ',' + toString(A_AD);


		} else {
			//both are ref -> let's find the second most likely genotype and alternative allele
			std::string genoSecond;
			findSecondMostProbableGenotype(randomGenerator, postProb, genoMap, MAPGenotype, genoSecond);
			int altAllele = getAltAlleleBasedOnSecondMostLikelyGenotype(geno, genoSecond, randomGenerator);

			//calculate AD for R and A
			for(unsigned int i=0; i<basesString.size(); ++i){
				if(basesString[i] == referenceBase) ++R_AD;
				else if(basesString[i] == genoSecond[altAllele]) ++A_AD;
			}
			AD = toString(R_AD);
			if (!noAltIfHomoRef) AD += ',' + toString(A_AD);

			//write alternative allele
			if(noAltIfHomoRef)
				out << "\t.";
			else
				out << "\t" << genoSecond[altAllele];

			//add to GP
			if(referenceBase != 'N' && !noAltIfHomoRef){
				if(onlyPhredGP){
					GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(referenceBase, genoSecond[altAllele])]));
					GP +=  "," + toString(makePhred(postProb[genoMap.getGenotype(genoSecond[altAllele], genoSecond[altAllele])]));
				} else {
					GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(referenceBase, genoSecond[altAllele])])));
					GP +=  "," + toString(round(makePhred(1 - postProb[genoMap.getGenotype(genoSecond[altAllele], genoSecond[altAllele])])));

				}
			}
			genoVCF = "0/0";
		}

		//print (no) variant quality and (no) filter
		out << "\t.\t.";

		//print info fields: coverage
		out << "\tDP=" << bases.size();
		if(printPP){
			//all posterior probabilities
			out << ";PP=" << round(makePhred(1 - postProb[0]));
			for(int i=1; i<numGenotypes; ++i)
				out << "," << round(makePhred(1 - postProb[i]));
		}

//		//all likelihoods
//		out << ";GG=" << round(emissionProbabilitiesPhredScaled[0] - maxGenotypeProb);
//		for(int i=1; i<numGenotypes; ++i)
//			out << "," << round(emissionProbabilitiesPhredScaled[i] - maxGenotypeProb);

		//print format and genotype field
		if(referenceBase != 'N') out << "\tGT:DP:AD:GQ:GP\t" <<  genoVCF << ":" <<  bases.size()<< ":" << AD << ":" << round(makePhred(1.0 - postProb[MAPGenotype])) << ':' << GP;
		else out << "\tGT:DP:GQ\t" << genoVCF << ":" <<  bases.size() << ':' << round(makePhred(1.0 - postProb[MAPGenotype]));
	} else {
		out << "\t.\t" << referenceBase << "\t.\t.\t.\tDP=0\tGT:DP:GQ\t./.:0:0";
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
		out << "\t" << referenceBase;
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

void TSiteDiploid::callAllelePresence(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out){
	out << "\t" << referenceBase;
	if(hasData){
		//print ref base, coverage (and read bases)
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
	} else 	out << "\t0\t-\t-\t-\t-\t-\t0";
}

void TSiteDiploid::callAllelePresenceVCF(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, bool noAltIfHomoRef, std::string basesString){
	if(hasData){
		//print reference allele
		out << "\t.\t" << referenceBase;
		//out << "\t(" << getBases() << ")"; //printing data for debugging

		//calculate posterior probability for each genotype
		double postProbAllele[4];
		int MAPAllele, R_AD=0, A_AD=0;
		char alt = 'X';
		calculatePosteriorOnAllelePresence(pGenotype, genoMap, randomGenerator, postProbAllele, MAPAllele);

		//print alternative allele
		std::string genoVCF;
		char base = genoMap.getBaseAsChar(MAPAllele);

		if(base == referenceBase){
			genoVCF = "0";
			if(noAltIfHomoRef) out << "\t.";
			else{
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
				alt = genoMap.getBaseAsChar(secondBase[randomGenerator.pickOne(secondBase.size())]);
				out << "\t" << alt;
			}
		} else {
			out << "\t" << base;
			genoVCF = "1";
			alt=base;
		}
		//calculate AD
		for(unsigned int i=0; i<basesString.size(); ++i){
			if(basesString[i] == referenceBase) ++R_AD;
			else if(basesString[i] == alt) ++A_AD;
		}

		//print (no) quality
		out << "\t.";
	//	std::cout << postProbAllele[MAPAllele] << " " << 1.0 - postProbAllele[MAPAllele] << " "<< makePhred(1.0-postProbAllele[MAPAllele]) << std::endl;

		//print (no) filter
		out << "\t.";

		//print (no) info
		if(bases.size() > 0) out << "\tDP=" << bases.size();
		else out << "\t.";

		//print format field and genotype, coverage and posterior probabilities field
		out << "\tGT:AD:DP:GQ:PP\t" << genoVCF << ":" << R_AD;
		if(!(noAltIfHomoRef && genoVCF == "0")) out << "," << A_AD;
		out << ":" << bases.size() << ":" << round(makePhred(1.0 - postProbAllele[MAPAllele])) << ":" <<round(makePhred(postProbAllele[0]));
		for(int i=1; i<4; ++i){
			out << "," << round(makePhred(postProbAllele[i]));
		}
	} else {
		out << "\t.\t" << referenceBase << "\t.\t.\t.\t.\tGT:DP:PP\t.:0:0,0,0,0";
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
		out << "\t" << referenceBase;
		out << "\t0\t-\t-\t-\t-\t-\t0";
	}
}

void TSiteDiploid::callAllelePresenceVCFKnownAlleles(double* pGenotype, TGenotypeMap & genoMap, TRandomGenerator & randomGenerator, gz::ogzstream & out, char & alt, bool noAltIfHomoRef, std::string basesString){
	if(hasData){
		//print reference and alternative allele
		out << "\t.\t" << referenceBase;
		//out << "\t(" << getBases() << ")"; //printing data for debugging

		//calculate posterior probability for each genotype
		double postProbAllele[2];
		int MAPAllele, R_AD=0, A_AD=0;
		calculatePosteriorOnAllelePresenceKnownAlleles(pGenotype, alt, genoMap, randomGenerator, postProbAllele, MAPAllele);

		if(noAltIfHomoRef && MAPAllele ==0 ) out << "\t.";
		else out << "\t" << alt;

		//calculate AD
		for(unsigned int i=0; i<basesString.size(); ++i){
			if(basesString[i] == referenceBase) ++R_AD;
			else if(basesString[i] == alt) ++A_AD;
		}
		//print (no) quality
		out << "\t.";

		//print (no) filter
		out << "\t.";

		//print (no) info fields: coverage
		if(bases.size() > 0) out << "\tDP=" << bases.size();
		else out << "\t.";

		//print chosen genotype and coverage and all posterior probabilities
		std::string genoVCF;
		if(MAPAllele == 0) out << "\tGT:AD:DP:GQ:PP\t0:" << R_AD << "," << A_AD << ":" << bases.size() << ':' << round(makePhred(1.0 - postProbAllele[0])) << ":" << round(makePhred(postProbAllele[0])) << "," << round(makePhred(postProbAllele[1]));
		else out << "\tGT:AD:DP:GQ:PP\t1:" << R_AD << "," << A_AD << ":" << bases.size() << ':' << round(makePhred(1.0 - postProbAllele[MAPAllele])) << ":" << round(makePhred(postProbAllele[0])) << "," << round(makePhred(postProbAllele[1]));
	} else {
		out << "\t.\t" << referenceBase << "\t" << "." << "\t.\t.\t.\tGT:DP\t.:0";
	}
}

void TSiteDiploid::callRandomBase(TRandomGenerator & randomGenerator, gz::ogzstream & out){
	if(hasData){
		//print ref base, alt base, coverage (and read bases)
		out << "\t" << referenceBase << "\t" << bases.size();
		out << "\t";
		for(unsigned int i = 0; i<bases.size(); ++i){
			out << bases[i]->getBase();
		}
		out << "\t" << bases[randomGenerator.pickOne(bases.size())]->getBase();
 	} else {
		out << "\t" << referenceBase << "\t0\t-\t-";
	}
}

void TSiteDiploid::majorityCall(TRandomGenerator & randomGenerator, gz::ogzstream & out){
	if(hasData){
		//print ref base, alt base, coverage (and read bases)
		out << "\t" << referenceBase << "\t" << bases.size();
		out << "\t";

		//count bases
		int counts[5] = {0};
		for(unsigned int i = 0; i<bases.size(); ++i){
			out << bases[i]->getBase();
			++counts[bases[i]->getBaseAsEnum()];
		}

		//find majority
		char b[5] = {'A','C','G','T','N'};
		int max = 0;
		std::vector<char> maxBase;
		for(int i=0; i<5; ++i){
			if(counts[i] > max){
				max = counts[i];
				maxBase.clear();
				maxBase.push_back(b[i]);
			} else if(counts[i] == max){
				maxBase.push_back(b[i]);
			}
		}
		out << "\t" << maxBase[randomGenerator.pickOne(maxBase.size())];
 	} else {
		out << "\t" << referenceBase << "\t0\t-\t-";
	}
}


double TSiteDiploid::calculatePHomozygous(double* pGenotype){
	//calculate posterior probability for each genotype
	double postProb[numGenotypes];
	double tot = 0.0;

	for(int i=0; i<numGenotypes; ++i){
		postProb[i] = emissionProbabilities[i] * pGenotype[i];
		tot += postProb[i];
	}

	//make sum for all homozygous genotypes
	return (postProb[AA] + postProb[CC] + postProb[GG] + postProb[TT]) / tot;
}


void TSiteHaploid::addToExpectedBaseCounts(TBaseFrequencies & baseFreq, double* expectedCounts){
	double* tmp = new double[4];
	for(int b=0; b<4; ++b) tmp[b]=0.0;
	for(std::vector<TBase*>::iterator it = bases.begin(); it != bases.end(); ++it){
		(*it)->addToExpectedBaseCounts(baseFreq, tmp);
	}
	for(int b=0; b<4; ++b) expectedCounts[b] += tmp[b] / (double) bases.size();
}

void TSiteHaploid::calculatePoolFreqLikelihoods(int & numChromosomes, TGenotypeMap & genoMap, Base & allele1, Base & allele2, gz::ogzstream & out){
	//write coverage
	out << "\t" << genoMap.getBaseAsChar(allele1) << "\t" << genoMap.getBaseAsChar(allele2) << "\t" << bases.size();

	if(hasData){
		//calculate likelihood for sample frequencies from 0 to num chromosomes for allele 1
		double LL;
		double f;
		for(int y = 0; y < (numChromosomes + 1); ++y){
			//calculate likelihood
			LL = 0.0;
			f = (double) y / (double) numChromosomes;
			for(std::vector<TBase*>::iterator it = bases.begin(); it != bases.end(); ++it){
				LL += log((*it)->getEmissionProbability(allele1) * f + (*it)->getEmissionProbability(allele2) * (1.0 -f));
			}

			//write it to file
			out << "\t" << LL;
		}
	} else {
		for(int y = 0; y < (numChromosomes + 1); ++y){
			out << "\t0.0";
		}
	}
}



























