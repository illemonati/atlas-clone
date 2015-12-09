/*
 * TPostMortemDamage.cpp
 *
 *  Created on: Nov 27, 2015
 *      Author: wegmannd
 */


/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#include "TPostMortemDamage.h"

//---------------------------------------------------------------
//TPMDTable
//---------------------------------------------------------------
void TPMDTable::calculateSums(){
	if(!sumsCalculated){
		sums = new long*[maxLength];
		for(int p=0; p<maxLength; ++p){
			sums[p] = new long[4];
			for(int i=0; i<4; ++i){
				sums[p][i] = counts[p][i][0];
				for(int j=1; j<4; ++j){
					sums[p][i] += counts[p][i][j];
				}
			}
		}
		sumsCalculated = true;
	}
};

void TPMDTable::deleteSums(){
	for(int p=0; p<maxLength; ++p){
		delete[] sums[p];
	}
	delete[] sums;
	sumsCalculated = false;
}

TPMDTable::TPMDTable(int MaxLength){
	maxLength = MaxLength;
	sumsCalculated = false;
	sums = NULL;

	//init counts
	counts = new long**[maxLength];
	for(int p=0; p<maxLength; ++p){
		counts[p] = new long*[5];
		for(int i=0; i<5; ++i){
			counts[p][i] = new long[5];
			for(int j=0; j<5; ++j){
				counts[p][i][j] = 0;
			}
		}
	}
};

TPMDTable::~TPMDTable(){
	for(int p=0; p<maxLength; ++p){
		for(int i=0; i<5; ++i){
			delete[] counts[p][i];
		}
		delete[] counts[p];
	}
	delete counts;
	if(sumsCalculated) deleteSums();
};

void TPMDTable::empty(){
	for(int p=0; p<maxLength; ++p){
		for(int i=0; i<5; ++i){
			for(int j=0; j<5; ++j){
				counts[p][i][j] = 0;
			}
		}
	}
};

void TPMDTable::add(int pos, Base ref, Base read){
	if(pos < maxLength)
		++counts[pos][ref][read];
};

void TPMDTable::writeTable(std::ofstream & out, std::string prefix){
	calculateSums();

	//now print out
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			out <<	prefix << genoMap.getBaseAsChar(i) << "->" << genoMap.getBaseAsChar(j);
			for(int p=0; p<maxLength; ++p){
				if(sums[p][i] < 1) out << "\t0.0";
				else out << "\t" << (double) counts[p][i][j] / (double) sums[p][i];
			}
			out << "\n";
		}
	}
};

void TPMDTable::writeTableWithCounts(std::ofstream & out, std::string prefix){
	for(int i=0; i<4; ++i){
		for(int j=0; j<4; ++j){
			out <<	prefix << genoMap.getBaseAsChar(i) << "->" << genoMap.getBaseAsChar(j);
			for(int p=0; p<maxLength; ++p){
				out << "\t" << (double) counts[p][i][j];
			}
			out << "\n";
		}
	}
};

std::string TPMDTable::getPMDStringCT(){
	calculateSums();
	std::string s = "Empiric[";
	double tmpCT, tmpTC;  //tmpRefRead
	for(int p=0; p<maxLength; ++p){
		if(p>0) s += ",";
		if(sums[p][T] < 1 || sums[p][C] < 1) s += "0.0";
		else {
			tmpCT = (double) counts[p][C][T] / (double) sums[p][C];
			tmpTC = (double) counts[p][T][C] / (double) sums[p][T];
			s += toString(std::max(0.0, (tmpCT - tmpTC)));
		}
	}
	return s + "]";
};

std::string TPMDTable::getPMDStringGA(){
	calculateSums();
	std::string s = "Empiric[";
	double tmpGA, tmpAG;  //tmpRefRead
	for(int p=0; p<maxLength; ++p){
		if(p>0) s += ",";
		if(sums[p][A] < 1 || sums[p][G] < 1) s += "0.0";
		else {
			tmpGA = (double) counts[p][G][A] / (double) sums[p][G];
			tmpAG = (double) counts[p][A][G] / (double) sums[p][A];
			s += toString(std::max(0.0, (tmpGA - tmpAG)));
		}
	}
	return s + "]";
};

//---------------------------------------------------------------
//TPMDTables
//---------------------------------------------------------------
TPMDTables::TPMDTables(TReadGroups* ReadGroups, int maxLength){
	readGroups = ReadGroups;
	forward = new TPMDTable*[readGroups->numGroups];
	backward = new TPMDTable*[readGroups->numGroups];
	for(int i=0; i<readGroups->numGroups; ++i){
		forward[i] = new TPMDTable(maxLength);
		backward[i] = new TPMDTable(maxLength);
	}
};


void TPMDTables::add(TSite & site){
	for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
		forward[(*it)->readGroup]->add((*it)->posInRead, site.getRefBaseAsEnum(), (*it)->getBaseAsEnum());
		backward[(*it)->readGroup]->add((*it)->posInReadRev, site.getRefBaseAsEnum(), (*it)->getBaseAsEnum());
	}
};


void TPMDTables::writePMDFile(std::string filename){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "'!";

	//loop over all read groups
	for(int i=0; i<readGroups->numGroups; ++i){
		out << readGroups->getName(i) << "\t" << forward[i]->getPMDStringCT() << "\t" << backward[i]->getPMDStringGA() << "\n";
	}
	out.close();
}

void TPMDTables::writeTable(std::string filename){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "'!";

	//loop over all read groups
	for(int i=0; i<readGroups->numGroups; ++i){
		forward[i]->writeTable(out, readGroups->getName(i) + "\tforward\t");
		backward[i]->writeTable(out, readGroups->getName(i) + "\tbackward\t");
	}
	out.close();
};

void TPMDTables::writeTableWithCounts(std::string filename){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "'!";

	//loop over all read groups
	for(int i=0; i<readGroups->numGroups; ++i){
		forward[i]->writeTableWithCounts(out, readGroups->getName(i) + "\tforward\t");
		backward[i]->writeTableWithCounts(out, readGroups->getName(i) + "\tbackward\t");
	}
	out.close();
};
//---------------------------------------------------------------
//TPMDFunction
//---------------------------------------------------------------
TPMDSkoglund::TPMDSkoglund(double & Lambda, double & C){
	lambda = Lambda; c = C;
};
double TPMDSkoglund::getProb(int & pos){
	//Note: distance is zero based!
	return c + pow(1.0 - lambda, (double) pos) * lambda;
};
std::string TPMDSkoglund::getString(){
	return "P(pmd|pos) = p * (1 - p)^pos + c = " + toString(lambda) + " * (1 - " + toString(lambda) + ")^pos + " + toString(c);
};

//---------------------------------------------------------------
TPMDVeeramah::TPMDVeeramah(double & A, double & B, double & C){
	a = A; b = B; c = C;
}
double TPMDVeeramah::getProb(int & pos){
	//Note: distance is zero based!
	return a * exp(- (double) pos * b) + c;
};
std::string TPMDVeeramah::getString(){
	return "P(pmd|pos) = a * exp(- pos * b) + c = " + toString(a) + "* exp(- pos * " + toString(b) + ") + " + toString(c);
};

//---------------------------------------------------------------
TPMDEmpiric::TPMDEmpiric(std::string & values, std::string & example){
	std::vector<double> vec;
	fillVectorFromString(values, vec, ',');
	length = vec.size();
	if(length < 1) throw "Can not initialize post mortem damage function 'empiric[" + values + "]': wrong format!\n" + example;
	probs = new double[length];
	for(int i=0; i<length; ++i){
		probs[i] = vec[i];
	}
	last = probs[length-1];
};

double TPMDEmpiric::getProb(int & pos){
	if(pos < length) return probs[pos];
	else return last;
};

std::string TPMDEmpiric::getString(){
	std::string s = "P(pmd|pos) = {" + toString(probs[0]);
	for(int i=1; i<length; ++i) s += "," + toString(probs[i]);
	return  s +"}";
};

//------------------------------------------------------
//TPMD
//------------------------------------------------------
void TPMD::initializeFunction(std::string & pmdString, PMDType type){
	//parse string to get model.  options are
	// none
	// Skoglund[lambda,c]
	// Veeramah[a,b,c]
	std::string example = "Use either Skoglund[p,c], Veeramah[a,b,c] or Empiric[0.2,0.3,...]";

	//check if function was initialized abefore
	if(functionsInitialized[type]) throw "PMD function has been initialized previously!";

	//check if it is none
	if(pmdString == "none"){
		myFunctions[type] = new TPMDFunction();
	} else {
		std::string::size_type pos = pmdString.find_first_of('[');
		if(pos == std::string::npos) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
		std::string name = pmdString.substr(0,pos);

		//switch between functions
		if(name == "Empiric"){
			std::string::size_type endPos = pmdString.find_first_of(']');
			if(endPos == std::string::npos || endPos != pmdString.length()-1) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
			std::string list = pmdString.substr(pos+1, endPos-pos-1);
			myFunctions[type] = new TPMDEmpiric(list, example);
		} else {
			//prepare first value
			std::string tmp = pmdString.substr(pos+1, pmdString.length() - pos - 1);
			pos = tmp.find_first_of(',');
			if(pos == std::string::npos) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
			double first = atof(tmp.substr(0, pos).c_str());

			if(name == "Skoglund"){
				//get lambda and
				if(first < 0.0) throw "Can not initialize Skoglund function with lambda < 0!";
				if(first > 1.0) throw "Can not initialize Skoglund function with lambda > 1!";
				double c = atof(tmp.substr(pos+1).c_str());
				if(c < 0.0) throw "Can not initialize Skoglund function with c < 0!";
				myFunctions[type] = new TPMDSkoglund(first, c);
			} else if(name == "Veeramah"){
				//get a, b and c
				if(first < 0.0) throw "Can not initialize Veeramah function with a < 0!";

				//get b
				tmp = tmp.substr(pos+1);
				pos = tmp.find_first_of(',');
				if(pos == std::string::npos) throw "Can not initialize post mortem damage function '" + pmdString + "': wrong format!\n" + example;
				double b = atof(tmp.substr(0, pos).c_str());
				if(b < 0.0) throw "Can not initialize Veeramah function with b < 0!";

				//get c
				double c = atof(tmp.substr(pos+1).c_str());
				if(c < 0.0) throw "Can not initialize Veeramah function with c < 0!";

				//test if it can be > 1
				//if(first + c > 1) throw "Can not initialize Veeramah function with a + c > 1!";

				//initialze
				myFunctions[type] = new TPMDVeeramah(first, b, c);
			} else throw "Can not initialize post mortem damage function '" + pmdString + "': wrong name!\n" + example;
		}
	}
	functionsInitialized[type] = true;
};



