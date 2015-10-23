/*
 * TPostMortemDamage.h
 *
 *  Created on: Oct 17, 2015
 *      Author: wegmannd
 */

#ifndef TPOSTMORTEMDAMAGE_H_
#define TPOSTMORTEMDAMAGE_H_

#include <math.h>

enum PMDType {pmdCT=0, pmdGA};
//---------------------------------------------------------------
//TPMD
//---------------------------------------------------------------
//Note: Base class is to be used when there is no PMD!
class TPMDFunction{
public:
	TPMDFunction(){};
	virtual ~TPMDFunction(){};
	virtual double getProb(double & pos){
		return 0.0;
	};
	virtual std::string getString(){ return "P(pmd|pos) = 0.0"; };
};

class TPMDSkoglund:public TPMDFunction{
private:
	double lambda, c;

public:
	TPMDSkoglund(double & Lambda, double & C){
		lambda = Lambda; c = C;
	};
	~TPMDSkoglund(){};
	double getProb(double & pos){
		return c + pow(1.0 - lambda, (double) pos - 1.0) * lambda;
	};
	std::string getString(){ return "P(pmd|pos) = p * (1 - p)^pos + c = " + toString(lambda) + " * (1 - " + toString(lambda) + ")^pos + " + toString(c); };
};

class TPMDVeeramah:public TPMDFunction{
private:
	double a,b,c;

public:
	TPMDVeeramah(double & A, double & B, double & C){
		a = A; b = B; c = C;
	}
	double getProb(double & pos){
		return a * exp(-pos * b) + c;
	};
	std::string getString(){ return "P(pmd|pos) = a * exp(- pos * b) + c = " + toString(a) + "* exp(- pos * " + toString(b) + ") + " + toString(c); };
};

class TPMD{
private:
	TPMDFunction* myFunctions[2];
	bool functionsInitialized[2];

public:
	TPMD(){
		myFunctions[pmdCT] = NULL;
		myFunctions[pmdGA] = NULL;
		functionsInitialized[pmdCT] = false;
		functionsInitialized[pmdGA] = false;
	};

	~TPMD(){
		if(functionsInitialized[pmdCT]) delete myFunctions[pmdCT];
		if(functionsInitialized[pmdGA]) delete myFunctions[pmdGA];
	};
	void initializeFunction(std::string & pmdString, PMDType type){
		//parse string to get model.  options are
		// none
		// Skoglund[lambda,c]
		// Veeramah[a,b,c]
		std::string example = "Use either Skoglund[p,c] or Veeramah[a,b,c]";

		//check if it is none
		if(pmdString == "none"){
			myFunctions[type] = new TPMDFunction();
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
				myFunctions[type] = new TPMDSkoglund(first, c);
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
				myFunctions[type] = new TPMDVeeramah(first, b, c);
			} else throw "Can not initialize post mortem damage function: wrong name!\n" + example;
		}
		functionsInitialized[type] = true;
	};
	double getProb(double pos, PMDType type){ return myFunctions[type]->getProb(pos); };
	double getProbCT(double pos){ return myFunctions[pmdCT]->getProb(pos); };
	double getProbGA(double pos){ return myFunctions[pmdGA]->getProb(pos); };
	std::string getFunctionString(PMDType type){ return myFunctions[type]->getString(); };
	bool functionInitialized(PMDType type){
		return functionsInitialized[type];
	};
};



#endif /* TPOSTMORTEMDAMAGE_H_ */
