/*
 * TAtlasTesting.h
 *
 *  Created on: Dec 11, 2017
 *      Author: phaentu
 */

#ifndef TATLASTESTING_H_
#define TATLASTESTING_H_

#include "TAtlasTestList.h"

//-----------------------------------
//TAtlasTesting
//-----------------------------------
class TTesting{
private:
	TLog* logfile;
	std::string outputName;
	TAtlasTestList testList;
	TParameters* params;

	void parseTests(TParameters & params);

public:
	TTesting(TParameters & params, TLog* Logfile);
	~TTesting(){};
	void printTests();
	void addTest(std::string & name, TParameters & params);
	void runTests();
};


#endif /* TATLASTESTING_H_ */
