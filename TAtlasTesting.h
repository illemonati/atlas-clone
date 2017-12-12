/*
 * TAtlasTesting.h
 *
 *  Created on: Dec 11, 2017
 *      Author: phaentu
 */

#ifndef TATLASTESTING_H_
#define TATLASTESTING_H_

#include "TAtlasTest.h"
#include <vector>

//-----------------------------------
//TAtlasTesting
//-----------------------------------
class TAtlasTesting{
private:
	TLog* logfile;
	std::string outputName;

	std::vector<TAtlasTest> tests;

	void parseTests(TParameters & params);
	void addTest(std::string & name);

public:
	TAtlasTesting(TParameters & params, TLog* Logfile);

	void runTests();
};


#endif /* TATLASTESTING_H_ */
