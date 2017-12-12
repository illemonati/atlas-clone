/*
 * TAtlasTest.h
 *
 *  Created on: Dec 11, 2017
 *      Author: phaentu
 */

#ifndef TATLASTEST_H_
#define TATLASTEST_H_

#include "TLog.h"
#include "TParameters.h"

//------------------------------------------
//TAtlasTest
//------------------------------------------
//Base class for individual tests.
//Tests can be combined to suites in TATlasTesting

class TAtlasTest{
protected:
	TLog* logfile;
	TParameters params;
	std::string _name;

public:
	TAtlasTest(TLog* Logfile){
		logfile = Logfile;
		_name = "empty";
	};
	virtual ~TAtlasTest(){};

	std::string name(){return _name;};

	virtual bool run(){
		return true;
	};
};


#endif /* TATLASTEST_H_ */
