/*
 * TAtlasTestList.h
 *
 *  Created on: Dec 13, 2017
 *      Author: phaentu
 */

#ifndef TATLASTESTLIST_H_
#define TATLASTESTLIST_H_

#include "TAtlasTest.h"
#include "TAtlasTestRecalibration.h"
#include "TAtlasTestPMD.h"
#include "TParameters.h"
#include "TVcfTest.h"
#include "TAtlasTestMergePairs.h"
#include "TAtlasTestFilter.h"
#include "TAtlasTestLibraries.h"
#include "TAtlasTestLibrariesSpeed.h"


//------------------------------------------
//TAtlasTestList
//------------------------------------------
//class holding a map of names to constructors

//function to create an instance of TAtlasTest.
//Used to create map linking names to constructors.
template<typename T> TAtlasTest* createInstance(TParameters & params, TLog* logfile){return new T(params, logfile);};


class TAtlasTestList{
private:
	//maps test names to constuctor functions
	std::vector< std::pair<std::string, TAtlasTest*(*)(TParameters &, TLog*)> > testMap;
	std::vector< std::pair<std::string, TAtlasTest*(*)(TParameters &, TLog*)> >::iterator testMapIt;

	//maps test suit to test names
	std::map<std::string, std::vector<std::string> > testSuites;
	std::map<std::string, std::vector<std::string> >::iterator testSuiteIt;

	//vector of pointers to initialized tests
	std::vector<TAtlasTest*> initializedTests;
	std::vector<TAtlasTest*>::iterator testIt;


public:
    TAtlasTestList();

    ~TAtlasTestList();

    bool initializeTest(const std::string & name, TParameters & params, TLog* logfile);
    void parseTests(TParameters & params, TLog* logfile);
    bool parseSuites(TParameters & params, TLog* logfile);
    size_t size();
    bool testExists(const std::string name);
    bool testInitialized(const std::string name);
    void printTestToLogfile(TLog* logfile);
    TAtlasTest* operator[](size_t num);
};


#endif /* TATLASTESTLIST_H_ */
