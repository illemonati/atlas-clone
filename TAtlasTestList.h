/*
 * TAtlasTestList.h
 *
 *  Created on: Dec 13, 2017
 *      Author: phaentu
 */

#ifndef TATLASTESTLIST_H_
#define TATLASTESTLIST_H_

#include "TAtlasTest.h"
#include "TParameters.h"


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
	TAtlasTestList(){
		//fill map of tests
		//NOTE: order will be the order in which test will be run if initialized from TParameters
		testMap.emplace_back("empty", &createInstance<TAtlasTest>);
		testMap.emplace_back("pileup", &createInstance<TAtlasTest_pileup>);
		testMap.emplace_back("BQSRSimulation", &createInstance<TAtlasTest_BQSRSimulation>);


		//fill map of test suites
		//Note: suites and tests within suites will be initialized in this order!
		testSuites["exampleSuite"] = {"empty", "pileup"};

		//automatically create a test suite "all"
		testSuites["all"] = {};
		for(testMapIt = testMap.begin(); testMapIt != testMap.end(); ++testMapIt)
			testSuites["all"].push_back(testMapIt->first);

	}

	~TAtlasTestList(){
		for(std::vector<TAtlasTest*>::iterator it=initializedTests.begin(); it!=initializedTests.end(); ++it)
			delete (*it);
	};

	bool initializeTest(const std::string & name, TParameters & params, TLog* logfile){
		//check if test exists
		if(testInitialized(name))
			return true;

		//if not create it
		for(testMapIt = testMap.begin(); testMapIt != testMap.end(); ++testMapIt){
			if(testMapIt->first == name){
				initializedTests.push_back( testMapIt->second(params, logfile) );
				return true;
			}
		}

		//only reached if test does not exist
		return false;
	};

	void parseTests(TParameters & params, TLog* logfile){
		for(testMapIt = testMap.begin(); testMapIt != testMap.end(); ++testMapIt){
			if(params.parameterExists(testMapIt->first))
				initializeTest(testMapIt->first, params, logfile);
		}
	};

	bool parseSuites(TParameters & params, TLog* logfile){
		bool returnVal = false;
		for(std::map<std::string, std::vector<std::string> >::iterator it = testSuites.begin(); it != testSuites.end(); ++it){
			if(params.parameterExists(it->first)){
				for(size_t i=0; i<it->second.size(); ++i)
					returnVal &= initializeTest(it->second[i], params, logfile);
			}
		}
		return returnVal;
	};

	size_t size(){
		return initializedTests.size();
	};

	bool testExists(const std::string name){
		for(testMapIt = testMap.begin(); testMapIt != testMap.end(); ++testMapIt){
			if(testMapIt->first == name)
				return true;
		}
		return false;
	};

	bool testInitialized(const std::string name){
		for(testIt=initializedTests.begin(); testIt!=initializedTests.end(); ++testIt){
			if((*testIt)->name() == name)
				return true;
		}
		return false;
	};

	void printTestToLogfile(TLog* logfile){
		for(testIt=initializedTests.begin(); testIt!=initializedTests.end(); ++testIt){
			logfile->list((*testIt)->name());
		}
	};

	TAtlasTest* operator[](size_t num){
		if(num < initializedTests.size())
			return initializedTests[num];
		else throw "Test number out of range!";
	};
};


#endif /* TATLASTESTLIST_H_ */
