/*
 * TAtlasTesting.cpp
 *
 *  Created on: Dec 11, 2017
 *      Author: phaentu
 */



#include "TAtlasTesting.h"



TAtlasTesting::TAtlasTesting(TParameters & params, TLog* Logfile){
	logfile = Logfile;

	//get name of report file
	outputName = params.getParameterStringWithDefault("out", "ATLAS_testing_report.txt");
	logfile->list("Writing testing report to '" + outputName + "'.");

	//add tests
	testList.parseSuites(params, logfile);
	testList.parseTests(params, logfile);
	printTests();
};

void TAtlasTesting::printTests(){
	if(testList.size() > 1)
		logfile->startIndent("Will run the following " + toString(testList.size()) + " tests:");
	else if(testList.size() == 1)
		logfile->startIndent("Will run the following test:");
	else throw "No tests requested!";

	testList.printTestToLogfile(logfile);
	logfile->endIndent();
};

void TAtlasTesting::addTest(std::string & name, TParameters & params){
	if(!testList.initializeTest(name, params, logfile))
		throw "Failed to initialize test '" + name + "': test does not exist!";
};

void TAtlasTesting::runTests(){
	//open report file
	std::ofstream out(outputName.c_str());
	if(!out)
		throw "Failed to open file '" + outputName + "' for writing!";

	//prepare test runs
	int numSuccess = 0;
	bool success;
	if(testList.size() < 1)
		throw "No tests requested!";

	//now run all tests
	if(testList.size() > 1)
		logfile->startNumbering("Running " + toString(testList.size()) + " tests:");
	else
		logfile->startNumbering("Running 1 test:");

	for(size_t testNum=0; testNum<testList.size(); ++testNum){
		//report test number and name
		logfile->numberWithIndent("Running test '" + testList[testNum]->name() + "' (test " + toString(testNum+1) + " of " + toString(testList.size()) + "):");
		out << testNum << '\t' << testList[testNum]->name() << '\t';

		//run test
		success = testList[testNum]->run();
		numSuccess += success;

		//report
		logfile->removeIndent();
		if(success){
			logfile->conclude("Test '" +testList[testNum]->name()  + "' passed!");
			out << "passed\n";
		} else {
			logfile->conclude("Test '" +testList[testNum]->name()  + "' failed!");
			out << "failed\n";
		}
	}
	logfile->endNumbering();

	//close report file
	out.close();
};


