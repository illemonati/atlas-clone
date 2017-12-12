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
	tests.emplace_back(logfile);
};

void TAtlasTesting::parseTests(TParameters & params){

}

void TAtlasTesting::addTest(std::string & name){
	//check if test already exists

};

void TAtlasTesting::runTests(){
	//open report file
	std::ofstream out(outputName.c_str());
	if(!out)
		throw "Failed to open file '" + outputName + "' for writing!";

	//prepare test runs
	int testNum = 1;
	int numSuccess = 0;
	bool success;
	if(tests.size() > 1)
		logfile->startNumbering("Running " + toString(tests.size()) + " tests:");
	else if(tests.size() == 1)
		logfile->startNumbering("Running " + toString(tests.size()) + " test:");
	else throw "No tests requested!";

	//now run all tests
	for(std::vector<TAtlasTest>::iterator t=tests.begin(); t!=tests.end(); ++t, ++testNum){
		//report test number and name
		logfile->startIndent("Running test '" + t->name() + "' (test " + toString(testNum) + " of " + toString(tests.size()) + "):");
		out << testNum << '\t' << t->name() << '\t';

		//run test
		success = t->run();
		numSuccess += t->run();

		//report
		logfile->removeIndent();
		if(success){
			logfile->conclude("Test passed!");
			out << "passed\n";
		} else {
			logfile->conclude("Test failed!");
			out << "failed\n";
		}
	}
	logfile->endIndent();

	//close report file
	out.close();
};


