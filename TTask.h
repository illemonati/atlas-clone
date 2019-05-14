/*
 * TTask.h
 *
 *  Created on: Mar 31, 2019
 *      Author: phaentu
 */

#ifndef TTASK_H_
#define TTASK_H_

#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"

//---------------------------------------------------------------------------
// TTask
//---------------------------------------------------------------------------
class TTask{
protected:
	std::string _explanation;
	std::vector<std::string> _citations;
	TRandomGenerator* randomGenerator;
	bool randomGeneratorInitialized;

	void initializeRandomGenerator(TParameters & parameters, TLog* logfile){
		if(!randomGeneratorInitialized){
			logfile->listFlush("Initializing random generator ...");

			if(parameters.parameterExists("fixedSeed")){
				randomGenerator=new TRandomGenerator(parameters.getParameterLong("fixedSeed"), true);
			} else if(parameters.parameterExists("addToSeed")){
				randomGenerator=new TRandomGenerator(parameters.getParameterLong("addToSeed"), false);
			} else {
				randomGenerator=new TRandomGenerator();
			}
			logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
			randomGeneratorInitialized = true;
		}
	};

public:
	TTask(){
		randomGenerator = nullptr;
		randomGeneratorInitialized = false;
	};
	virtual ~TTask(){
		if(randomGeneratorInitialized){
			delete randomGenerator;
		}
	};

	void run(std::string taskName, TParameters & parameters, TLog* logfile){
		logfile->startIndent(_explanation + " (task = " + taskName + "):");

		//print citations
		printCitations(logfile);

		initializeRandomGenerator(parameters, logfile);

		//now run task
		run(parameters, logfile);
		logfile->endIndent();
	};

	virtual void run(TParameters & parameters, TLog* logfile){
		throw "Base task can not be run!";
	};

	std::string explanation(){ return _explanation; };

	void printCitations(TLog* logfile){
		//write citations, if there are any
		if(!_citations.empty()){
			logfile->startIndent("Relevant citations:");
			for(std::vector<std::string>::iterator it = _citations.begin(); it != _citations.end(); ++it){
				logfile->list(*it);
			}
			logfile->endIndent();
		}
	};
};



#endif /* TTASK_H_ */
