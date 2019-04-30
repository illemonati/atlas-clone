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

//---------------------------------------------------------------------------
// TTask
//---------------------------------------------------------------------------
class TTask{
protected:
	std::string _explanation;
	std::vector<std::string> _citations;

public:
	TTask(){};
	virtual ~TTask(){};

	void run(std::string taskName, TParameters & parameters, TLog* logfile){
		logfile->startIndent(_explanation + " (task = " + taskName + "):");

		//print citations
		printCitations(logfile);

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
