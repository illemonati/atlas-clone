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

    void initializeRandomGenerator(TParameters & parameters, TLog* logfile);

public:
    TTask();

    virtual ~TTask();

    virtual void run(TParameters & parameters, TLog* logfile);

    void run(std::string taskName, TParameters & parameters, TLog* logfile);
    std::string explanation();
    void printCitations(TLog* logfile);
};



#endif /* TTASK_H_ */
