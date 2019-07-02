/*
 * TTaskList.h
 *
 *  Created on: Mar 31, 2019
 *      Author: phaentu
 */

#ifndef TTASKLIST_H_
#define TTASKLIST_H_

#include <algorithm>
#include "TParameters.h"
#include "TLog.h"

//---------------------------------------------------------------------------
// Include list of tasks specific to the application
//---------------------------------------------------------------------------
#include "allTasks.h"

//---------------------------------------------------------------------------
// TTaskList
//---------------------------------------------------------------------------
class TTaskList{
private:
	std::map< std::string, TTask* > taskMap_all;
	std::map< std::string, TTask* > taskMap_regular;
	std::map< std::string, TTask* > taskMap_debug;

    void fillTaskMap_all();
    TTask* getTask(const std::string task, TLog* logfile);
    void throwErrorUnknownTask(const std::string & Task);

public:
    TTaskList();

    ~TTaskList();

    bool taskExists(const std::string task);
    void run(TParameters & parameters, TLog* logfile);
    void run(const std::string taskName, TParameters & parameters, TLog* logfile);
    void printAvailableTasks(TLog* logfile);
};



#endif /* TTASKLIST_H_ */
