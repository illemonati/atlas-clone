/*
 * estimHet.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#include "Tests/TAtlasTesting.h"
#include "TVcfDiagnostics.h"
#include "gitversion.h"

#include "TTaskList.h"

//---------------------------------------------------------------------------
//Main function
//---------------------------------------------------------------------------
int main(int argc, char* argv[]){
	struct timeval start, end;
    gettimeofday(&start, NULL);

	TLog logfile;
	logfile.newLine();
	logfile.flush(" ATLAS 1.0, commit ");
	std::string version = GITVERSION;
	logfile.write(version.substr(0,7));
	logfile.write("***************************");
    try{
		//read parameters from the command line
    	TParameters parameters(argc, argv, &logfile);

		//verbose?
		bool silent = parameters.parameterExists("silent");
		if(silent) logfile.listNoFile("Running in silent mode (omit argument 'silent' to get a status report on screen)");
		logfile.setVerbose(!silent);

		//warnings?
		bool suppressWarnings=parameters.parameterExists("suppressWarnings");
		if(suppressWarnings){
			logfile.list("Suppressing Warnings");
			logfile.suppressWarings();
		}

		//open log file that handles the output
		std::string  logFilename=parameters.getParameterString("logFile", false);
		if(logFilename.length()>0){
			logfile.openFile(logFilename.c_str());
			logfile.flush(" ATLAS 0.9, commit ");
			logfile.write(version.substr(0,7));
			logfile.writeFileOnly("***************************");
		}

		//is task provided?
		if(!parameters.parameterExists("task")){
			TTaskList taskList;
			logfile.setVerbose(true);
			taskList.printAvailableTasks(&logfile);
			throw "The parameter 'task' is not defined!";
		}

		//what to do?
		std::string task = parameters.getParameterString("task");
		if(task == "test"){
			//run testing utilities
			logfile.startIndent("Unit testing of atlas functionalities:");
			TAtlasTesting test(parameters, &logfile);
			test.runTests();
		} else {
			//run requested task
			TTaskList taskList;
			taskList.run(task, parameters, &logfile);
		}
		logfile.clearIndent();

		//write unsused parameters
		std::string unusedParams=parameters.getListOfUnusedParameters();
		if(unusedParams!=""){
			logfile.newLine();
			logfile.warning("The following parameters were not used: " + unusedParams + "!");
		}
    }
	catch (std::string & error){
		logfile.error(error);
	}
	catch (const char* error){
		logfile.error(error);
	}
	catch(std::exception & error){
		logfile.error(error.what());
	}
	catch (...){
		logfile.error("unhandled error!");
	}
	logfile.clearIndent();
	gettimeofday(&end, NULL);
	float runtime=(end.tv_sec  - start.tv_sec)/60.0;
	logfile.list("Program terminated in ", runtime, " min!");
	logfile.close();

    return 0;
}

