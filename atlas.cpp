/*
 * estimHet.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#include "TAtlasTesting.h"
#include "atlasTaskSwitcher.h"
#include "gitversion.h"

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
    	TParameters myParameters(argc, argv, &logfile);

		//verbose?
		bool verbose = myParameters.parameterExists("verbose");
		if(!verbose) logfile.listNoFile("Running in silent mode (use 'verbose' to get a status report on screen)");
		logfile.setVerbose(verbose);
		//warnings?
		bool suppressWarnings=myParameters.parameterExists("suppressWarnings");
		if(suppressWarnings){
			logfile.list("Suppressing Warnings");
			logfile.suppressWarings();
		}

		//open log file that handles the output
		std::string  logFilename=myParameters.getParameterString("logFile", false);
		if(logFilename.length()>0){
			logfile.openFile(logFilename.c_str());
			logfile.flush(" ATLAS 1.0, commit ");
			logfile.write(version.substr(0,7));
			logfile.writeFileOnly("***************************");
		}

		//what to do?
		std::string task = myParameters.getParameterString("task");
		if(task == "test"){
			//run testing utilities
			logfile.startIndent("Unit testing of atlas functionalities:");
			TAtlasTesting test(myParameters, &logfile);
			test.runTests();
		} else {
			//run requested task
			atlasTaskSwitcher taskSwitcher(&myParameters, &logfile);
			taskSwitcher.runTask(task);
		}
		logfile.clearIndent();

		//write unsused parameters
		std::string unusedParams=myParameters.getListOfUnusedParameters();
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

