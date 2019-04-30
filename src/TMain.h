/*
 * TMAin.h
 *
 *  Created on: Apr 1, 2019
 *      Author: phaentu
 */

#ifndef TMAIN_H_
#define TMAIN_H_

#include "TTaskList.h"
#include "gitversion.h"
#include "Tests/TTesting.h"

//---------------------------------------------------------------------------
// TMain
//---------------------------------------------------------------------------
class TMain{
private:
	struct timeval start;
	std::string applicationName;
	std::string version;
	std::string web;
	std::string title;
	TLog logfile;

	void fillTitleString(){
		title.clear();
		std::vector<std::string> tmp;

		//application name and version (plus underline and extra space)
		tmp.push_back(applicationName + " " + version);
		tmp.push_back(std::string(tmp[0].length() + 2, '-'));
		tmp.push_back("");

		//add lab
		tmp.push_back("Wegmann Lab - University of Fribourg, Switzerland");

		//add website
		if(!web.empty()) tmp.push_back(web);

		//commit (first extra space)
		tmp.push_back("");
		tmp.push_back(std::string("Commit ") + getGitVersion());

		//get max length
		size_t len = tmp[0].length();
		for(size_t i=1; i<tmp.size(); ++i){
			if(len < tmp[i].length())
				len = tmp[i].length();
		}
		len += 2; //add one space on either side

		//now compile title. Add spaces in front to align center
		for(size_t i=0; i<tmp.size(); ++i){
			if(tmp[i].empty()){
				title += '\n';
			} else {
				//add extra spaces
				std::string add((len - tmp[i].length()) / 2, ' ');
				title += add + tmp[i] + '\n';
			}
		}
	};

public:
	TMain(std::string ApplicationName, std::string Version, std::string Web){
		//start time conounting
		gettimeofday(&start, NULL);

		//store name
		applicationName = ApplicationName;
		version = Version;
		web = Web;

		//open logfile
		logfile.newLine();

		//compile and write title
		fillTitleString();
		logfile.write(title);
	};

	int run(int argc, char* argv[]){
		bool hadErrors = false;
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
				logfile.list("Suppressing Warnings.");
				logfile.suppressWarings();
			}

			//open log file, if requested
			std::string logFilename = parameters.getParameterString("logFile", false);
			if(logFilename.length() > 0){
				logfile.list("Will write log to file '" + logFilename + "'.");
				logfile.openFile(logFilename);
				logfile.writeFileOnly(title);
			}

			//create task list
			TTaskList taskList;

			//is a task provided?
			if(!parameters.parameterExists("task")){
				logfile.setVerbose(true);
				taskList.printAvailableTasks(&logfile);
				throw "The parameter 'task' is not defined!";
			}

			//what to do?
			std::string task = parameters.getParameterString("task");
			if(task == "test"){
				//run testing utilities
				logfile.startIndent("Testing of " + applicationName + " functionalities:");
				TTesting test(parameters, &logfile);
				test.runTests();
			} else {
				//run requested task
				taskList.run(task, parameters, &logfile);
			}
			logfile.clearIndent();

			//report not used arguments
			std::string unusedArguments = parameters.getListOfUnusedParameters();
			if(unusedArguments!=""){
				logfile.newLine();
				logfile.warning("The following arguments were not used: " + unusedArguments + "!");
			}
	    }
		catch (std::string & error){
			logfile.error(error);
			hadErrors = true;
		}
		catch (const char* error){
			logfile.error(error);
			hadErrors = true;
		}
		catch(std::exception & error){
			logfile.error(error.what());
			hadErrors = true;
		}
		catch (...){
			logfile.error("unhandled error!");
			hadErrors = true;
		}
		logfile.clearIndent();

		//calculate runtime
		struct timeval end;
		gettimeofday(&end, NULL);
		float runtime = (end.tv_sec  - start.tv_sec)/60.0;

		//report end of program and return exit status
		if(hadErrors){
			logfile.list(applicationName + " terminated with errors in ", runtime, " min!");
			logfile.close();
			return 1;
		} else {
			logfile.list(applicationName + " terminated successfully in ", runtime, " min!");
			logfile.close();
			return 0;
		}
	};
};


#endif /* TMAIN_H_ */
