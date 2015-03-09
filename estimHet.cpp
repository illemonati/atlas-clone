/*
 * estimHet.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#include "TLog.h"
#include "TParameters.h"
#include "TGenome.h"

//---------------------------------------------------------------------------
//Main function
//---------------------------------------------------------------------------
int main(int argc, char* argv[]){
	struct timeval start, end;
    gettimeofday(&start, NULL);

	TLog logfile;
	logfile.newLine();
	logfile.write(" estimHet 1.0 ");
	logfile.write("**************");
    try{
		//read parameters from the command line
    	TParameters myParameters(argc, argv, &logfile);

		//verbose?
		bool verbose=myParameters.parameterExists("verbose");
		if(!verbose) logfile.listNoFile("Running in silent mode (use 'verbose' to get a status report on screen)");
		logfile.setVerbose(verbose);


		//open log file that handles the output
		std::string  logFilename=myParameters.getParameterString("logFile", false);
		if(logFilename.length()>0){
			logfile.openFile(logFilename.c_str());
			logfile.writeFileOnly(" estimHet 1.0 ");
			logfile.writeFileOnly("**************");
		}

		//switch among tasks
		std::string task=myParameters.getParameterString("task");
		if(task=="EM"){
			//run EM to estimate heterozygosity
			//read variables
			logfile.startIndent("Running an EM algorithm to estimate heterozygosity:");
			int maxIter = myParameters.getParameterIntWithDefault("maxIter", 100);
			logfile.list("Running the EM for up to " + toString(maxIter) + " iterations");
			double tol = myParameters.getParameterDoubleWithDefault("tol", 0.01);
			logfile.list("Stopping the EM if the LL has improved by less than " + toString(tol));
			double initHet = myParameters.getParameterDoubleWithDefault("init", 0.01);
			logfile.list("Starting the EM at H = " + toString(initHet));

			//read probabilities
			TGenome genome(&logfile);
			genome.readProbabilities(myParameters.getParameterString("probs"), myParameters.getParameterLongWithDefault("max", -1));

			//run EM
			genome.runEM(maxIter, tol, initHet);
		} else if(task=="test"){
			TGenome genome(&logfile, myParameters.getParameterString("bam"));
		}
		else throw "Unknown task!";

		//write unsused parameters
		std::string unusedParams=myParameters.getListOfUnusedParameters();
		if(unusedParams!="") logfile.warning("The following parameters were not used: " + unusedParams + "!");
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
		logfile.error("unhandeled error!");
	}
	logfile.clearIndent();
	gettimeofday(&end, NULL);
	float runtime=(end.tv_sec  - start.tv_sec)/60.0;
	logfile.list("Program terminated in ", runtime, " min!");
	logfile.close();
    return 0;
}

