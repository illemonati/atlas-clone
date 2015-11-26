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

		//create genome object
		TGenome genome(&logfile, myParameters);

		//what to do?
		std::string task = myParameters.getParameterStringWithDefault("task", "estimate");
		if(task == "estimate"){
			logfile.startIndent("Running an EM algorithm to estimate heterozygosity:");
			genome.estimateTheta(myParameters);
		} else if(task == "LLsurface"){
			logfile.startIndent("Calculating the LL surface for each window:");
			genome.calcLikelihoodSurfaces(myParameters);
		} else if(task == "pileup"){
			logfile.startIndent("Printing pileup for each window:");
			genome.printPileup();
		//} else if(task=="calibration"){
		//	logfile.startIndent("Estimating eror calibration function:");
		//	genome.estimateErrorCalibration(myParameters);
		} else if(task == "calibration"){
			logfile.startIndent("Estimating error calibration function with EM:");
			genome.estimateErrorCalibrationEM(myParameters);
		} else if(task == "calibrationLLsurface"){
			logfile.startIndent("Estimating LL surface for error calibration function:");
			genome.calculateLikelihoodSurfaceErrorCalibrationEM(myParameters);
		} else if(task == "BQSR"){
			logfile.startIndent("Estimating recalibration parameters (BQSR):");
			genome.BQSR(myParameters);
		} else if(task == "callMLE"){
			logfile.startIndent("Calling MLE Genotypes:");
			genome.callMLEGenotypes(myParameters);
		} else if(task == "callBayes"){
			logfile.startIndent("Calling Bayesian Genotypes:");
			genome.callBayesianGenotypes(myParameters);
		} else if(task == "allelePresence"){
			logfile.startIndent("Calling Allele Presence:");
			genome.callAllelePresence(myParameters);
		} else if(task == "qualityTransformation"){
			logfile.startIndent("Printing Quality Transformation:");
			genome.printQualityTransformation(myParameters);
		} else if(task == "recalBAM"){
			logfile.startIndent("Recalibrating a BAM file:");
			genome.recalibrateBamFile(myParameters);
    	} else throw "Unknown task '" + task + "'!";
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
		logfile.error("unhandeled error!");
	}
	logfile.clearIndent();
	gettimeofday(&end, NULL);
	float runtime=(end.tv_sec  - start.tv_sec)/60.0;
	logfile.list("Program terminated in ", runtime, " min!");
	logfile.close();
    return 0;
}

