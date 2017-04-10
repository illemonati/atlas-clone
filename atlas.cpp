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
	logfile.write(" ATLAS 1.0 ");
	logfile.write("***********");
    try{
		//read parameters from the command line
    	TParameters myParameters(argc, argv, &logfile);

		//verbose?
		bool verbose=myParameters.parameterExists("verbose");
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
			logfile.writeFileOnly(" ATLAS 1.0 ");
			logfile.writeFileOnly("***********");
		}

		//create genome object
		TGenome genome(&logfile, myParameters);

		//what to do?
		std::string task = myParameters.getParameterString("task");
		if(task == "estimateTheta"){
			logfile.startIndent("Running an EM algorithm to estimate heterozygosity (task = estimateTheta):");
			genome.estimateTheta(myParameters);
		} else if(task == "LLsurface"){
			logfile.startIndent("Calculating the LL surface for each window (task = LLSurface):");
			genome.calcLikelihoodSurfaces(myParameters);
		} else if(task == "pileup"){
			logfile.startIndent("Printing pileup for each window (task = pileup):");
			genome.printPileup(myParameters);
		} else if(task == "recal"){
			logfile.startIndent("Estimating error calibration function with EM (task = recal):");
			genome.estimateErrorCalibrationEM(myParameters);
		} else if(task == "recalLLsurface"){
			logfile.startIndent("Estimating LL surface for error calibration function (task = recalLLsurface):");
			genome.calculateLikelihoodSurfaceErrorCalibrationEM(myParameters);
		} else if(task == "BQSR"){
			logfile.startIndent("Estimating recalibration parameters (task = BQSR):");
			genome.BQSR(myParameters);
		} else if(task == "callMLE"){
			logfile.startIndent("Calling MLE Genotypes (task = callMLE):");
			genome.callMLEGenotypes(myParameters);
		} else if(task == "callBayes"){
			logfile.startIndent("Calling Bayesian Genotypes (task = callBayes):");
			genome.callBayesianGenotypes(myParameters);
		} else if(task == "allelePresence"){
			logfile.startIndent("Calling Allele Presence (task = allelePresence):");
			genome.callAllelePresence(myParameters);
		} else if(task == "combineBeagleFiles"){
			logfile.startIndent("combining beagle files (task = combineBeagleFiles):");
			genome.combineBeagleFiles(myParameters);
		} else if(task == "qualityTransformation"){
			logfile.startIndent("Printing Quality Transformation (task = qualityTransformation):");
			genome.printQualityTransformation(myParameters);
		} else if(task == "recalBAM"){
			logfile.startIndent("Recalibrating a BAM file (task = recalBAM):");
			genome.recalibrateBamFile(myParameters);
    	} else if(task == "splitRGbyLength"){
			logfile.startIndent("Splitting single end read groups in a BAM file (task = splitRGbyLength):");
			genome.splitSingleEndReadGroups(myParameters);
    	} else if(task == "mergeReadGroups"){
			logfile.startIndent("Merging read groups in a BAM file (task = mergeReadGroups):");
			genome.mergeReadGroups(myParameters);
    	} else if(task == "estimatePMD"){
    		logfile.startIndent("Estimating Post-Mortem Damage (PMD) patterns (task = estimatePMD):");
    		genome.estimatePMD(myParameters);
    	} else if(task == "PMDS"){
    		logfile.startIndent("Filtering for ancient reads using PMDS (Skoglund et al. 2014, task = PMDS):");
    		genome.runPMDS(myParameters);
    	} else if(task == "mergeReads"){
    		logfile.startIndent("Merging paired-end reads (task = mergeReads):");
    		genome.mergePairedEndReads(myParameters);
    	} else if(task == "PSMC"){
    	    logfile.startIndent("Generating a PSMC Input file probabilistically (task = PSMC):");
    	    genome.generatePSMCInput(myParameters);
    	} else if(task == "downsample"){
    		logfile.startIndent("Downsampling a bam file (task = downsample):");
    	    genome.downSampleBamFile(myParameters);
    	} else if(task == "downSampleReads"){
			logfile.startIndent("Downsampling a bam file (task = downSampleReads):");
			genome.downSampleReads(myParameters);
    	} else if(task == "BAMDiagnostics"){
    		logfile.startIndent("Estimating approximate coverage, read length frequencies and mapping quality frequencies (task = BAMDiagnostics):");
    	    genome.estimateApproximateCoverage(myParameters);
    	} else if(task == "coveragePerWindow"){
    		logfile.startIndent("Estimating coverage per window (task = coveragePerWindow):");
    	    genome.estimateApproximateCoveragePerWindow(myParameters);
    	} else if(task == "coveragePerSite"){
    		logfile.startIndent("Estimating coverage per site (task = coveragePerSitecoveragePerWindow:");
    		genome.estimateCoveragePerSite(myParameters);

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

