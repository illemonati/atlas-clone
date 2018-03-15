/*
 * atlasTaskSwitcher.h
 *
 *  Created on: Dec 17, 2017
 *      Author: phaentu
 */

#ifndef ATLASTASKSWITCHER_H_
#define ATLASTASKSWITCHER_H_

#include "commonutilities/TLog.h"
#include "commonutilities/TParameters.h"
#include "TGenome.h"
#include "TDistanceEstimator.h"
#include "runSimulations.h"

//---------------------------------------------------------------------------
//Switch task
//---------------------------------------------------------------------------
class atlasTaskSwitcher{
private:
	 TParameters* parameters;
	 TLog* logfile;

public:
	atlasTaskSwitcher(TParameters* Parameters, TLog* Logfile){
		parameters = Parameters;
		logfile = Logfile;
	};

	void runTask(std::string task){
		//first all task that do not require TGenome
		if(task == "simulate"){
			logfile->startIndent("Generating simulations (task = simulate):");
			runSimulations(*parameters, logfile);
		} else if(task == "printGLF"){
			logfile->startIndent("Printing a GLF file to screen (task=printGLF):");
			TDistanceEstimator distEst(logfile);
			distEst.printGLF(*parameters);
		} else if(task == "estimateDist"){
			logfile->startIndent("Estimating the genetic distance between individuals (task=estimateDist):");
			TDistanceEstimator distEst(logfile);
			distEst.estimateDistances(*parameters);
		} else {
			//now all task that DO require TGenome
			TGenome genome(logfile, *parameters);
			if(task == "estimateTheta"){
				logfile->startIndent("Running an EM algorithm to estimate heterozygosity (task = estimateTheta):");
				genome.estimateTheta(*parameters);
			} else if(task == "LLsurface"){
				logfile->startIndent("Calculating the LL surface for each window (task = LLSurface):");
				genome.calcLikelihoodSurfaces(*parameters);
			} else if(task == "pileup"){
				logfile->startIndent("Printing pileup for each window (task = pileup):");
				genome.printPileup(*parameters);
			} else if(task == "recal"){
				logfile->startIndent("Estimating error calibration function with EM (task = recal):");
				genome.estimateErrorCalibrationEM(*parameters);
			} else if(task == "recalLL"){
				logfile->startIndent("Calculating LL for error calibration function (task = recalLL):");
				genome.calculateLikelihoodErrorCalibrationEM(*parameters);
			} else if(task == "BQSR"){
				logfile->startIndent("Estimating recalibration parameters (task = BQSR):");
				genome.BQSR(*parameters);
			} else if(task == "callMLE"){
				logfile->startIndent("Calling MLE Genotypes (task = callMLE):");
				genome.callMLEGenotypes(*parameters);
			} else if(task == "callBayes"){
				logfile->startIndent("Calling Bayesian Genotypes (task = callBayes):");
				genome.callBayesianGenotypes(*parameters);
			} else if(task == "allelePresence"){
				logfile->startIndent("Calling Allele Presence (task = allelePresence):");
				genome.callAllelePresence(*parameters);
			} else if(task == "randomBaseCaller"){
				logfile->startIndent("Calling random bases (task = randomBaseCaller");
				genome.randomBaseCaller(*parameters);
			} else if(task == "majorityBaseCaller"){
				logfile->startIndent("Calling random bases (task = majorityBaseCaller");
				genome.majorityBaseCaller(*parameters);
			} else if(task == "glf"){
				logfile->startIndent("Writing genotype likelihoods to a GLF file (task = glf):");
				genome.writeGLF(*parameters);
			} else if(task == "combineBeagleFiles"){
				logfile->startIndent("combining beagle files (task = combineBeagleFiles):");
				genome.combineBeagleFiles(*parameters);
			} else if(task == "qualityDist"){
				logfile->startIndent("Printing Quality Distribution (task = qualityDist):");
				genome.printQualityDistribution(*parameters);
			} else if(task == "qualityTransformation"){
				logfile->startIndent("Printing Quality Transformation (task = qualityTransformation):");
				genome.printQualityTransformation(*parameters);
			} else if(task == "recalBAM"){
				logfile->startIndent("Recalibrating a BAM file (task = recalBAM):");
				genome.recalibrateBamFile(*parameters);
			} else if(task == "binQualityScores"){
				logfile->startIndent("Binning quality scores (task = binQualityScores");
				genome.binQualityScores(*parameters);
			} else if(task == "assessSoftClipping"){
				logfile->startIndent("Assessing level of soft clipping in BAM file:");
				genome.assessSoftClipping(*parameters);
			} else if(task == "splitRGbyLength"){
				logfile->startIndent("Splitting single end read groups in a BAM file (task = splitRGbyLength):");
				genome.splitSingleEndReadGroups(*parameters);
			} else if(task == "mergeReadGroups"){
				logfile->startIndent("Merging read groups in a BAM file (task = mergeReadGroups):");
				genome.mergeReadGroups(*parameters);
			} else if(task == "estimatePMD"){
				logfile->startIndent("Estimating Post-Mortem Damage (PMD) patterns (task = estimatePMD):");
				genome.estimatePMD(*parameters);
			} else if(task == "PMDS"){
				logfile->startIndent("Filtering for ancient reads using PMDS (Skoglund et al. 2014, task = PMDS):");
				genome.runPMDS(*parameters);
			} else if(task == "mergeReads"){
				logfile->startIndent("Merging paired-end reads (task = mergeReads):");
				genome.mergePairedEndReads(*parameters);
			} else if(task == "PSMC"){
				logfile->startIndent("Generating a PSMC Input file probabilistically (task = PSMC):");
				genome.generatePSMCInput(*parameters);
			} else if(task == "downsample"){
				logfile->startIndent("Downsampling a bam file (task = downsample):");
				genome.downSampleBamFile(*parameters);
			} else if(task == "downSampleReads"){
				logfile->startIndent("Downsampling a bam file (task = downSampleReads):");
				genome.downSampleReads(*parameters);
			} else if(task == "BAMDiagnostics"){
				logfile->startIndent("Estimating approximate depth, read length frequencies and mapping quality frequencies (task = BAMDiagnostics):");
				genome.diagnoseBamFile(*parameters);
			} else if(task == "depthPerWindow"){
				logfile->startIndent("Estimating depth per window (task = depthPerWindow):");
				genome.estimateApproximateDepthPerWindow(*parameters);
			} else if(task == "depthPerSite"){
				logfile->startIndent("Estimating depth per site (task = depthPerSite):");
				genome.estimateDepthPerSite(*parameters);
			} else if(task == "writeDepthPerSite"){
				logfile->startIndent("Writing depth per site (task = writeDepthPerSite):");
				genome.writeDepthPerSite(*parameters);
			} else if(task=="createDepthMask"){
				logfile->startIndent("Creating depth mask BED file (task=createDepthMask:");
				genome.createDepthMask(*parameters);
			} else throw "Unknown task '" + task + "'!";
		}
		logfile->endIndent();
	};

	void runTask(){
		std::string task = parameters->getParameterString("task");
		runTask(task);
	};
};


#endif /* ATLASTASKSWITCHER_H_ */
