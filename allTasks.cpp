/*
 * allTasks.cpp
 *
 *  Created on: Mar 31, 2019
 *      Author: phaentu
 */

#include "allTasks.h"

void fillTaskMaps(std::map< std::string, TTask* > & taskMap_regular, std::map< std::string, TTask* > & taskMap_debug){
	taskMap_regular.insert(std::pair<std::string, TTask*>("splitRGbyLength", new TTask_splitRGbyLength()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("mergeReadGroups", new TTask_mergeReadGroups()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("pileup", new TTask_pileup()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("BAMDiagnostics", new TTask_BAMDiagnostics()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("assessReadOverlap", new TTask_assessReadOverlap()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("mergeReads", new TTask_mergeReads()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("assessDuplication", new TTask_assessDuplication()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("mateInfo", new TTask_mateInfo()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("assessSoftClipping", new TTask_assessSoftClipping()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("removeSoftClippedBases", new TTask_removeSoftClippedBasesFromReads()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("writeDepthPerWindow", new TTask_writeDepthPerWindow()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("depthPerSiteDist", new TTask_depthPerSiteDist()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("writeDepthPerSite", new TTask_writeDepthPerSite()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("createDepthMask", new TTask_createDepthMask()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("allelicDepth", new TTask_allelicDepth()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("downsample", new TTask_downsample()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("downSampleReads", new TTask_downSampleReads()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("estimatePMD", new TTask_estimatePMD()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("PMDS", new TTask_PMDS()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("PSMC", new TTask_PSMC()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("recal", new TTask_recal()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("recalLL", new TTask_recalLL()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("BQSR", new TTask_BQSR()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("recalBAM", new TTask_recalBAM()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("qualityDist", new TTask_qualityDist()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("qualityTransformation", new TTask_qualityTransformation()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("binQualityScores", new TTask_binQualityScores()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("call", new TTask_call()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("estimateTheta", new TTask_estimateTheta()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("estimateThetaRatio", new TTask_estimateThetaRatio()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("GLF", new TTask_GLF()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("printGLF", new TTask_printGLF()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("majorMinor", new TTask_majorMinor()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("estimateDist", new TTask_estimateDist()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("estimateAlleleCounts", new TTask_estimateAlleleCounts()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("estimateAlleleFreq", new TTask_estimateAlleleFreq()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("estimateInbreeding", new TTask_estimateInbreeding()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("simulate", new TTask_simulate()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("VCFDiagnostics", new TTask_VCFDiagnostics()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("VCFToInvariantBed", new TTask_VCFToInvariantBed()));
	taskMap_regular.insert(std::pair<std::string, TTask*>("VCFFixInt", new TTask_VCFFixInt()));
	//taskMap_regular.insert(std::pair<std::string, TTask*>("filterVCF", new TTask_filterVCF()));

	//and debug tasks
	taskMap_debug.insert(std::pair<std::string, TTask*>("inbreedingLikelihood", new TTask_inbreedingLikelihood()));
	taskMap_debug.insert(std::pair<std::string, TTask*>("thetaLLSurface", new TTask_thetaLLSurface()));
};


