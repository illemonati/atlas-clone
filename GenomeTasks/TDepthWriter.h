/*
 * TDepthWriter.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TDEPTHWRITER_H_
#define GENOMETASKS_TDEPTHWRITER_H_

#include <string>

#include "TFile.h"
#include "TGenome.h"
#include "TTask.h"
#include "counters.h"

namespace GenomeTasks{

//----------------------------------------
// TDepthWriter
//----------------------------------------
class TDepthWriter:public TGenome_windows{
private:
	coretools::TOutputFile _out;
	coretools::TCountDistribution _distPerSite;

	void _handleWindow() override;
	void _handleAlignment() override {}
public:
	void writeDepth();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_depthWriter:public coretools::TTask{
public:
	TTask_depthWriter(){ _explanation = "Estimating the distribution of depth among sites and writing depth per window"; };

	void run(){
		TDepthWriter depthWriter;
		depthWriter.writeDepth();
	};
};

}; // end namespace



#endif /* GENOMETASKS_TDEPTHWRITER_H_ */
