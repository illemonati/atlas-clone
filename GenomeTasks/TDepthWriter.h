/*
 * TDepthWriter.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TDEPTHWRITER_H_
#define GENOMETASKS_TDEPTHWRITER_H_

#include <string>

#include "coretools/Files/TOutputFile.h"
#include "TGenome_OLD.h"
#include "coretools/Main/TTask.h"
#include "coretools/Math/counters.h"

namespace GenomeTasks{

//----------------------------------------
// TDepthWriter
//----------------------------------------
class TDepthWriter:public old::TGenome_windows{
private:
	coretools::TOutputFile _out;
	coretools::TCountDistribution<> _distPerSite;

	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _handleAlignment(BAM::TAlignment&) override {}
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
