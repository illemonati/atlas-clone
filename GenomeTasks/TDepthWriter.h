/*
 * TDepthWriter.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TDEPTHWRITER_H_
#define GENOMETASKS_TDEPTHWRITER_H_

#include "TGenome.h"
#include "TTask.h"

namespace GenomeTasks{

//----------------------------------------
// TDepthWriter
//----------------------------------------
class TDepthWriter:public TGenome_windows{
private:
	coretools::TOutputFile _out;
	coretools::TCountDistribution _distPerSite;

	void _handleWindow();

public:
	TDepthWriter(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void writeDepth();
};


//--------------------------------------
// Tasks
//--------------------------------------
class TTask_depthWriter:public coretools::TTask{
public:
	TTask_depthWriter(){ _explanation = "Estimating the distribution of depth among sites and writing depth per window"; };

	void run(){
		using namespace coretools::instances;
		TDepthWriter depthWriter(parameters(), &logfile(), &randomGenerator());
		depthWriter.writeDepth();
	};
};

}; // end namespace



#endif /* GENOMETASKS_TDEPTHWRITER_H_ */
