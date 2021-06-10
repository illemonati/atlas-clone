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
	TDepthWriter(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void writeDepth();
};


//--------------------------------------
// Tasks
//--------------------------------------
class TTask_depthWriter:public coretools::TTask{
public:
	TTask_depthWriter(){ _explanation = "Estimating the distribution of depth among sites and writing depth per window"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TDepthWriter depthWriter(Parameters, Logfile, _randomGenerator);
		depthWriter.writeDepth();
	};
};

}; // end namespace



#endif /* GENOMETASKS_TDEPTHWRITER_H_ */
