/*
 * TDepthWriter.h
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#ifndef GENOMETASKS_TDEPTHWRITER_H_
#define GENOMETASKS_TDEPTHWRITER_H_

#include "TGenome.h"

//----------------------------------------
// TDepthWriter
//----------------------------------------
class TDepthWriter:public TGenome_windows{
private:
	TOutputFile _out;
	TCountDistribution _distPerSite;

	void _handleWindow();

public:
	TDepthWriter(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void writeDepth();
};





#endif /* GENOMETASKS_TDEPTHWRITER_H_ */
