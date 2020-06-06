/*
 * TDepthWriter.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#include "TDepthWriter.h"

namespace GenomeTasks{

//----------------------------------------
// TDepthWriter
//----------------------------------------
TDepthWriter::TDepthWriter(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Parameters, Logfile, RandomGenerator){

};

void TDepthWriter::_handleWindow(){
	_logfile->listFlush("Writing sequencing depth estimates to file ...");
	_window.writeCoordinates(_out);
	_out <<  _window.depth() << std::endl;
	_logfile->done();

	_logfile->listFlush("Adding per site depth to distribution ...");
	for(auto& s : _window){
		_distPerSite.add(s.depth());
	}
	_logfile->done();
};

void TDepthWriter::writeDepth(){
	std::string filename = _outputName + "_depthPerWindow.txt.gz";
	_logfile->list("Writing per window depth estimates to '" + filename + "'.");
	_out.open(filename, {"chr", "start", "end", "depth"});

	_traverseBAMWindows();

	//write distribution
	filename = _outputName + "_depthPerSiteHistogram.txt";
	_logfile->list("Writing depth per site distribution to file '" + filename + "' ...");
	_distPerSite.write(filename, "depth");
};


}; // end namespace


