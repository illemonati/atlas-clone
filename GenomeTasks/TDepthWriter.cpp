/*
 * TDepthWriter.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#include "TDepthWriter.h"

#include <ostream>
#include <vector>

#include "TLog.h"
#include "TSite.h"
#include "TWindow.h"

namespace GenomeTasks{
using coretools::instances::logfile;

//----------------------------------------
// TDepthWriter
//----------------------------------------
void TDepthWriter::_handleWindow(){
	logfile().listFlush("Writing sequencing depth estimates to file ...");
	_out.writeln(_window, _window.depth());
	logfile().done();

	logfile().listFlush("Adding per site depth to distribution ...");
	for(auto& s : _window){
		_distPerSite.add(s.depth());
	}
	logfile().done();
};

void TDepthWriter::writeDepth(){
	std::string filename = _outputName + "_depthPerWindow.txt.gz";
	logfile().list("Writing per window depth estimates to '" + filename + "'.");
	const std::vector<std::string> header = {"window", "depth"};
	_out.open(filename, header);

	_traverseBAMWindows();

	//write distribution
	filename = _outputName + "_depthPerSiteHistogram.txt";
	logfile().list("Writing depth per site distribution to file '" + filename + "' ...");
	_distPerSite.write(filename, "depth");
};


}; // end namespace


