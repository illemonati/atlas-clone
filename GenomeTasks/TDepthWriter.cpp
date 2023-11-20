/*
 * TDepthWriter.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#include "TDepthWriter.h"

#include <ostream>
#include <vector>

#include "coretools/Main/TLog.h"
#include "TSite.h"
#include "TWindow.h"

namespace GenomeTasks{
using coretools::instances::logfile;

//----------------------------------------
// TDepthWriter
//----------------------------------------
void TDepthWriter::_handleWindow(){
	logfile().listFlush("Writing sequencing depth estimates to file ...");
	_out.writeNoDelim(_window.chrName(), ':', _window.from().position() + 1, '-', _window.to().position()).writeDelim();
	_out.writeln(_window.depth());
	logfile().done();

	logfile().listFlush("Adding per site depth to distribution ...");
	for(auto& s : _window){
		_distPerSite.add(s.depth());
	}
	logfile().done();
};

void TDepthWriter::writeDepth(){
	const std::string filename = _genome.outputName() + "_depthPerWindow.txt.gz";
	logfile().list("Writing per window depth estimates to '", filename, "'.");
	_out.open(filename, {"window", "depth"});

	_traverseBAMWindows();

	//write distribution
	logfile().list("Writing depth per site distribution to file '", _genome.outputName(), "_depthPerSiteHistogram.txt' ...");
	_distPerSite.write(_genome.outputName() + "_depthPerSiteHistogram.txt", "depth");
};


}; // end namespace


