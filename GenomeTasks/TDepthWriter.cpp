/*
 * TDepthWriter.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */





void TGenomeWindows::estimateApproximateDepthPerWindow(TParameters & params){
	//open output file
	std::ofstream output;
	std::string outputFileName = _outputName + "_depthPerWindow.txt";
	_logfile->list("Writing sequencing depth estimates to '" + outputFileName + "'");
	output.open(outputFileName.c_str());
	if(!output) throw "Failed to open output file '" + outputFileName + "'!";
	int nCharOnLine = 0;

	//write header
	output << "chr\tstart\tend\tdepth" << std::endl;

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			//write to file
			_logfile->listFlush("Writing sequencing depth estimates to file ...");
			if(window.depth == -1.0) output << alignmentParser.getCurChrName() << "\t" << window.start << "\t" << window.end << "\t" << "0" << "\n";
			else output << alignmentParser.getCurChrName() << "\t" << window.start << "\t" << window.end << "\t" << window.depth << "\n";
			_logfile->done();
		}
	}

	//clean up
	if(nCharOnLine > 0) output << '\n';
	output.close();
};

void TGenomeWindows::estimateDepthPerSite(TParameters & params){
	//initialize count object
	int maxDepth = params.getParameterIntWithDefault("maxDepth", 20);
	TDistributionOfCounts counts(maxDepth, "depth");

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			_logfile->listFlush("Adding depth to table ...");
			window.countDepthPerSite(counts);
			_logfile->done();
		}
	}

	//write
	std::string filename = _outputName + "_distDepthPerSite.txt";
	_logfile->listFlush("Writing depth distribution to '" + filename + "' ...");
	counts.writeCounts(filename);
	_logfile->done();

	filename = _outputName + "_cumulativeDepthPerSite.txt";
	_logfile->listFlush("Writing normalized cumulative depth distribution to '" + filename + "' ...");
	counts.writeNormalizedCumulativeCounts(filename);
	_logfile->done();

	filename = _outputName + "_quantilesDepthPerSite.txt";
	_logfile->listFlush("Writing quantiles of depth distribution '" + filename + "' ...");
	counts.writeQuantiles(filename);
	_logfile->done();
};

void TGenomeWindows::writeDepthPerSite(TParameters & params){
	gz::ogzstream out;

	std::string outputFileName = _outputName + "_depthPerSite.txt.gz";
	_logfile->list("Writing per site depth to '" + outputFileName + "'");

	out.open(outputFileName.c_str());
	if(!out) throw "Failed to open output file '" + outputFileName + "'!";

	//write header
	out << "chr\tpos\tdepth" << std::endl;

	//prepare windows
	TWindow window;

	//iterate through windows
	while(alignmentParser.readDataInNextWindow(window)){
		//write chromosome to file
		if(window.passedFilters){
			_logfile->listFlush("Writing depth per site ...");
			window.printDepthPerSite(out, alignmentParser.getCurChrName());
			_logfile->done();
		}
	}

	//clean up
	out.close();
};


