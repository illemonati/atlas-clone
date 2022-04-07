/*
 * TBAmDownsampler.h
 *
 *  Created on: May 31, 2020
 *      Author: phaentu
 */

#ifndef TBAMDOWNSAMPLER_H_
#define TBAMDOWNSAMPLER_H_

#include <string>
#include <vector>

#include "TBamFile.h"
#include "TBamFilter.h"
#include "TGenome.h"
#include "TLog.h"
#include "TParameters.h"
#include "TRandomGenerator.h"
#include "TTask.h"
#include "probability.h"

namespace BAM {
class TAlignment;
}

namespace GenomeTasks{

//-----------------------------------------
// TBamSample
//-----------------------------------------
class TBamSample{
private:
	coretools::Probability _prob;
	std::string _outName;
	BAM::TOutputBamFile _out;

	//lists to keep track of mates
	BAM::TAlignmentList _keep;
	BAM::TAlignmentList _discard;

public:
	TBamSample(const coretools::Probability & Prob, const std::string & OutName);

	void open(BAM::TBamFile & bamFile);
	void close(coretools::TLog* logfile);

	void sample(BAM::TBamFile & bamfile, coretools::TRandomGenerator & randomGenerator);
	void downsampleRead(BAM::TAlignment & alignment, coretools::TRandomGenerator & randomGenerator);
};

//-----------------------------------------
// TBamDownsampler_base
//-----------------------------------------
class TBamDownsampler_base:public TGenome_basic{
protected:
	std::vector<coretools::Probability> _probs;
	std::vector<std::string> _names;

	void _readVectorOfDownsamplingProbabilities(coretools::TParameters & Parameters);

public:
	TBamDownsampler_base(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
};


//-----------------------------------------
// TBamDownsampler
//-----------------------------------------
class TBamDownsampler:public TBamDownsampler_base{
protected:
	std::vector<TBamSample> _bamSamples;
public:
	TBamDownsampler(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	virtual ~TBamDownsampler(){};

	virtual void downsample();
};

//-----------------------------------------
// TBamReadDownsampler
//-----------------------------------------
class TBamReadDownsampler:public TBamDownsampler{
private:

public:
	TBamReadDownsampler(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);

	void downsample();
};

//-----------------------------------------
// TBamSeparator
//-----------------------------------------
class TBamSeparator:public TBamDownsampler_base{
private:
	std::vector<double> _cumulProbs;

public:
	TBamSeparator(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);

	void separate();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_downsample:public coretools::TTask{
public:
	TTask_downsample(){ _explanation = "Downsampling a BAM file by removing reads"; };

	void run(){
		using namespace coretools::instances;
		TBamDownsampler downsampler(parameters(), &logfile(), &randomGenerator());
		downsampler.downsample();
	};
};

class TTask_downSampleReads:public coretools::TTask{
public:
	TTask_downSampleReads(){ _explanation = "Downsampling a BAM file by setting bases to N"; };

	void run(){
		using namespace coretools::instances;
		TBamReadDownsampler downsampler(parameters(), &logfile(), &randomGenerator());
		downsampler.downsample();
	};
};

class TTask_separateReads:public coretools::TTask{
public:
	TTask_separateReads(){ _explanation = "Separating reads into different BAM files"; };

	void run(){
		using namespace coretools::instances;
		TBamSeparator separator(parameters(), &logfile(), &randomGenerator());
		separator.separate();
	};
};




}; // end namespace

#endif /* TBAMDOWNSAMPLER_H_ */
