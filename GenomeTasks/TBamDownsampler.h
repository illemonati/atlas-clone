/*
 * TBAmDownsampler.h
 *
 *  Created on: May 31, 2020
 *      Author: phaentu
 */

#ifndef TBAMDOWNSAMPLER_H_
#define TBAMDOWNSAMPLER_H_

#include "TGenome.h"
#include "TTask.h"

namespace GenomeTasks{

using coretools::Probability;
using coretools::TTask;
using coretools::TParameters;
using coretools::TLog;
using coretools::TRandomGenerator;

//-----------------------------------------
// TBamSample
//-----------------------------------------
class TBamSample{
private:
	Probability _prob;
	std::string _outName;
	BAM::TOutputBamFile _out;

	//lists to keep track of mates
	BAM::TAlignmentList _keep;
	BAM::TAlignmentList _discard;

public:
	TBamSample(const Probability & Prob, const std::string & OutName);

	void open(BAM::TBamFile & bamFile);
	void close(TLog* logfile);

	void sample(BAM::TBamFile & bamfile, TRandomGenerator & randomGenerator);
	void downsampleRead(BAM::TAlignment & alignment, TRandomGenerator & randomGenerator);
};

//-----------------------------------------
// TBamDownsampler_base
//-----------------------------------------
class TBamDownsampler_base:public TGenome_basic{
protected:
	std::vector<Probability> _probs;
	std::vector<std::string> _names;

	void _readVectorOfDownsamplingProbabilities(TParameters & Parameters);

public:
	TBamDownsampler_base(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
};


//-----------------------------------------
// TBamDownsampler
//-----------------------------------------
class TBamDownsampler:public TBamDownsampler_base{
protected:
	std::vector<TBamSample> _bamSamples;
public:
	TBamDownsampler(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TBamDownsampler(){};

	virtual void downsample();
};

//-----------------------------------------
// TBamReadDownsampler
//-----------------------------------------
class TBamReadDownsampler:public TBamDownsampler{
private:

public:
	TBamReadDownsampler(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);

	void downsample();
};

//-----------------------------------------
// TBamSeparator
//-----------------------------------------
class TBamSeparator:public TBamDownsampler_base{
private:
	std::vector<double> _cumulProbs;

public:
	TBamSeparator(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator);

	void separate();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_downsample:public TTask{
public:
	TTask_downsample(){ _explanation = "Downsampling a BAM file by removing reads"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TBamDownsampler downsampler(Parameters, Logfile, _randomGenerator);
		downsampler.downsample();
	};
};

class TTask_downSampleReads:public TTask{
public:
	TTask_downSampleReads(){ _explanation = "Downsampling a BAM file by setting bases to N"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TBamReadDownsampler downsampler(Parameters, Logfile, _randomGenerator);
		downsampler.downsample();
	};
};

class TTask_separateReads:public TTask{
public:
	TTask_separateReads(){ _explanation = "Separating reads into different BAM files"; };

	void run(TParameters & Parameters, TLog* Logfile){
		TBamSeparator separator(Parameters, Logfile, _randomGenerator);
		separator.separate();
	};
};




}; // end namespace

#endif /* TBAMDOWNSAMPLER_H_ */
