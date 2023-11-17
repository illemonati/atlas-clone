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
#include "coretools/Main/TTask.h"
#include "coretools/Types/probability.h"

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
	TBamSample(const coretools::Probability & Prob, const std::string & OutName, BAM::TBamFile & bamFile);

	void close();

	void sample(BAM::TBamFile & bamfile);
	void downsampleRead(BAM::TAlignment & alignment);
	
friend class TBamDownsampler;
};

//-----------------------------------------
// TBamDownsampler_base
//-----------------------------------------
class TBamDownsampler_base:public TGenome_basic{
protected:
	std::vector<coretools::Probability> _probs;
	std::vector<std::string> _names;

	void _readVectorOfDownsamplingProbabilities();
};


//-----------------------------------------
// TBamDownsampler
//-----------------------------------------
class TBamDownsampler:public TBamDownsampler_base{
private:
	bool separateReads = false;
	bool _writeN        = false;
	std::vector<double> _cumulProbs;
	//lists to keep track of mates
	std::map<std::string, uint16_t> _mateWasWritten;
	BAM::TAlignmentList _discard;

protected:
	std::vector<TBamSample> _bamSamples;
public:
	TBamDownsampler();
	virtual void run();
	void sample();
};

//-----------------------------------------
// TBamSeparator
//-----------------------------------------
class TBamSeparator:public TBamDownsampler_base{
private:
	std::vector<double> _cumulProbs;
public:
	TBamSeparator();
	void run();
};

}; // end namespace

#endif /* TBAMDOWNSAMPLER_H_ */
