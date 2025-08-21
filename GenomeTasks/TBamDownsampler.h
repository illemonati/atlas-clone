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

#include "TGenome.h"
#include "TOutputBamFile.h"
#include "coretools/Types/probability.h"

namespace BAM {
class TAlignment;
class TBamFile;
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

	void sample(BAM::TBamFile & bamfile);
	void downsampleRead(BAM::TAlignment & alignment);
	
friend class TBamDownsampler;
};

//-----------------------------------------
// TBamDownsampler
//-----------------------------------------
class TBamDownsampler {
	TGenome _genome{false};
	bool separateReads = false;
	bool _writeN       = false;
	std::vector<double> _cumulProbs;
	//lists to keep track of mates
	std::map<std::string, uint16_t> _mateWasWritten;
	BAM::TAlignmentList _discard;
	std::vector<coretools::Probability> _probs;
	std::vector<std::string> _names;
	std::vector<TBamSample> _bamSamples;

	void _readVectorOfDownsamplingProbabilities();

public:
	TBamDownsampler();
	void run();
	void sample();
};

}; // end namespace

#endif /* TBAMDOWNSAMPLER_H_ */
