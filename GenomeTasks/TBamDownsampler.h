/*
 * TBAmDownsampler.h
 *
 *  Created on: May 31, 2020
 *      Author: phaentu
 */

#ifndef TBAMDOWNSAMPLER_H_
#define TBAMDOWNSAMPLER_H_

#include "TGenome.h"

namespace GenomeTasks{

//-----------------------------------------
// TBamSample
//-----------------------------------------
class TBamSample{
private:
	double _prob;
	std::string _outName;
	BAM::TOutputBamFile _out;

	//lists to keep track of mates
	BAM::TAlignmentList _keep;
	BAM::TAlignmentList _discard;

public:
	TBamSample(const double Prob, const std	::string OutName);

	void open(BAM::TBamFile & bamFile);
	void close(TLog* logfile);

	void sample(BAM::TBamFile & bamfile, TRandomGenerator & randomGenerator);
	void downsampleRead(BAM::TAlignment & alignment, TRandomGenerator & randomGenerator, const TGenotypeMap & genoMap, const TQualityMap & qualMap);
};

//-----------------------------------------
// TBamDownsampler_base
//-----------------------------------------
class TBamDownsampler_base:public TGenome_basic{
private:
	std::vector<double> _probs;
	std::vector<std::string> _names;

	void _readVectorOfDownsamplingProbabilities(TParameters & Params);

public:
	TBamDownsampler_base(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
};


//-----------------------------------------
// TBamDownsampler
//-----------------------------------------
class TBamDownsampler:public TBamDownsampler_base{
protected:
	std::vector<TBamSample> _bamSamples;
public:
	TBamDownsampler(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	virtual ~TBamDownsampler(){};

	virtual void downsample();
};

//-----------------------------------------
// TBamReadDownsampler
//-----------------------------------------
class TBamReadDownsampler:TBamDownsampler{
private:

public:
	TBamReadDownsampler(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);

	void downsample();
};

//-----------------------------------------
// TBamSeparator
//-----------------------------------------
class TBamSeparator:public TBamDownsampler_base{
private:
	std::vector<double> _cumulProbs;

public:
	TBamSeparator(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);

	void separate();
};

}; // end namespace

#endif /* TBAMDOWNSAMPLER_H_ */
